#include "./nvme.h"

static uint16_t nvme_io_cmd(FemuCtrl *n, NvmeCmd *cmd, NvmeRequest *req);

static void nvme_update_sq_eventidx(const NvmeSQueue *sq)
{
    if (sq->eventidx_addr_hva) {
        *((uint32_t *)(sq->eventidx_addr_hva)) = sq->tail;
        return;
    }

    if (sq->eventidx_addr) {
        nvme_addr_write(sq->ctrl, sq->eventidx_addr, (void *)&sq->tail,
                        sizeof(sq->tail));
    }
}

static inline void nvme_copy_cmd(NvmeCmd *dst, NvmeCmd *src)
{
#if defined(__AVX__)
    __m256i *d256 = (__m256i *)dst;
    const __m256i *s256 = (const __m256i *)src;

    _mm256_store_si256(&d256[0], _mm256_load_si256(&s256[0]));
    _mm256_store_si256(&d256[1], _mm256_load_si256(&s256[1]));
#elif defined(__SSE2__)
    __m128i *d128 = (__m128i *)dst;
    const __m128i *s128 = (const __m128i *)src;

    _mm_store_si128(&d128[0], _mm_load_si128(&s128[0]));
    _mm_store_si128(&d128[1], _mm_load_si128(&s128[1]));
    _mm_store_si128(&d128[2], _mm_load_si128(&s128[2]));
    _mm_store_si128(&d128[3], _mm_load_si128(&s128[3]));
#else
    *dst = *src;
#endif
}

static void nvme_process_sq_io(void *opaque, int index_poller)
{
    NvmeSQueue *sq = opaque;
    FemuCtrl *n = sq->ctrl;

    uint16_t status;
    hwaddr addr;
    NvmeCmd cmd;
    NvmeRequest *req;
    int processed = 0;

    nvme_update_sq_tail(sq);
    while (!(nvme_sq_empty(sq))) {
        if (sq->phys_contig) {
            addr = sq->dma_addr + sq->head * n->sqe_size;
            nvme_copy_cmd(&cmd, (void *)&(((NvmeCmd *)sq->dma_addr_hva)[sq->head]));
        } else {
            addr = nvme_discontig(sq->prp_list, sq->head, n->page_size,
                                  n->sqe_size);
            nvme_addr_read(n, addr, (void *)&cmd, sizeof(cmd));
        }
        nvme_inc_sq_head(sq);

        req = QTAILQ_FIRST(&sq->req_list);
        QTAILQ_REMOVE(&sq->req_list, req, entry);
        memset(&req->cqe, 0, sizeof(req->cqe));
        /* Coperd: record req->stime at earliest convenience */
        req->expire_time = req->stime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        req->cqe.cid = cmd.cid;
        req->cmd_opcode = cmd.opcode;
        memcpy(&req->cmd, &cmd, sizeof(NvmeCmd));

        if (n->print_log) {
            femu_debug("%s,cid:%d\n", __func__, cmd.cid);
        }

        status = nvme_io_cmd(n, &cmd, req);
        if (1 && status == NVME_SUCCESS) {
            req->status = status;

            int rc = femu_ring_enqueue(n->to_ftl[index_poller], (void *)&req, 1);
            if (rc != 1) {
                femu_err("enqueue failed, ret=%d\n", rc);
            }
        } else if (status == NVME_SUCCESS) {
            /* Normal I/Os that don't need delay emulation */
            req->status = status;
        } else {
            femu_err("Error IO processed!\n");
        }

        processed++;
    }

    nvme_update_sq_eventidx(sq);
    sq->completed += processed;
}

static void nvme_post_cqe(NvmeCQueue *cq, NvmeRequest *req)
{
    FemuCtrl *n = cq->ctrl;
    NvmeSQueue *sq = req->sq;
    NvmeCqe *cqe = &req->cqe;
    uint8_t phase = cq->phase;
    hwaddr addr;

    if (n->print_log) {
        femu_debug("%s,req,lba:%lu,lat:%lu\n", n->devname, req->slba, req->reqlat);
    }
    cqe->status = cpu_to_le16((req->status << 1) | phase);
    cqe->sq_id = cpu_to_le16(sq->sqid);
    cqe->sq_head = cpu_to_le16(sq->head);

    if (cq->phys_contig) {
        addr = cq->dma_addr + cq->tail * n->cqe_size;
        ((NvmeCqe *)cq->dma_addr_hva)[cq->tail] = *cqe;
    } else {
        addr = nvme_discontig(cq->prp_list, cq->tail, n->page_size, n->cqe_size);
        nvme_addr_write(n, addr, (void *)cqe, sizeof(*cqe));
    }

    nvme_inc_cq_tail(cq);
}

static void nvme_process_cq_cpl(void *arg, int index_poller)
{
    FemuCtrl *n = (FemuCtrl *)arg;
    NvmeCQueue *cq = NULL;
    NvmeRequest *req = NULL;
    struct rte_ring *rp = n->to_ftl[index_poller];
    pqueue_t *pq = n->pq[index_poller];
    uint64_t now;
    int processed = 0;
    int rc;
    int i;

    if (BBSSD(n)) {
        rp = n->to_poller[index_poller];
    }

    while (femu_ring_count(rp)) {
        req = NULL;
        rc = femu_ring_dequeue(rp, (void *)&req, 1);
        if (rc != 1) {
            femu_err("dequeue from to_poller request failed\n");
        }
        assert(req);

        pqueue_insert(pq, req);
    }

    while ((req = pqueue_peek(pq))) {
        now = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        if (now < req->expire_time) {
            break;
        }

        cq = n->cq[req->sq->sqid];
        if (!cq->is_active)
            continue;
        nvme_post_cqe(cq, req);
        QTAILQ_INSERT_TAIL(&req->sq->req_list, req, entry);
        pqueue_pop(pq);
        processed++;
        n->nr_tt_ios++;
        if (now - req->expire_time >= 20000) {
            n->nr_tt_late_ios++;
            if (n->print_log) {
                femu_debug("%s,diff,pq.count=%lu,%" PRId64 ", %lu/%lu\n",
                           n->devname, pqueue_size(pq), now - req->expire_time,
                           n->nr_tt_late_ios, n->nr_tt_ios);
            }
        }
        n->should_isr[req->sq->sqid] = true;
    }

    if (processed == 0)
        return;

    switch (n->multipoller_enabled) {
    case 1:
        nvme_isr_notify_io(n->cq[index_poller]);
        break;
    default:
        for (i = 1; i <= n->nr_io_queues; i++) {
            if (n->should_isr[i]) {
                nvme_isr_notify_io(n->cq[i]);
                n->should_isr[i] = false;
            }
        }
        break;
    }
}

void *nvme_poller(void *arg)
{
    FemuCtrl *n = ((NvmePollerThreadArgument *)arg)->n;
    int index = ((NvmePollerThreadArgument *)arg)->index;
    int i;

    switch (n->multipoller_enabled) {
    case 1:
        while (1) {
            if ((!n->dataplane_started)) {
                usleep(1000);
                continue;
            }

            NvmeSQueue *sq = n->sq[index];
            NvmeCQueue *cq = n->cq[index];
            if (sq && sq->is_active && cq && cq->is_active) {
                nvme_process_sq_io(sq, index);
            }
            nvme_process_cq_cpl(n, index);
        }
        break;
    default:
        while (1) {
            if ((!n->dataplane_started)) {
                usleep(1000);
                continue;
            }

            for (i = 1; i <= n->nr_io_queues; i++) {
                NvmeSQueue *sq = n->sq[i];
                NvmeCQueue *cq = n->cq[i];
                if (sq && sq->is_active && cq && cq->is_active) {
                    nvme_process_sq_io(sq, index);
                }
            }
            nvme_process_cq_cpl(n, index);
        }
        break;
    }

    return NULL;
}

static inline uint16_t nvme_make_pid(NvmeEnduranceGroup *endgrp, uint16_t rg, uint16_t ph) //update~
{
	uint16_t rgif = endgrp->fdp.rgif; 

	if (rgif == 0)
		return ph;

	return (rg << (16 - rgif) | ph);												
}																				//~update

static inline bool nvme_ph_valid(NvmeNamespace *ns, uint16_t ph)				//update~
{
	return ph < ns->fdp.nphs;
}																				//~update

static inline bool nvme_rg_valid(NvmeEnduranceGroup *endgrp, uint16_t rg)		//update~
{
	return rg < endgrp->fdp.nrg;
}																				//~update

static inline uint16_t nvme_pid2ph(NvmeNamespace *ns, uint16_t pid)				//update~
{
	uint16_t rgif = ns->endgrp->fdp.rgif; 

	if (rgif == 0)
		return pid;

	return pid & ((1 << (15 - rgif))- 1);
}																				//~update

static inline uint16_t nvme_pid2rg(NvmeNamespace *ns, uint16_t pid)				//update~
{
	uint16_t rgif = ns->endgrp->fdp.rgif;

	if (rgif == 0)
		return 0;

	return pid >> (16 - rgif);
}																				//~update 

static inline bool nvme_parse_pid(NvmeNamespace *ns, uint16_t pid,
							uint16_t *ph, uint16_t *rg)							//update~
{
	*ph = nvme_pid2ph(ns, pid);
	*rg = nvme_pid2rg(ns, pid); 

	return nvme_ph_valid(ns, *ph) && nvme_rg_valid(ns->endgrp, *rg);
}																				//~update 

static inline void nvme_fdp_stat_inc(uint64_t *a, uint64_t b)					//update~
{
    uint64_t ret = *a + b;
    *a = ret < *a ? UINT64_MAX : ret;
} 																				//~update

static uint16_t nvme_io_mgmt_recv_ruhs(FemuCtrl *n, NvmeNamespace* ns, NvmeCmd *cmd, //update~
											NvmeRequest *req, size_t len) 
{ 
	printf("nvme_io_mgmt_recv_ruhs() called\n");
    NvmeEnduranceGroup *endgrp;
    NvmeRuhStatus *hdr;
    NvmeRuhStatusDescr *ruhsd;
    unsigned int nruhsd;
    uint16_t rg, ph, *ruhid;
    size_t trans_len;
    uint8_t *buf = NULL;  // A buffer to be transmitted to host
	uint64_t prp1 = le64_to_cpu(cmd->dptr.prp1); //update
	uint64_t prp2 = le64_to_cpu(cmd->dptr.prp2); //update

	printf("here0\n");
    if (ns->id == 0 || ns->id == 0xffffffff) {
        return NVME_INVALID_NSID | NVME_DNR;
    }

	printf("here0\n");
    if (!ns->endgrp->fdp.enabled) {
        return NVME_FDP_DISABLED | NVME_DNR;
    }

	printf("here0\n");
    endgrp = ns->endgrp;

	printf("here1\n");
    nruhsd = ns->fdp.nphs * endgrp->fdp.nrg; // The number of streams in the endurance group
    trans_len = sizeof(NvmeRuhStatus) + nruhsd * sizeof(NvmeRuhStatusDescr);
    buf = g_malloc(trans_len);

	printf("here2\n");
    trans_len = MIN(trans_len, len);

	printf("here3\n");
    hdr = (NvmeRuhStatus *)buf; // Start Address of RUHS
    ruhsd = (NvmeRuhStatusDescr *)(buf + sizeof(NvmeRuhStatus)); // Start Address of RUHSD

	printf("here4\n");
	// header buffering
    hdr->nruhsd = cpu_to_le16(nruhsd); 

	printf("here5\n");
    ruhid = ns->fdp.phs;

	printf("here6\n");
	// ruhsd buffering
    for (ph = 0; ph < ns->fdp.nphs; ph++, ruhid++) {
        NvmeRuHandle *ruh = &endgrp->fdp.ruhs[*ruhid];

        for (rg = 0; rg < endgrp->fdp.nrg; rg++, ruhsd++) {
            uint16_t pid = nvme_make_pid(endgrp, rg, ph);

            ruhsd->pid = cpu_to_le16(pid); 
            ruhsd->ruhid = *ruhid;
            ruhsd->earutr = 0;
            ruhsd->ruamw = cpu_to_le64(ruh->rus[rg].ruamw);
        }
    }

	printf("here7\n");
    return dma_read_prp(n, buf, trans_len, prp1, prp2);
}																		//~update 

static uint16_t nvme_io_mgmt_recv(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,	
													NvmeRequest *req)		//udpate~	
{
    uint32_t cdw10 = le32_to_cpu(cmd->cdw10);
    uint32_t numd = le32_to_cpu(cmd->cdw11);
    uint8_t mo = (cdw10 & 0xff); // Management operation
    size_t len = (numd + 1) << 2; // unit: byte

    switch (mo) {
    case NVME_IOMR_MO_NOP:
        return 0;
    case NVME_IOMR_MO_RUH_STATUS:
        return nvme_io_mgmt_recv_ruhs(n, ns, cmd, req, len);
    default:
        return NVME_INVALID_FIELD | NVME_DNR;
    };
}																	//~update

static NvmeFdpEvent *nvme_fdp_alloc_event(FemuCtrl *n, NvmeFdpEventBuffer *ebuf)//update~
{
	NvmeFdpEvent *ret = NULL;
	bool is_full = ebuf->next == ebuf->start && ebuf->nelems;

	ret = &ebuf->events[ebuf->next++]; // allocation
	if (ebuf->next == NVME_FDP_MAX_EVENTS)
		ebuf->next = 0;
	if (is_full)
		ebuf->start = ebuf->next;
	else
		ebuf->nelems++;

	memset(ret, 0, sizeof(NvmeFdpEvent)); // fill the 'ret' with 0
	ret->timestamp = n->features.time_stamp; //update

	return ret;
}																				//~update

static inline int log_event(NvmeRuHandle *ruh, uint8_t event_type)				//update~
{
    return (ruh->event_filter >> nvme_fdp_evf_shifts[event_type]) & 0x1;
}																				//~update

static bool nvme_update_ruh(FemuCtrl *n, NvmeNamespace *ns, uint16_t pid)		//update~
{
	NvmeEnduranceGroup *endgrp = ns->endgrp;
	NvmeRuHandle *ruh;
    const uint8_t lba_index = NVME_ID_NS_FLBAS_INDEX(ns->id_ns.flbas);
    const uint8_t data_shift = ns->id_ns.lbaf[lba_index].lbads;
	NvmeReclaimUnit *ru;
	NvmeFdpEvent *e = NULL;
    uint64_t data_size;
	uint16_t ph, rg, ruhid;

	if (!nvme_parse_pid(ns, pid, &ph, &rg)) 
		return false;

	/* A stage that changes PHNDL to RUH */
	ruhid = ns->fdp.phs[ph];				
	ruh = &endgrp->fdp.ruhs[ruhid];	
	ru = &ruh->rus[rg]; 
	
    data_size = (uint64_t)ru->ruamw << data_shift;

	if (ru->ruamw) 
	{
	// Logging only when the device updates an RUH implicitly,
	// Not logging when the host updates an RUH explicitly by i/o mgmt cmd.
		if (log_event(ruh, FDP_EVT_RU_NOT_FULLY_WRITTEN))
		{
			e = nvme_fdp_alloc_event(n, &endgrp->fdp.host_events);
			e->type = FDP_EVT_RU_NOT_FULLY_WRITTEN;
			e->flags = FDPEF_PIV | FDPEF_NSIDV | FDPEF_LV;
			e->pid = cpu_to_le16(pid);
			e->nsid = cpu_to_le32(ns->id);
			e->rgid = cpu_to_le16(rg);
			e->ruhid = cpu_to_le16(ruhid);
		}
		
		nvme_fdp_stat_inc(&endgrp->fdp.mbmw, data_size);
	}

	ru->ruamw = ruh->ruamw; // Reset

	return true;
}	 																				//~update

static uint16_t nvme_io_mgmt_send_ruh_update(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd, //update~
															NvmeRequest *req)			
{
    uint32_t cdw10 = le32_to_cpu(cmd->cdw10);
    uint16_t ret = NVME_SUCCESS;
    uint32_t npid = (cdw10 >> 1) + 1;
    int i;
    uint16_t *pids = NULL;
    uint32_t maxnpid = n->endgrps->fdp.nrg * n->endgrps->fdp.nruh;
	uint64_t prp1 = le64_to_cpu(cmd->dptr.prp1);
	uint64_t prp2 = le64_to_cpu(cmd->dptr.prp2);

    if (unlikely(npid >= MIN(NVME_FDP_MAXPIDS, maxnpid))) {
        return NVME_INVALID_FIELD | NVME_DNR;
    }

    pids = g_new(uint16_t, npid);

	ret = dma_write_prp(n, (uint8_t *)pids, npid * sizeof(uint16_t), prp1, prp2);

    if (ret) {
        return ret;
    }

    for (i = 0; i < npid; i++) {
        if (!nvme_update_ruh(n, ns, pids[i])) {
            return NVME_INVALID_FIELD | NVME_DNR;
        }
    }

    return ret;
}																					//~update

static uint16_t nvme_io_mgmt_send(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd, //update~
												NvmeRequest *req)	
{
    uint32_t cdw10 = le32_to_cpu(cmd->cdw10);
    uint8_t mo = (cdw10 & 0xff);

    switch (mo) {
    case NVME_IOMS_MO_NOP:
        return 0;
    case NVME_IOMS_MO_RUH_UPDATE:
        return nvme_io_mgmt_send_ruh_update(n, ns, cmd, req);
    default:
        return NVME_INVALID_FIELD | NVME_DNR;
    };
}																	//~update 

static void nvme_do_write_fdp(FemuCtrl *n, NvmeRequest *req, uint32_t nlb) //update~
{
    NvmeNamespace *ns = req->ns;
    NvmeRwCmd *rw = (NvmeRwCmd *)&req->cmd;
    const uint8_t lba_index = NVME_ID_NS_FLBAS_INDEX(ns->id_ns.flbas);
    const uint8_t data_shift = ns->id_ns.lbaf[lba_index].lbads;
    uint64_t data_size = (uint64_t)nlb << data_shift;
    uint32_t dw12 = le32_to_cpu(req->cmd.cdw12);
    uint8_t dtype = (dw12 >> 20) & 0xf;
    uint16_t pid = le16_to_cpu(rw->dspec);
    uint16_t ph, rg, ruhid;
    NvmeReclaimUnit *ru;

	/* Normal write */
    if (dtype != NVME_DIRECTIVE_DATA_PLACEMENT || !nvme_parse_pid(ns, pid, &ph, &rg))
	{
        ph = 0;
        rg = 0;
    } 

    ruhid = ns->fdp.phs[ph];
    ru = &ns->endgrp->fdp.ruhs[ruhid].rus[rg];

    nvme_fdp_stat_inc(&ns->endgrp->fdp.hbmw, data_size);
    nvme_fdp_stat_inc(&ns->endgrp->fdp.mbmw, data_size);

    while (nlb) {
        if (nlb < ru->ruamw) {
            ru->ruamw -= nlb;
            break;
        }

        nlb -= ru->ruamw;
        nvme_update_ruh(n, ns, pid);
    }
} 																			//~update

uint16_t nvme_rw(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd, NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint16_t ctrl = le16_to_cpu(rw->control);
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint64_t prp1 = le64_to_cpu(rw->prp1);
    uint64_t prp2 = le64_to_cpu(rw->prp2);
    const uint8_t lba_index = NVME_ID_NS_FLBAS_INDEX(ns->id_ns.flbas);
    const uint16_t ms = le16_to_cpu(ns->id_ns.lbaf[lba_index].ms);
    const uint8_t data_shift = ns->id_ns.lbaf[lba_index].lbads;
    uint64_t data_size = (uint64_t)nlb << data_shift;
    uint64_t data_offset = slba << data_shift;
    uint64_t meta_size = nlb * ms;
    uint64_t elba = slba + nlb;
    uint16_t err;
    int ret;

    req->is_write = (rw->opcode == NVME_CMD_WRITE) ? 1 : 0;

    err = femu_nvme_rw_check_req(n, ns, cmd, req, slba, elba, nlb, ctrl,
                                 data_size, meta_size);
    if (err) {
        return err;
	}

	if (ns->endgrp && ns->endgrp->fdp.enabled) { //update~
		nvme_do_write_fdp(n, req, nlb);  		
	}											 //~update

    if (nvme_map_prp(&req->qsg, &req->iov, prp1, prp2, data_size, n)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, prp1), 0, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }

    assert((nlb << data_shift) == req->qsg.size);

    req->slba = slba;
    req->status = NVME_SUCCESS;
    req->nlb = nlb;

    ret = backend_rw(n->mbe, &req->qsg, &data_offset, req->is_write);
    if (!ret) {
        return NVME_SUCCESS;
    }

    return NVME_DNR;
}

static uint16_t nvme_dsm(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                         NvmeRequest *req)
{
    uint32_t dw10 = le32_to_cpu(cmd->cdw10);
    uint32_t dw11 = le32_to_cpu(cmd->cdw11);
    uint64_t prp1 = le64_to_cpu(cmd->dptr.prp1);
    uint64_t prp2 = le64_to_cpu(cmd->dptr.prp2);
    int i;

    if (dw11 & NVME_DSMGMT_AD) {
        uint16_t nr = (dw10 & 0xff) + 1;

        uint64_t slba;
        uint32_t nlb;
        NvmeDsmRange range[nr];

        if (dma_write_prp(n, (uint8_t *)range, sizeof(range), prp1, prp2)) {
            nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                                offsetof(NvmeCmd, dptr.prp1), 0, ns->id);
            return NVME_INVALID_FIELD | NVME_DNR;
        }

        req->status = NVME_SUCCESS;
        for (i = 0; i < nr; i++) {
            slba = le64_to_cpu(range[i].slba);
            nlb = le32_to_cpu(range[i].nlb);
            if (slba + nlb > le64_to_cpu(ns->id_ns.nsze)) {
                nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                                    offsetof(NvmeCmd, cdw10), slba + nlb, ns->id);
                return NVME_LBA_RANGE | NVME_DNR;
            }

            bitmap_clear(ns->util, slba, nlb);
        }
    }

    return NVME_SUCCESS;
}

static uint16_t nvme_compare(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                             NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint64_t prp1 = le64_to_cpu(rw->prp1);
    uint64_t prp2 = le64_to_cpu(rw->prp2);
    int i;

    uint64_t elba = slba + nlb;
    uint8_t lba_index = NVME_ID_NS_FLBAS_INDEX(ns->id_ns.flbas);
    uint8_t data_shift = ns->id_ns.lbaf[lba_index].lbads;
    uint64_t data_size = nlb << data_shift;
    uint64_t offset  = ns->start_block + (slba << data_shift);

    if ((slba + nlb) > le64_to_cpu(ns->id_ns.nsze)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), elba, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }
    if (n->id_ctrl.mdts && data_size > n->page_size * (1 << n->id_ctrl.mdts)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, nlb), nlb, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if (nvme_map_prp(&req->qsg, &req->iov, prp1, prp2, data_size, n)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, prp1), 0, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if (find_next_bit(ns->uncorrectable, elba, slba) < elba) {
        return NVME_UNRECOVERED_READ;
    }

    for (i = 0; i < req->qsg.nsg; i++) {
        uint32_t len = req->qsg.sg[i].len;
        uint8_t tmp[2][len];

        nvme_addr_read(n, req->qsg.sg[i].base, tmp[1], len);
        if (memcmp(tmp[0], tmp[1], len)) {
            qemu_sglist_destroy(&req->qsg);
            return NVME_CMP_FAILURE;
        }
        offset += len;
    }

    qemu_sglist_destroy(&req->qsg);

    return NVME_SUCCESS;
}

static uint16_t nvme_flush(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                           NvmeRequest *req)
{
    return NVME_SUCCESS;
}

static uint16_t nvme_write_zeros(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                                 NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;

    if ((slba + nlb) > ns->id_ns.nsze) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), slba + nlb, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }

    return NVME_SUCCESS;
}

static uint16_t nvme_write_uncor(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                                 NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;

    if ((slba + nlb) > ns->id_ns.nsze) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), slba + nlb, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }

    bitmap_set(ns->uncorrectable, slba, nlb);

    return NVME_SUCCESS;
}

static uint16_t nvme_io_cmd(FemuCtrl *n, NvmeCmd *cmd, NvmeRequest *req)
{
    NvmeNamespace *ns;
    uint32_t nsid = le32_to_cpu(cmd->nsid);

    if (nsid == 0 || nsid > n->num_namespaces) {
        femu_err("%s, NVME_INVALID_NSID %" PRIu32 "\n", __func__, nsid);
        return NVME_INVALID_NSID | NVME_DNR;
    }

    req->ns = ns = &n->namespaces[nsid - 1];

    switch (cmd->opcode) {
    case NVME_CMD_FLUSH:
        if (!n->id_ctrl.vwc || !n->features.volatile_wc) {
            return NVME_SUCCESS;
        }
        return nvme_flush(n, ns, cmd, req);
    case NVME_CMD_DSM:
        if (NVME_ONCS_DSM & n->oncs) {
            return nvme_dsm(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    case NVME_CMD_COMPARE:
        if (NVME_ONCS_COMPARE & n->oncs) {
            return nvme_compare(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    case NVME_CMD_WRITE_ZEROES:
        if (NVME_ONCS_WRITE_ZEROS & n->oncs) {
            return nvme_write_zeros(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    case NVME_CMD_WRITE_UNCOR:
        if (NVME_ONCS_WRITE_UNCORR & n->oncs) {
            return nvme_write_uncor(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
	case NVME_CMD_IO_MGMT_RECV: //update~
		return nvme_io_mgmt_recv(n, ns, cmd, req); // FIXME: No enable condition
	case NVME_CMD_IO_MGMT_SEND:					
		return nvme_io_mgmt_send(n, ns, cmd, req); // FIXME: No enable condition //~update
    default:
        if (n->ext_ops.io_cmd) {
            return n->ext_ops.io_cmd(n, ns, cmd, req);
        }

        femu_err("%s, NVME_INVALID_OPCODE\n", __func__);
        return NVME_INVALID_OPCODE | NVME_DNR;
    }
}

void nvme_post_cqes_io(void *opaque)
{
    NvmeCQueue *cq = opaque;
    NvmeRequest *req, *next;
    int64_t cur_time, ntt = 0;
    int processed = 0;

    QTAILQ_FOREACH_SAFE(req, &cq->req_list, entry, next) {
        if (nvme_cq_full(cq)) {
            break;
        }

        cur_time = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        if (cq->cqid != 0 && cur_time < req->expire_time) {
            ntt = req->expire_time;
            break;
        }

        nvme_post_cqe(cq, req);
        processed++;
    }

    if (ntt == 0) {
        ntt = qemu_clock_get_ns(QEMU_CLOCK_REALTIME) + CQ_POLLING_PERIOD_NS;
    }

    /* Only interrupt guest when we "do" complete some I/Os */
    if (processed > 0) {
        nvme_isr_notify_io(cq);
    }
}
