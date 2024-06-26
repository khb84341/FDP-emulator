#include "ftl.h"

//#define FEMU_DEBUG_FTL
//#define FDP_DEBUG

static void *ftl_thread(void *arg);

static inline bool should_gc(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines);
}

static inline bool should_gc_high(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines_high);
}

static inline bool should_fdp_gc(struct ssd *ssd, uint16_t rg) 	//update~
{
	return (ssd->rums[rg].free_ru_cnt <= ssd->sp.gc_thres_rus);
}																

static inline bool should_fdp_gc_high(struct ssd *ssd, uint16_t rg)
{
	return (ssd->rums[rg].free_ru_cnt <= ssd->sp.gc_thres_rus_high); 
}																//~update

static inline struct ppa get_maptbl_ent(struct ssd *ssd, uint64_t lpn)
{
    return ssd->maptbl[lpn];
}

static inline void set_maptbl_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    ftl_assert(lpn < ssd->sp.tt_pgs);
    ssd->maptbl[lpn] = *ppa;
}

static uint64_t ppa2pgidx(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t pgidx;

    pgidx = ppa->g.ch  * spp->pgs_per_ch  + \
            ppa->g.lun * spp->pgs_per_lun + \
            ppa->g.pl  * spp->pgs_per_pl  + \
            ppa->g.blk * spp->pgs_per_blk + \
            ppa->g.pg;

    ftl_assert(pgidx < spp->tt_pgs);

    return pgidx;
}

static inline uint64_t get_rmap_ent(struct ssd *ssd, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    return ssd->rmap[pgidx];
}

/* set rmap[page_no(ppa)] -> lpn */
static inline void set_rmap_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    ssd->rmap[pgidx] = lpn;
}

static inline int victim_line_cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
    return (next > curr);
}

static inline pqueue_pri_t victim_line_get_pri(void *a)
{
    return ((struct line *)a)->vpc;
}

static inline void victim_line_set_pri(void *a, pqueue_pri_t pri)
{
    ((struct line *)a)->vpc = pri;
}

static inline size_t victim_line_get_pos(void *a)
{
    return ((struct line *)a)->pos;
}

static inline void victim_line_set_pos(void *a, size_t pos)
{
    ((struct line *)a)->pos = pos;
}

static void ssd_init_lines(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *line;

    lm->tt_lines = spp->blks_per_pl;
    ftl_assert(lm->tt_lines == spp->tt_lines);
    lm->lines = g_malloc0(sizeof(struct line) * lm->tt_lines);

    QTAILQ_INIT(&lm->free_line_list);
    lm->victim_line_pq = pqueue_init(spp->tt_lines, victim_line_cmp_pri,
            victim_line_get_pri, victim_line_set_pri,
            victim_line_get_pos, victim_line_set_pos);
    QTAILQ_INIT(&lm->full_line_list);

    lm->free_line_cnt = 0;
    for (int i = 0; i < lm->tt_lines; i++) {
        line = &lm->lines[i];
        line->id = i;
        line->ipc = 0;
        line->vpc = 0;
        line->pos = 0;
        /* initialize all the lines as free lines */
        QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
        lm->free_line_cnt++;
    }

    ftl_assert(lm->free_line_cnt == lm->tt_lines);
    lm->victim_line_cnt = 0;
    lm->full_line_cnt = 0;
}

static inline int victim_ru_cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)//update~
{
    return (next > curr);
}

static inline pqueue_pri_t victim_ru_get_pri(void *a)		
{
    return ((struct ru *)a)->vpc;
}

static inline void victim_ru_set_pri(void *a, pqueue_pri_t pri)
{
    ((struct ru *)a)->vpc = pri;
}

static inline size_t victim_ru_get_pos(void *a)
{
    return ((struct ru *)a)->pos;
}

static inline void victim_ru_set_pos(void *a, size_t pos)
{
    ((struct ru *)a)->pos = pos;
}																		//~update 

static void ssd_init_fdp_ru_mgmts(struct ssd *ssd) 		//update~
{
	struct ssdparams *spp = &ssd->sp;
    struct fdp_ru_mgmt *rum = NULL;
    struct ru *ru = NULL;
	int nrg = spp->tt_luns / RG_DEGREE;
	int blkoff;
	
	ssd->rums = g_malloc(sizeof(struct fdp_ru_mgmt) * nrg);
	for (int i = 0; i < nrg; i++) {
		rum = &ssd->rums[i];
		rum->tt_rus = spp->blks_per_lun;
		assert(rum->tt_rus == spp->blks_per_lun);
		rum->rus = g_malloc0(sizeof(struct ru) * rum->tt_rus);

		QTAILQ_INIT(&rum->free_ru_list);
		rum->victim_ru_pq = pqueue_init(spp->tt_blks, victim_ru_cmp_pri,
            victim_ru_get_pri, victim_ru_set_pri,
            victim_ru_get_pos, victim_ru_set_pos);
		QTAILQ_INIT(&rum->full_ru_list);

		rum->free_ru_cnt = 0;

		for (int j = 0; j < rum->tt_rus; j++) {
			ru = &rum->rus[j];
			ru->id = j;
			ru->wp.ch = i * RG_DEGREE / spp->luns_per_ch;
			ru->wp.lun = i * RG_DEGREE % spp->luns_per_ch; 
			ru->wp.pl = 0;
			ru->wp.blk = j;
			ru->wp.pg = 0;
			ru->ipc = 0;
			ru->vpc = 0;
			ru->pos = 0; 
			ru->rut = RU_TYPE_NORMAL;
			/* ruh history for gc */	
			if (j < MAX_RUHS)
				ru->ruhid = j;

			blkoff = j % spp->blks_per_pl;
			//ru->blks = gmalloc0(sizeof(struct nand_block*) * RG_DEGREE); 
			for (int k = 0; k < RG_DEGREE; k++) {
				int cur_ch = ru->wp.ch + k / spp->luns_per_ch;
				int cur_lun = (ru->wp.lun + k) % spp->luns_per_ch;

				ru->blks[k] = &ssd->ch[cur_ch].lun[cur_lun].pl[0].blk[blkoff]; 
			} 

			/* initialize all the reclaim units as free reclaim units */
			QTAILQ_INSERT_TAIL(&rum->free_ru_list, ru, entry);
			rum->free_ru_cnt++;
		}

		ftl_assert(rum->free_ru_cnt == rum->tt_rus);
		rum->victim_ru_cnt = 0;
		rum->full_ru_cnt = 0; 
	} 
}														//~update

static int get_next_free_ruid(struct ssd *ssd, struct fdp_ru_mgmt *rum) //update~
{
#ifdef FDP_DEBUG
	printf("get_next_free_ruid() called -> ");
#endif
	struct ru *retru = NULL;

	retru = QTAILQ_FIRST(&rum->free_ru_list);
	if (!retru) {
		ftl_err("No free reclaim units left in [%s] !!!!\n", ssd->ssdname);
		return -1;
	}
#ifdef FDP_DEBUG 
	printf("new ru: %d\n", retru->id);
#endif
	QTAILQ_REMOVE(&rum->free_ru_list, retru, entry);
	rum->free_ru_cnt--;

	return retru->id; 
}																		//~update

static void ssd_init_fdp_ruhtbl(struct FemuCtrl *n, struct ssd *ssd) 	//update~
{
	NvmeEnduranceGroup *endgrp = &n->endgrps[0];
	struct ruh *ruh = NULL;
	struct fdp_ru_mgmt *rum = NULL;

	ssd->fdp_enabled = endgrp->fdp.enabled;
	ssd->ruhtbl = g_malloc0(sizeof(struct ruh) * endgrp->fdp.nruh); 
	for (int i = 0; i < endgrp->fdp.nruh; i++) {
		ruh = &ssd->ruhtbl[i];
		ruh->ruht = endgrp->fdp.ruhs[i].ruht;
		ruh->cur_ruids = g_malloc0(sizeof(int) * endgrp->fdp.nrg);
		ruh->pi_gc_ruids = g_malloc0(sizeof(int) * endgrp->fdp.nrg);
		for (int j = 0; j < endgrp->fdp.nrg; j++)  {
			rum = &ssd->rums[j];
			ruh->cur_ruids[j] = get_next_free_ruid(ssd, rum);
		} 
	} 

	/* reserve one ru for initially isolated gc */
	for (int i = 0; i < endgrp->fdp.nrg; i++) {
		rum = &ssd->rums[i];
		rum->ii_gc_ruid = get_next_free_ruid(ssd, rum);
		rum->rus[rum->ii_gc_ruid].rut = RU_TYPE_II_GC;
	}

	/* reserve rus for persistently isolated gc */
	/* FIXME: if you want to set the ruhs as PI, activate the code below */
	
	/*
	for (int i = 0; i < endgrp->fdp.nrg; i++) {
		int pi_gc_ruid;
		rum = &ssd->rums[i]; 
		
		for (int j = 0; j < MAX_RUHS; j++) {
			pi_gc_ruid = get_next_free_ruid(ssd, rum);
			ssd->ruhtbl[j].pi_gc_ruids[i] = pi_gc_ruid;
			rum->rus[pi_gc_ruid].rut = RU_TYPE_PI_GC;
		} 
	} */
}																		//~update

static void ssd_init_write_pointer(struct ssd *ssd)
{
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);
    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;

    /* wpp->curline is always our next-to-write super-block */
    wpp->curline = curline;
    wpp->ch = 0;
    wpp->lun = 0;
    wpp->pg = 0;
    wpp->blk = 0;
    wpp->pl = 0;
}

static inline void check_addr(int a, int max)
{
    ftl_assert(a >= 0 && a < max);
}

static struct line *get_next_free_line(struct ssd *ssd)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);
    if (!curline) {
        ftl_err("No free lines left in [%s] !!!!\n", ssd->ssdname);
        return NULL;
    }

    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;
    return curline;
}

static void ssd_advance_write_pointer(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;

    check_addr(wpp->ch, spp->nchs);
    wpp->ch++;
    if (wpp->ch == spp->nchs) {
        wpp->ch = 0;
        check_addr(wpp->lun, spp->luns_per_ch);
        wpp->lun++;
        /* in this case, we should go to next lun */
        if (wpp->lun == spp->luns_per_ch) {
            wpp->lun = 0;
            /* go to next page in the block */
            check_addr(wpp->pg, spp->pgs_per_blk);
            wpp->pg++;
            if (wpp->pg == spp->pgs_per_blk) {
                wpp->pg = 0;
                /* move current line to {victim,full} line list */
                if (wpp->curline->vpc == spp->pgs_per_line) {
                    /* all pgs are still valid, move to full line list */
                    ftl_assert(wpp->curline->ipc == 0);
                    QTAILQ_INSERT_TAIL(&lm->full_line_list, wpp->curline, entry);
                    lm->full_line_cnt++;
                } else {
                    ftl_assert(wpp->curline->vpc >= 0 && wpp->curline->vpc < spp->pgs_per_line);
                    /* there must be some invalid pages in this line */
                    ftl_assert(wpp->curline->ipc > 0);
                    pqueue_insert(lm->victim_line_pq, wpp->curline);
                    lm->victim_line_cnt++;
                }
                /* current line is used up, pick another empty line */
                check_addr(wpp->blk, spp->blks_per_pl);
                wpp->curline = NULL;
                wpp->curline = get_next_free_line(ssd);
                if (!wpp->curline) {
                    /* TODO */
                    abort();
                }
                wpp->blk = wpp->curline->id;
                check_addr(wpp->blk, spp->blks_per_pl);
                /* make sure we are starting from page 0 in the super block */
                ftl_assert(wpp->pg == 0);
                ftl_assert(wpp->lun == 0);
                ftl_assert(wpp->ch == 0);
                /* TODO: assume # of pl_per_lun is 1, fix later */
                ftl_assert(wpp->pl == 0);
            }
        }
    }
}

// #define SMALL_RG
static void ssd_advance_fdp_write_pointer(struct ssd *ssd, uint16_t rgid, int lpn, uint16_t ruhid, bool for_gc) //update~
{
	struct ssdparams *spp = &ssd->sp;
	struct fdp_ru_mgmt *rum = &ssd->rums[rgid];
	struct ruh *ruh = &ssd->ruhtbl[ruhid];
	int max_ch = (rgid + 1) * (RG_DEGREE / spp->luns_per_ch);
	int ruid;
	struct ru *ru = NULL;
	if (for_gc) {
		if (ruh->ruht == NVME_RUHT_INITIALLY_ISOLATED) {
			ruid = rum->ii_gc_ruid; 
		}
		else if (ruh->ruht == NVME_RUHT_PERSISTENTLY_ISOLATED) {
			ruid = ruh->pi_gc_ruids[rgid];
		}
		else if (ruh->ruht == NVME_RUHT_PERSISTENTLY_ISOLATED_NO_SEPARATED_BLK)
			ruid = ruh->cur_ruids[rgid];
		else {
			if (ssd->gc_cnt[lpn] == 0)
				ruid = ruh->pi_gc_ruids[rgid];
			else 
				ruid = rum->ii_gc_ruid;
		}
	}
	else
		ruid = ruh->cur_ruids[rgid];

	ru = &rum->rus[ruid]; 

	/* Case that an RG has more than two channels -> same with the origin */
#ifdef SMALL_RG
	if (RG_DEGREE > spp->luns_per_ch)
	{ 
#endif
		check_addr(ru->wp.ch, max_ch);
		ru->wp.ch++;
		if (ru->wp.ch == max_ch) {
			ru->wp.ch = rgid * (RG_DEGREE / spp->luns_per_ch);
			check_addr(ru->wp.lun, spp->luns_per_ch);
			ru->wp.lun++;
			/* in this case, we should go to next lun */
			if (ru->wp.lun == spp->luns_per_ch) {
				ru->wp.lun = 0;
				/* go to next page in the block */
				check_addr(ru->wp.pg, spp->pgs_per_blk);
				ru->wp.pg++;
				if (ru->wp.pg == spp->pgs_per_blk) {
					ru->wp.pg = 0;
					/* move current ru to {victim,full} ru list */
					if (ru->vpc == spp->pgs_per_ru) {
						/* all pgs are still valid, move to full ru list */
						ftl_assert(ru->ipc == 0);
						QTAILQ_INSERT_TAIL(&rum->full_ru_list, ru, entry);
						rum->full_ru_cnt++;
					} else {
						ftl_assert(ru->vpc >= 0 && ru->vpc < spp->pgs_per_ru);
						ftl_assert(ru->ipc > 0);
						pqueue_insert(rum->victim_ru_pq, ru);
						rum->victim_ru_cnt++;
					}
					/* current ru is used up, pick another empty ru */
					check_addr(ru->wp.blk, spp->blks_per_pl); 
					/* ruh history for gc later */
					ru->ruhid = ruhid; 
					if (ru->rut == RU_TYPE_II_GC) {
						rum->ii_gc_ruid = get_next_free_ruid(ssd, rum);
						rum->rus[rum->ii_gc_ruid].rut = RU_TYPE_II_GC;
					}
					else if (ru->rut == RU_TYPE_PI_GC) {
						ruh->pi_gc_ruids[rgid] = get_next_free_ruid(ssd, rum);
						rum->rus[ruh->pi_gc_ruids[rgid]].rut = RU_TYPE_PI_GC;
					}
					else {
						/* update ruhtbl */
						ruh->cur_ruids[rgid] = get_next_free_ruid(ssd, rum);
						rum->rus[ruh->cur_ruids[rgid]].rut = RU_TYPE_NORMAL;
					} 
					check_addr(ru->wp.blk, spp->blks_per_pl);
					/* make sure we are starting from page 0 in the ru */
					ftl_assert(ru->wp.pg == 0);
					ftl_assert(ru->wp.lun == 0);
					ftl_assert(ru->wp.ch == rgid * (RG_DEGREE / spp->luns_per_ch));
					/* TODO: assume # of pl_per_lun is 1, fix later */
					ftl_assert(ru->wp.pl == 0);
				}
			}
		} 
#ifdef SMALL_RG
	} 
#endif

#ifdef SMALL_RG
	/* Case that an RG is included in one channel */
	else 
	{
		check_addr(ru->wp.lun, spp->luns_per_ch); 
		ru->wp.lun++;
		/* move to next page */
		if (ru->wp.lun % RG_DEGREE == 0)
		{
			ru->wp.lun = 0;
			check_addr(ru->wp.pg, spp->pgs_per_ch);
			ru->wp.pg++;
			if (ru->wp.pg == spp->pgs_per_blk)
			{
				ru->wp.pg = 0;
				/* move current ru to {victim,full} ru list */
				if (ru->vpc == spp->pgs_per_blk * RG_DEGREE)
				{
					/* all pgs are still valid, move to full ru list */
					ftl_assert(ru->ipc == 0);
					QTAILQ_INSERT_TAIL(&rum->full_ru_list, ru, entry);
					rum->full_ru_cnt++;
				}
				else
				{
					/* there must be some invalid pages in this ru */
					ftl_assert(ru->vpn >= 0 && ru->vpc < RG_DEGREE * spp->pgs_per_blk);
					ftl_assert(ru->ipc > 0);
					pqueue_insert(rum->victim_ru_pq, ru);
					rum->victim_ru_cnt++;
				}
				/* current ru is used up, pick another empty ru */ 
				ru->ruhid = ruhid;
				if (ru->rut == RU_TYPE_II_GC) {
					rum->ii_gc_ruid = get_next_free_ruid(ssd, rum);
					rum->rus[rum->ii_gc_ruid].rut = RU_TYPE_II_GC;
				}
				else if (ru->rut == RU_TYPE_PI_GC) {
					ruh->pi_gc_ruids[rgid] = get_next_free_ruid(ssd, rum);
					rum->rus[ruh->pi_gc_ruids[rgid]].rut = RU_TYPE_PI_GC;
				}
				else {
					/* update ruhtbl */
					ruh->cur_ruids[rgid] = get_next_free_ruid(ssd, rum);
					rum->rus[ruh->cur_ruids[rgid]].rut = RU_TYPE_NORMAL;
				} 
				check_addr(ru->wp.blk, spp->blks_per_pl);
				/* make sure we are starting from page 0 in the ru */
				ftl_assert(ru->wp.pg == 0);
				ftl_assert(ru->wp.lun == 0);
				/* TODO: assume # of pl_per_lun is 1, fix later */
				ftl_assert(ru->wp.pl == 0);
			}
		}
	}
#endif
}																					

static struct ppa fdp_get_new_page(struct ssd *ssd, uint16_t rgid, 
		int lpn, uint16_t ruhid, bool for_gc) //update
{
	struct fdp_ru_mgmt *rum = &ssd->rums[rgid];
	struct ruh* ruh = &ssd->ruhtbl[ruhid];
	int ruid;
	struct ru *ru;
    struct ppa ppa; 
	if (for_gc) {
		if (ruh->ruht == NVME_RUHT_INITIALLY_ISOLATED)
			ruid = rum->ii_gc_ruid;
		else if (ruh->ruht == NVME_RUHT_PERSISTENTLY_ISOLATED)
			ruid = ruh->pi_gc_ruids[rgid];
		else if (ruh->ruht == NVME_RUHT_PERSISTENTLY_ISOLATED_NO_SEPARATED_BLK) {
			ruid = ruh->cur_ruids[rgid];
		}
		else { 
			if (ssd->gc_cnt[lpn] == 0)
				ruid = ruh->pi_gc_ruids[rgid];
			else 
				ruid = rum->ii_gc_ruid;
		}
	}
	else // normal
		ruid = ruh->cur_ruids[rgid];

	ppa.ppa = 0;

	rum = &ssd->rums[rgid];
	ru = &rum->rus[ruid];
	ppa.g.ch = ru->wp.ch;
	ppa.g.lun = ru->wp.lun;
	ppa.g.pg = ru->wp.pg;
	ppa.g.blk = ru->wp.blk;
	ppa.g.pl = ru->wp.pl; 

    ftl_assert(ppa.g.pl == 0);

    return ppa;
}																		//~update

static struct ppa get_new_page(struct ssd *ssd)
{
    struct write_pointer *wpp;
    struct ppa ppa;

    ppa.ppa = 0;

	wpp = &ssd->wp; 
	ppa.g.ch = wpp->ch;
	ppa.g.lun = wpp->lun;
	ppa.g.pg = wpp->pg;
	ppa.g.blk = wpp->blk;
	ppa.g.pl = wpp->pl;

    ftl_assert(ppa.g.pl == 0);

    return ppa;
}

static void check_params(struct ssdparams *spp)
{
    /*
     * we are using a general write pointer increment method now, no need to
     * force luns_per_ch and nchs to be power of 2
     */

    //ftl_assert(is_power_of_2(spp->luns_per_ch));
    //ftl_assert(is_power_of_2(spp->nchs));
}

static void ssd_init_params(struct ssdparams *spp, FemuCtrl *n)
{
    spp->secsz = n->bb_params.secsz; // 512
    spp->secs_per_pg = n->bb_params.secs_per_pg; // 8
    spp->pgs_per_blk = n->bb_params.pgs_per_blk; //256
    spp->blks_per_pl = n->bb_params.blks_per_pl; /* 256 16GB */
    spp->pls_per_lun = n->bb_params.pls_per_lun; // 1
    spp->luns_per_ch = n->bb_params.luns_per_ch; // 8
    spp->nchs = n->bb_params.nchs; // 8

    spp->pg_rd_lat = n->bb_params.pg_rd_lat;
    spp->pg_wr_lat = n->bb_params.pg_wr_lat;
    spp->blk_er_lat = n->bb_params.blk_er_lat;
    spp->ch_xfer_lat = n->bb_params.ch_xfer_lat;

    /* calculated values */
    spp->secs_per_blk = spp->secs_per_pg * spp->pgs_per_blk;
    spp->secs_per_pl = spp->secs_per_blk * spp->blks_per_pl;
    spp->secs_per_lun = spp->secs_per_pl * spp->pls_per_lun;
    spp->secs_per_ch = spp->secs_per_lun * spp->luns_per_ch;
    spp->tt_secs = spp->secs_per_ch * spp->nchs;

    spp->pgs_per_pl = spp->pgs_per_blk * spp->blks_per_pl;
    spp->pgs_per_lun = spp->pgs_per_pl * spp->pls_per_lun;
    spp->pgs_per_ch = spp->pgs_per_lun * spp->luns_per_ch;
    spp->tt_pgs = spp->pgs_per_ch * spp->nchs;
#ifdef DEVICE_UTIL_DEBUG
	spp->tt_valid_pgs = 0; //update
#endif

    spp->blks_per_lun = spp->blks_per_pl * spp->pls_per_lun;
    spp->blks_per_ch = spp->blks_per_lun * spp->luns_per_ch;
    spp->tt_blks = spp->blks_per_ch * spp->nchs;

    spp->pls_per_ch =  spp->pls_per_lun * spp->luns_per_ch;
    spp->tt_pls = spp->pls_per_ch * spp->nchs;

    spp->tt_luns = spp->luns_per_ch * spp->nchs;

    /* line is special, put it at the end */
    spp->blks_per_line = spp->tt_luns; /* TODO: to fix under multiplanes */
    spp->pgs_per_line = spp->blks_per_line * spp->pgs_per_blk;
    spp->secs_per_line = spp->pgs_per_line * spp->secs_per_pg;
    spp->tt_lines = spp->blks_per_lun; /* TODO: to fix under multiplanes */

    spp->blks_per_ru = RG_DEGREE; 									//update~
    spp->pgs_per_ru = spp->blks_per_ru * spp->pgs_per_blk;
    spp->secs_per_ru = spp->pgs_per_ru * spp->secs_per_pg;
	spp->chs_per_ru = RG_DEGREE / spp->luns_per_ch;
	spp->luns_per_ru = spp->blks_per_ru;
    spp->tt_rus = spp->blks_per_lun;								//~update
	
    spp->gc_thres_pcent = n->bb_params.gc_thres_pcent/100.0;
    spp->gc_thres_lines = (int)((1 - spp->gc_thres_pcent) * spp->tt_lines);
    spp->gc_thres_rus = (int)((1 - spp->gc_thres_pcent) * spp->tt_rus);
    spp->gc_thres_pcent_high = n->bb_params.gc_thres_pcent_high/100.0; // update
    spp->gc_thres_lines_high = (int)((1 - spp->gc_thres_pcent_high) * spp->tt_lines);
    spp->gc_thres_rus_high = (int)((1 - spp->gc_thres_pcent_high) * spp->tt_rus); // update
    spp->enable_gc_delay = true; 

    check_params(spp);
}

static void ssd_init_nand_page(struct nand_page *pg, struct ssdparams *spp)
{
    pg->nsecs = spp->secs_per_pg;
    pg->sec = g_malloc0(sizeof(nand_sec_status_t) * pg->nsecs);
    for (int i = 0; i < pg->nsecs; i++) {
        pg->sec[i] = SEC_FREE;
    }
    pg->status = PG_FREE;
}

static void ssd_init_nand_blk(struct nand_block *blk, struct ssdparams *spp)
{
    blk->npgs = spp->pgs_per_blk;
    blk->pg = g_malloc0(sizeof(struct nand_page) * blk->npgs);
    for (int i = 0; i < blk->npgs; i++) {
        ssd_init_nand_page(&blk->pg[i], spp);
    }
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt = 0;
    blk->wp = 0;
}

static void ssd_init_nand_plane(struct nand_plane *pl, struct ssdparams *spp)
{
    pl->nblks = spp->blks_per_pl;
    pl->blk = g_malloc0(sizeof(struct nand_block) * pl->nblks);
    for (int i = 0; i < pl->nblks; i++) {
        ssd_init_nand_blk(&pl->blk[i], spp);
    }
}

static void ssd_init_nand_lun(struct nand_lun *lun, struct ssdparams *spp)
{
    lun->npls = spp->pls_per_lun;
    lun->pl = g_malloc0(sizeof(struct nand_plane) * lun->npls);
    for (int i = 0; i < lun->npls; i++) {
        ssd_init_nand_plane(&lun->pl[i], spp);
    }
    lun->next_lun_avail_time = 0;
    lun->busy = false;
}

static void ssd_init_ch(struct ssd_channel *ch, struct ssdparams *spp)
{
    ch->nluns = spp->luns_per_ch;
    ch->lun = g_malloc0(sizeof(struct nand_lun) * ch->nluns);
    for (int i = 0; i < ch->nluns; i++) {
        ssd_init_nand_lun(&ch->lun[i], spp);
    }
    ch->next_ch_avail_time = 0;
    ch->busy = 0;
}

static void ssd_init_maptbl(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->maptbl = g_malloc0(sizeof(struct ppa) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->maptbl[i].ppa = UNMAPPED_PPA;
    }
}

static void ssd_init_rmap(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->rmap = g_malloc0(sizeof(uint64_t) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->rmap[i] = INVALID_LPN;
    }
}

void ssd_init(FemuCtrl *n)
{
    struct ssd *ssd = n->ssd;
    struct ssdparams *spp = &ssd->sp;

    ftl_assert(ssd);
#ifdef UPDATE_FREQ
	/*for (int i = 0; i < 4; i++)
		for (int j = 0; j < 24; j++)
			ssd->ten[i][j] = 0;*/
	memset(ssd->ten, 0x00, NR_TENANTS * sizeof(struct tenant));
#endif

    ssd_init_params(spp, n);

    /* initialize ssd internal layout architecture */
    ssd->ch = g_malloc0(sizeof(struct ssd_channel) * spp->nchs);
    for (int i = 0; i < spp->nchs; i++) {
        ssd_init_ch(&ssd->ch[i], spp);
    } 

    /* initialize maptbl */
    ssd_init_maptbl(ssd);

    /* initialize rmap */
    ssd_init_rmap(ssd);

#ifdef WAF_TEST
	n->host_writes = 0;
	n->gc_writes = 0;
#endif

	ssd->gc_cnt = g_malloc0(sizeof(int) * spp->tt_pgs);  // update

    /* initialize all the lines */
    ssd_init_lines(ssd);

    ssd_init_write_pointer(ssd);

	ssd_init_fdp_ru_mgmts(ssd); 						//update~

	ssd_init_fdp_ruhtbl(n, ssd);							//~update
    /* initialize write pointer, this is how we allocate new pages for writes */

    qemu_thread_create(&ssd->ftl_thread, "FEMU-FTL-Thread", ftl_thread, n,
                       QEMU_THREAD_JOINABLE);
}

static inline bool valid_ppa(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    int ch = ppa->g.ch;
    int lun = ppa->g.lun;
    int pl = ppa->g.pl;
    int blk = ppa->g.blk;
    int pg = ppa->g.pg;
    int sec = ppa->g.sec;

    if (ch >= 0 && ch < spp->nchs && lun >= 0 && lun < spp->luns_per_ch && pl >=
        0 && pl < spp->pls_per_lun && blk >= 0 && blk < spp->blks_per_pl && pg
        >= 0 && pg < spp->pgs_per_blk && sec >= 0 && sec < spp->secs_per_pg)
        return true;

    return false;
}

static inline bool valid_lpn(struct ssd *ssd, uint64_t lpn)
{
    return (lpn < ssd->sp.tt_pgs);
}

static inline bool mapped_ppa(struct ppa *ppa)
{
    return !(ppa->ppa == UNMAPPED_PPA);
}

static inline struct ssd_channel *get_ch(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->ch[ppa->g.ch]);
}

static inline struct nand_lun *get_lun(struct ssd *ssd, struct ppa *ppa)
{
    struct ssd_channel *ch = get_ch(ssd, ppa);
    return &(ch->lun[ppa->g.lun]);
}

static inline struct nand_plane *get_pl(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_lun *lun = get_lun(ssd, ppa);
    return &(lun->pl[ppa->g.pl]);
}

static inline struct nand_block *get_blk(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_plane *pl = get_pl(ssd, ppa);
    return &(pl->blk[ppa->g.blk]);
}

static inline struct line *get_line(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->lm.lines[ppa->g.blk]);
}

static inline struct ru *get_ru(struct ssd *ssd, struct ppa *ppa)
{
	struct ssdparams *spp = &ssd->sp;
	//uint16_t rgid = ppa->g.ch * (spp->luns_per_ch / RG_DEGREE) + (ppa->g.lun / RG_DEGREE);
	uint16_t rgid = (ppa->g.ch * spp->luns_per_ch + ppa->g.lun) / RG_DEGREE;
	struct fdp_ru_mgmt *rum = &ssd->rums[rgid];
	uint16_t ruid = ppa->g.blk;

    return &rum->rus[ruid];
}

static inline struct nand_page *get_pg(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = get_blk(ssd, ppa);
    return &(blk->pg[ppa->g.pg]);
}

static uint64_t ssd_advance_status(struct ssd *ssd, struct ppa *ppa, struct
        nand_cmd *ncmd)
{
    int c = ncmd->cmd;
    uint64_t cmd_stime = (ncmd->stime == 0) ? \
        qemu_clock_get_ns(QEMU_CLOCK_REALTIME) : ncmd->stime;
    uint64_t nand_stime;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lun = get_lun(ssd, ppa);
    uint64_t lat = 0;

    switch (c) {
    case NAND_READ:
        /* read: perform NAND cmd first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;
        lat = lun->next_lun_avail_time - cmd_stime;
#if 0
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;

        /* read: then data transfer through channel */
        chnl_stime = (ch->next_ch_avail_time < lun->next_lun_avail_time) ? \
            lun->next_lun_avail_time : ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        lat = ch->next_ch_avail_time - cmd_stime;
#endif
        break;

    case NAND_WRITE:
        /* write: transfer data through channel first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        if (ncmd->type == USER_IO) {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        } else {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        }
        lat = lun->next_lun_avail_time - cmd_stime;

#if 0
        chnl_stime = (ch->next_ch_avail_time < cmd_stime) ? cmd_stime : \
                     ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        /* write: then do NAND program */
        nand_stime = (lun->next_lun_avail_time < ch->next_ch_avail_time) ? \
            ch->next_ch_avail_time : lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
#endif
        break;

    case NAND_ERASE:
        /* erase: only need to advance NAND status */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->blk_er_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
        break;

    default:
        ftl_err("Unsupported NAND command: 0x%x\n", c);
    }

    return lat;
}

/* update SSD status about one page from PG_VALID -> PG_VALID */
static void mark_page_invalid(struct ssd *ssd, struct ppa *ppa, uint16_t rgid)
{
    struct line_mgmt *lm = &ssd->lm;
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    bool was_full_line = false;
    bool was_full_ru = false;	//update
	struct fdp_ru_mgmt *rum = &ssd->rums[rgid];//update
    struct line *line;
    struct ru *ru; //update

    /* update corresponding page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_VALID);
    pg->status = PG_INVALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    ftl_assert(blk->ipc >= 0 && blk->ipc < spp->pgs_per_blk);
    blk->ipc++;
    ftl_assert(blk->vpc > 0 && blk->vpc <= spp->pgs_per_blk);
    blk->vpc--;

	if (ssd->fdp_enabled)	//update
	{ 
		/* update corresponding ru status */
		ru = get_ru(ssd, ppa);
		ftl_assert(ru->ipc >= 0 && ru->ipc < spp->pgs_per_ru);
		if (ru->vpc == spp->pgs_per_ru) {
			ftl_assert(ru->ipc == 0);
			was_full_ru = true;
		}
		ru->ipc++;
		ftl_assert(ru->vpc > 0 && ru->vpc <= spp->pgs_per_ru);
		/* Adjust the position of the victime ru in the pq under over-writes */
		if (ru->pos) {
			/* Note that ru->vpc will be updated by this call */
			pqueue_change_priority(rum->victim_ru_pq, ru->vpc - 1, ru);
		} else {
			ru->vpc--;
		}

		if (was_full_ru) {
			/* move ru: "full" -> "victim" */
			QTAILQ_REMOVE(&rum->full_ru_list, ru, entry);
			rum->full_ru_cnt--;
			pqueue_insert(rum->victim_ru_pq, ru);
			rum->victim_ru_cnt++;
		}
	}
	else
	{
		/* update corresponding line status */
		line = get_line(ssd, ppa);
		ftl_assert(line->ipc >= 0 && line->ipc < spp->pgs_per_line);
		if (line->vpc == spp->pgs_per_line) {
			ftl_assert(line->ipc == 0);
			was_full_line = true;
		}
		line->ipc++;
		ftl_assert(line->vpc > 0 && line->vpc <= spp->pgs_per_line);
		/* Adjust the position of the victime line in the pq under over-writes */
		if (line->pos) {
			/* Note that line->vpc will be updated by this call */
			pqueue_change_priority(lm->victim_line_pq, line->vpc - 1, line);
		} else {
			line->vpc--;
		}

		if (was_full_line) {
			/* move line: "full" -> "victim" */
			QTAILQ_REMOVE(&lm->full_line_list, line, entry);
			lm->full_line_cnt--;
			pqueue_insert(lm->victim_line_pq, line);
			lm->victim_line_cnt++;
		} 
	}
}

static void mark_page_valid(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    struct line *line;
    struct ru *ru;

    /* update page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_FREE);
    pg->status = PG_VALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    ftl_assert(blk->vpc >= 0 && blk->vpc < ssd->sp.pgs_per_blk);
    blk->vpc++;

    /* update corresponding ru status */
	if (ssd->fdp_enabled) {
		ru = get_ru(ssd, ppa);
		ftl_assert(ru->vpc >= 0 && ru->vpc < ssd->sp.pgs_per_ru);
		ru->vpc++; 
	}
    /* update corresponding line status */
	else {
		line = get_line(ssd, ppa);
		ftl_assert(line->vpc >= 0 && line->vpc < ssd->sp.pgs_per_line);
		line->vpc++;
	}
} 

static void mark_block_free(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = get_blk(ssd, ppa);
    struct nand_page *pg = NULL;

    for (int i = 0; i < spp->pgs_per_blk; i++) {
        /* reset page status */
        pg = &blk->pg[i];
        ftl_assert(pg->nsecs == spp->secs_per_pg);
        pg->status = PG_FREE;
    }

    /* reset block status */
    ftl_assert(blk->npgs == spp->pgs_per_blk);
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt++;
}

static void gc_read_page(struct ssd *ssd, struct ppa *ppa)
{
    /* advance ssd status, we don't care about how long it takes */
    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcr;
        gcr.type = GC_IO;
        gcr.cmd = NAND_READ;
        gcr.stime = 0;
        ssd_advance_status(ssd, ppa, &gcr);
    }
}

/* move valid page data (already in DRAM) from victim line to a new page */
static uint64_t fdp_gc_write_page(struct ssd *ssd, struct ppa *old_ppa, uint16_t rgid, uint16_t ruhid)
{
    struct ppa new_ppa;
    struct nand_lun *new_lun;
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);

    ftl_assert(valid_lpn(ssd, lpn));

	new_ppa = fdp_get_new_page(ssd, rgid, lpn, ruhid, true);
    /* update maptbl */

#ifdef FDP_DEBUG
	printf("gc -> rgid: %d ruhid: %d new ch: %d new lun: %d new blk: %d new pg: %d\n", 
			rgid, ruhid, new_ppa.g.ch, new_ppa.g.lun, new_ppa.g.blk, new_ppa.g.pg);
#endif

    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

	mark_page_valid(ssd, &new_ppa);

    /* need to advance the write pointer here */
	ssd_advance_fdp_write_pointer(ssd, rgid, lpn, ruhid, true);

	ssd->gc_cnt[lpn]++;
	
    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcw;
        gcw.type = GC_IO;
        gcw.cmd = NAND_WRITE;
        gcw.stime = 0;
        ssd_advance_status(ssd, &new_ppa, &gcw);
    }

    /* advance per-ch gc_endtime as well */
#if 0
    new_ch = get_ch(ssd, &new_ppa);
    new_ch->gc_endtime = new_ch->next_ch_avail_time;
#endif

    new_lun = get_lun(ssd, &new_ppa);
    new_lun->gc_endtime = new_lun->next_lun_avail_time;

    return 0;
}

/* move valid page data (already in DRAM) from victim line to a new page */
static uint64_t gc_write_page(struct ssd *ssd, struct ppa *old_ppa)
{
    struct ppa new_ppa;
    struct nand_lun *new_lun;
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);

    ftl_assert(valid_lpn(ssd, lpn));

	new_ppa = get_new_page(ssd);
    /* update maptbl */

#ifdef FEMU_DEBUG_FTL
	printf("ch: %d, lun: %d, blk: %d, pg: %d\n", 
			new_ppa.g.ch, new_ppa.g.lun, new_ppa.g.blk, new_ppa.g.pg);
#endif

    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

	mark_page_valid(ssd, &new_ppa);

    /* need to advance the write pointer here */
	ssd_advance_write_pointer(ssd);

    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcw;
        gcw.type = GC_IO;
        gcw.cmd = NAND_WRITE;
        gcw.stime = 0;
        ssd_advance_status(ssd, &new_ppa, &gcw);
    }

    /* advance per-ch gc_endtime as well */
#if 0
    new_ch = get_ch(ssd, &new_ppa);
    new_ch->gc_endtime = new_ch->next_ch_avail_time;
#endif

    new_lun = get_lun(ssd, &new_ppa);
    new_lun->gc_endtime = new_lun->next_lun_avail_time;

    return 0;
}

static struct line *select_victim_line(struct ssd *ssd, bool force)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *victim_line = NULL;

    victim_line = pqueue_peek(lm->victim_line_pq);
    if (!victim_line) {
        return NULL;
    }

    if (!force && victim_line->ipc < ssd->sp.pgs_per_line / 8) {
        return NULL;
    }

    pqueue_pop(lm->victim_line_pq);
    victim_line->pos = 0;
    lm->victim_line_cnt--;

    /* victim_line is a danggling node now */
    return victim_line;
}

static struct ru *select_victim_ru(struct ssd *ssd, bool force, int rgid)
{
    struct fdp_ru_mgmt *rum = &ssd->rums[rgid];
    struct ru *victim_ru = NULL;

    victim_ru = pqueue_peek(rum->victim_ru_pq);
    if (!victim_ru) {
        return NULL;
    } 

    if (!force && victim_ru->ipc < ssd->sp.pgs_per_ru / 8) {
        return NULL;
    }

    pqueue_pop(rum->victim_ru_pq);
    victim_ru->pos = 0;
    rum->victim_ru_cnt--;

    /* victim_ru is a danggling node now */
    return victim_ru;
}

/* here ppa identifies the block we want to clean */
static int fdp_clean_one_block(struct ssd *ssd, struct ppa *ppa, uint16_t rgid, uint16_t ruhid)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_page *pg_iter = NULL;
    int cnt = 0;

    for (int pg = 0; pg < spp->pgs_per_blk; pg++) {
        ppa->g.pg = pg;
#ifdef FDP_DEBUG
	printf("old_ch: %d old_lun: %d old_pl: %d old_blk: %d old_pg: %d\n",
			ppa->g.ch, ppa->g.lun, ppa->g.pl, ppa->g.blk, ppa->g.pg);
#endif
        pg_iter = get_pg(ssd, ppa);
        /* there shouldn't be any free page in victim blocks */
        ftl_assert(pg_iter->status != PG_FREE);
        if (pg_iter->status == PG_VALID) {
            gc_read_page(ssd, ppa);
            /* delay the maptbl update until "write" happens */
            fdp_gc_write_page(ssd, ppa, rgid, ruhid);
            cnt++;
        }
    }

    ftl_assert(get_blk(ssd, ppa)->vpc == cnt);
	return cnt;
}

/* here ppa identifies the block we want to clean */
static void clean_one_block(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_page *pg_iter = NULL;
    int cnt = 0;

    for (int pg = 0; pg < spp->pgs_per_blk; pg++) {
        ppa->g.pg = pg;
        pg_iter = get_pg(ssd, ppa);
        /* there shouldn't be any free page in victim blocks */
        ftl_assert(pg_iter->status != PG_FREE);
        if (pg_iter->status == PG_VALID) {
            gc_read_page(ssd, ppa);
            /* delay the maptbl update until "write" happens */
            gc_write_page(ssd, ppa);
            cnt++;
        }
    }

    ftl_assert(get_blk(ssd, ppa)->vpc == cnt);
}

static void mark_line_free(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *line = get_line(ssd, ppa);
    line->ipc = 0;
    line->vpc = 0;
    /* move this line to free line list */
    QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
    lm->free_line_cnt++;
}

static int do_gc(struct ssd *ssd, bool force)
{
    struct line *victim_line = NULL;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lunp;
    struct ppa ppa;
    int ch, lun;

    victim_line = select_victim_line(ssd, force);
    if (!victim_line) {
        return -1;
    }

    ppa.g.blk = victim_line->id;
    ftl_debug("GC-ing line:%d,ipc=%d,victim=%d,full=%d,free=%d\n", ppa.g.blk,
              victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
              ssd->lm.free_line_cnt);

    /* copy back valid data */
    for (ch = 0; ch < spp->nchs; ch++) {
        for (lun = 0; lun < spp->luns_per_ch; lun++) {
            ppa.g.ch = ch;
            ppa.g.lun = lun;
            ppa.g.pl = 0;
            lunp = get_lun(ssd, &ppa);
            clean_one_block(ssd, &ppa);
            mark_block_free(ssd, &ppa);

            if (spp->enable_gc_delay) {
                struct nand_cmd gce;
                gce.type = GC_IO;
                gce.cmd = NAND_ERASE;
                gce.stime = 0;
                ssd_advance_status(ssd, &ppa, &gce);
            }

            lunp->gc_endtime = lunp->next_lun_avail_time;
        }
    }

    /* update line status */
    mark_line_free(ssd, &ppa);

    return 0;
}

bool gc = 0; 
static int do_fdp_gc(struct ssd *ssd, uint16_t rgid, bool force, NvmeRequest *req)
{
	if (!gc)
		printf("do_fdp_gc() called\n");
	gc = 1;

	struct ru *victim_ru = NULL;
	struct ssdparams *spp = &ssd->sp;
	struct nand_lun *lunp; 
	struct ppa ppa;
	struct fdp_ru_mgmt *rum = &ssd->rums[rgid];
	NvmeRuHandle *ruh;
	NvmeFdpEvent *e = NULL;
	NvmeNamespace *ns = req->ns;
	int start_lunidx = rgid * RG_DEGREE;
	uint16_t ruhid;

	int gc_pgs = 0;

	victim_ru = select_victim_ru(ssd, force, rgid);
	if (!victim_ru) {
		return -1;
	}

    ppa.g.blk = victim_ru->id;
	ruhid = victim_ru->ruhid; 
	ruh = &ns->endgrp->fdp.ruhs[ruhid];	

    ftl_debug("GC-ing line:%d,ipc=%d,victim=%d,full=%d,free=%d\n", ppa.g.blk,
              victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
              ssd->lm.free_line_cnt); 

#ifdef FDP_DEBUG
	printf("rgid: %d\n", rgid);
	printf("ruhid: %d\n", ruhid);
	printf("victim_ru id: %d\n", victim_ru->id);
	printf("victim_ru->ipc: %d\n", victim_ru->ipc);
	printf("victim_ru->vpc: %d\n", victim_ru->vpc);
#endif 

	for (int lunidx = start_lunidx; lunidx < start_lunidx + RG_DEGREE; lunidx++) {
		ppa.g.ch = lunidx / spp->luns_per_ch;
		ppa.g.lun = lunidx % spp->luns_per_ch;
		ppa.g.pl = 0;
		lunp = get_lun(ssd, &ppa);
		gc_pgs += fdp_clean_one_block(ssd, &ppa, rgid, ruhid);	//update
		mark_block_free(ssd, &ppa);

		if (spp->enable_gc_delay)
		{
			struct nand_cmd gce;
			gce.type = GC_IO;
			gce.cmd = NAND_ERASE;
			gce.stime = 0;
			ssd_advance_status(ssd, &ppa, &gce);
		} 

		lunp->gc_endtime = lunp->next_lun_avail_time;
	} 

	if (ruh->ruht == NVME_RUHT_INITIALLY_ISOLATED && log_event(ruh, FDP_EVT_MEDIA_REALLOC)) {
		printf("here0\n");
		struct nvme_fdp_event_realloc mr;
		e = nvme_fdp_alloc_event(req->ns->ctrl, &ns->endgrp->fdp.ctrl_events);
		e->type = FDP_EVT_MEDIA_REALLOC;
		e->flags = FDPEF_PIV | FDPEF_NSIDV | FDPEF_LV;
		e->pid = cpu_to_le16(ruhid);
		e->nsid = cpu_to_le32(ns->id);
		mr.flags = 1 << 0; // LIV on
		mr.nlbam = gc_pgs * 8;
		mr.lba = 0; // TODO
		memcpy(e->type_specific, &mr, sizeof(mr));
		e->rgid = cpu_to_le16(rgid);
		e->ruhid = cpu_to_le16(ruhid);
	}

#ifdef WAF_TEST
	req->ns->ctrl->gc_writes += gc_pgs * 8;
#endif
#ifdef DEVICE_UTIL_DEBUG
	spp->tt_valid_pgs += gc_pgs;
#endif
	/* reset wp of victim ru */
	victim_ru->wp.ch = start_lunidx / spp->luns_per_ch;
	victim_ru->wp.lun = start_lunidx % spp->luns_per_ch;
	victim_ru->wp.pl = 0;
	victim_ru->wp.blk = victim_ru->id;
	victim_ru->wp.pg = 0;

    /* update ru status */
	victim_ru->ipc = 0;
	victim_ru->vpc = 0;
	QTAILQ_INSERT_TAIL(&rum->free_ru_list, victim_ru, entry);
	rum->free_ru_cnt++;

	return 0;
}


static uint64_t ssd_read(struct ssd *ssd, NvmeRequest *req)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t lba = req->slba;
    int nsecs = req->nlb;
    struct ppa ppa;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + nsecs - 1) / spp->secs_per_pg;
    uint64_t lpn;
    uint64_t sublat, maxlat = 0;

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    /* normal IO read path */
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);
        if (!mapped_ppa(&ppa) || !valid_ppa(ssd, &ppa)) {
            //printf("%s,lpn(%" PRId64 ") not mapped to valid ppa\n", ssd->ssdname, lpn);
            //printf("Invalid ppa,ch:%d,lun:%d,blk:%d,pl:%d,pg:%d,sec:%d\n",
            //ppa.g.ch, ppa.g.lun, ppa.g.blk, ppa.g.pl, ppa.g.pg, ppa.g.sec);
            continue;
        }

        struct nand_cmd srd;
        srd.type = USER_IO;
        srd.cmd = NAND_READ;
        srd.stime = req->stime;
        sublat = ssd_advance_status(ssd, &ppa, &srd);
        maxlat = (sublat > maxlat) ? sublat : maxlat;
    }

#ifdef UPDATE_FREQ
	if (lba == 1000) {
		for (int i = 0; i < NR_TENANTS; i++)
			for (int j = 0; j < LGROUPS_PER_TENANT; j++)
				printf("tenant: %d lgroup: %d cnt: %d\n", i, j, ssd->ten[i].update_cnt[j]);
	}
#endif
#ifdef DEVICE_UTIL_DEBUG
	if (lba == 2000) 
		printf("util: %lf\n", (double) spp->tt_valid_pgs / spp->tt_pgs);
#endif
    return maxlat; 
}

static uint64_t ssd_write(struct ssd *ssd, NvmeRequest *req)
{
	uint64_t lba = req->slba;
	struct ssdparams *spp = &ssd->sp;
    int len = req->nlb;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + len - 1) / spp->secs_per_pg;
    struct ppa ppa;
    uint64_t lpn;
    uint64_t curlat = 0, maxlat = 0;
    int r = 0;

	NvmeRwCmd *rw = (NvmeRwCmd*)&req->cmd;				//update~
	NvmeNamespace *ns = req->ns;
	NvmeEnduranceGroup *endgrp = ns->endgrp;
	bool fdp_enabled = endgrp->fdp.enabled;			
    uint32_t dw12 = le32_to_cpu(req->cmd.cdw12);
    uint8_t dtype = (dw12 >> 20) & 0xf;
	uint16_t pid = le16_to_cpu(rw->dspec);
	uint16_t rgif = endgrp->fdp.rgif;						
	uint16_t rgid = pid >> (16 - rgif);
	uint16_t ph = pid & ((1 << (15 - rgif)) - 1);
	uint16_t ruhid;										

	if (dtype != NVME_DIRECTIVE_DATA_PLACEMENT) {
		ph = 0;
		rgid = 0; // TODO: consider striping later
	}
	ruhid = ns->fdp.phs[ph];							//~update

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    } 

	if (fdp_enabled) {	//update~
		/* perform GC here until !should_fdp_gc(ssd, rgid) */
		while (should_fdp_gc_high(ssd, rgid)) {
#ifdef FDP_DEBUG
			printf("do_fdp_gc() called in high\n");
#endif
			spp->tt_valid_pgs -= spp->pgs_per_ru;
			r = do_fdp_gc(ssd, rgid, true, req);
			ftl_assert(spp->tt_valid_pgs >= 0);
#ifdef DEVICE_UTIL_DEBUG
#endif
			if (r == -1)
				break;
		}
	}					//~update
	else {
		/* perform GC here until !should_gc(ssd) */
		while (should_gc_high(ssd)) {
			r = do_gc(ssd, true);
			if (r == -1)
				break;
		}
	}

    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);
        if (mapped_ppa(&ppa)) {
            /* update old page information first */
#ifdef UPDATE_FREQ 
			int tid = lpn / LPNS_PER_TENANT; 					// tenant id
			int lgroup_offset = (lpn % (int)LPNS_PER_TENANT) / LPNS_PER_LGROUP;
			/*
			printf("lpn: %ld\n", lpn);
			printf("tid: %d\n", tid);
			printf("lgroup_offset: %d\n", lgroup_offset);*/
			ssd->ten[tid].update_cnt[lgroup_offset]++;
#endif
			uint16_t old_rgid = (ppa.g.ch * spp->luns_per_ch + ppa.g.lun) / RG_DEGREE;
			mark_page_invalid(ssd, &ppa, old_rgid); 
            set_rmap_ent(ssd, INVALID_LPN, &ppa);
			ssd->gc_cnt[lpn] = 0;
        }

        /* new write */
		ppa = (fdp_enabled ? fdp_get_new_page(ssd, rgid, 0, ruhid, false) : get_new_page(ssd));

#ifdef FDP_DEBUG
		printf("pid: %10d lpn: %10ld rgid: %5d ruhid: %5d ch: %5d, lun: %5d, blk: %5d, pg: %5d\n", 
				pid, lpn, rgid, ruhid, ppa.g.ch, ppa.g.lun, ppa.g.blk, ppa.g.pg);
#endif

        /* update maptbl */
        set_maptbl_ent(ssd, lpn, &ppa);
        /* update rmap */
        set_rmap_ent(ssd, lpn, &ppa);

        mark_page_valid(ssd, &ppa);
#ifdef DEVICE_UTIL_DEBUG
		spp->tt_valid_pgs += 1;
		ftl_assert(spp->tt_valid_pgs <= spp->tt_pgs);
#endif

        /* need to advance the write pointer here */
		if (fdp_enabled)  {
			ssd_advance_fdp_write_pointer(ssd, rgid, 0, ruhid, false);
		}
		else
			ssd_advance_write_pointer(ssd);

        struct nand_cmd swr;
        swr.type = USER_IO;
        swr.cmd = NAND_WRITE;
        swr.stime = req->stime;
        /* get latency statistics */
        curlat = ssd_advance_status(ssd, &ppa, &swr);
        maxlat = (curlat > maxlat) ? curlat : maxlat;
    }

    return maxlat;
}

static void *ftl_thread(void *arg)
{
    FemuCtrl *n = (FemuCtrl *)arg;
    struct ssd *ssd = n->ssd;
    NvmeRequest *req = NULL;
    uint64_t lat = 0;
    int rc;
    int i;

	//NvmeRwCmd *rw;			//update~
	//NvmeNamespace *ns;
	//NvmeEnduranceGroup *endgrp;
	//uint16_t pid;
	/*uint16_t rgif;
	bool fdp_enabled; */
	//uint16_t rgid;			//~update

    while (!*(ssd->dataplane_started_ptr)) {
        usleep(100000);
    }

    /* FIXME: not safe, to handle ->to_ftl and ->to_poller gracefully */
    ssd->to_ftl = n->to_ftl;
    ssd->to_poller = n->to_poller;

    while (1) {
        for (i = 1; i <= n->nr_pollers; i++) {
            if (!ssd->to_ftl[i] || !femu_ring_count(ssd->to_ftl[i]))
                continue;

            rc = femu_ring_dequeue(ssd->to_ftl[i], (void *)&req, 1);
            if (rc != 1) {
                printf("FEMU: FTL to_ftl dequeue failed\n");
			}

			//rw = (NvmeRwCmd*)&req->cmd;			//update~
			//ns = req->ns;
			//endgrp = ns->endgrp;
			//pid = le16_to_cpu(rw->dspec);
			/*rgif = endgrp->fdp.rgif;							//update~
			fdp_enabled = endgrp->fdp.enabled;			
			rgid = pid >> (16 - rgif); */

			ftl_assert(req);
			switch (req->cmd.opcode) {
				case NVME_CMD_WRITE:
                lat = ssd_write(ssd, req);
#ifdef WAF_TEST 
				n->host_writes += req->nlb;
#endif
                break;
            case NVME_CMD_READ:
                lat = ssd_read(ssd, req);
                break;
            case NVME_CMD_DSM:
                lat = 0;
                break;
            default:
                //ftl_err("FTL received unkown request type, ERROR\n");
                ;
            }

            req->reqlat = lat;
            req->expire_time += lat;

            rc = femu_ring_enqueue(ssd->to_poller[i], (void *)&req, 1);
            if (rc != 1) {
                ftl_err("FTL to_poller enqueue failed\n");
            }

            /* clean one ru if needed (in the background) */
			/*
			if (fdp_enabled) {
				if (should_fdp_gc(ssd, rgid)) {
#ifdef FDP_DEBUG
					printf("do_fdp_gc() called in normal\n");
#endif
					do_fdp_gc(ssd, rgid, false, req);
				} 
			}
			else {
				if (should_gc(ssd)) {
					do_gc(ssd, false);
				}
			}*/
        }
    }

    return NULL;
}
