#include "test.h"

int main()
{
	long int a;
	a |= (long int)1 << 32;
	printf("%ld\n", a);
	return 0;
}
