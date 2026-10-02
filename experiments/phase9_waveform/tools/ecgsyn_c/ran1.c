/* Stand-in for Numerical Recipes ran1 (uniform deviate in (0,1)); used only by the RR
   process (irrelevant when hrstd = 0) and the additive noise (Anoise = 0 here). */
#include <stdlib.h>
float ran1(long *idum) { static int init=0; if(!init){ srand48((long)(-*idum)); init=1; } return (float)drand48(); }
