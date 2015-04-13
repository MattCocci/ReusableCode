#include<stdio.h>
#include "lapacke.h"

int main(int argc, char *argv[])
{
  double a[2*2] = {0,1,0,1};
  lapack_int lda, n, info;

  n = 2;
  lda = 2;


  info = LAPACKE_dpotri( LAPACK_COL_MAJOR, 'U', n, a, lda );

  return info;
}
