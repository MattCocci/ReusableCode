/***********************************************************************
 * LU.c -- LU Decomposition of a Matrix
 *
 * Author: Matthew D. Cocci
 *
 * Adapted from Numerical Recipes in C
 * ********************************************************************/

#include<math.h>
#include "utils.h"
#define TINY 1.0e-20

void ludcmp(float **a, int n, int *indx, float *d)
{
  int i, imax, j, k;
  float big, dum, sum, temp;
  float *vv;

  vv=vector(1,n);
  *d = 1.0;
}

int main(int argc, char *argv[])
{
  return 0;
}
