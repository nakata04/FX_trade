/*=================      note      =================

  This file is read by all code.

  FUNCTION
  pick_base : Create basic unit parameters on all codes in input redshift.
  rtoz      : Convert distance to redshift
  buffer    : Skip any lines when loading a file.

  moodycat

==================================================*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "fft_fx.h"

int main (int argc, char **argv)
{
  /*===== read "parameter.txt" file =====*/
  load_parameter();

  NDATA = file_number(histrical_data, 1);
  printf("%d\n", NDATA);

  fft_accuracy();

  printf("done\n");

  return 0;
}
