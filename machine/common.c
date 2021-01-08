/* ====================        note        ====================

  This file is read by all code.

  FUNCTION
  pick_base : Create basic unit parameters on all codes in input redshift.
  rtoz      : Convert distance to redshift
  buffer    : Skip any lines when loading a file.

============================================================ */

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <sys/types.h>
#include <unistd.h>
#include "common.h"


char fx_data_ori[256];
char fx_data_pre[256];
char fx_deep_res[256];

int rand_array[1000];

// int N_WEIGHT;
// int WEIGHT_DAY;
int DEEP_DAY;




/* ============================================================ */
/*                       buffer function                        */
/*  k lines skip in s file                                      */
/*                                                              */
/*  argument   k  : the number of skip lines                    */
/*             s  : file name                                   */
/* ============================================================ */

void buffer(int k, FILE *s)
{
  int i;
  char buf[256];
  for(i = 0; i < k; i++) fgets(buf, sizeof(buf), s);
}


/* ============================================================ */
/*                    count the file number                     */
/*                                                              */
/*  argument  f_name  : file name you want to count             */
/*            rm_line : number of removed line                  */
/*  return    line    : number of line                          */
/* ============================================================ */

int file_number(char *f_name, int rm_line)
{
  int c;
  int line = 0;
  FILE *fp;

  fp = fopen(f_name, "r");
  if (fp == NULL)
  {
    printf("Could not open file : %s ( common.c )\n", f_name);
    exit (1);
  }

  while ((c = fgetc(fp)) != EOF)
  {
    if (c == '\n') line++;
  }

  fclose(fp);

  line -= rm_line;

  return line;
}


/* ============================================================ */
/*                    find the Max value                        */
/*  One point value by spline function                          */
/*                                                              */
/*  argument   arr  : input array                               */
/*             init :                                           */
/*             size :                                           */
/*  return     buf  :                                           */
/* ============================================================ */

double array_Max(double* arr, int init, int size)
{
  int i;
  double buf = 0;

  for (i = init; i < init + size; i++)
  {
    if(arr[i] == 0)
    {
      printf("worning : array has arr[%d]=%lf value\n", i, arr[i]);
      exit (1);
    }
    if(buf < arr[i]) buf = arr[i];
  }

  return buf;
}


/* ============================================================ */
/*                    find the Min value                        */
/*                                                              */
/*  argument   arr  : input array                               */
/*             init :                                           */
/*             size :                                           */
/*  return     buf  :                                           */
/* ============================================================ */

double array_Min(double* arr, int init, int size)
{
  int i;
  double buf;

  buf = arr[0];
  for (i = init+1; i < init + size; i++)
  {
    if(arr[i] == 0)
    {
      printf("worning : array has arr[%d]=%lf value\n", i, arr[i]);
      exit (1);
    }
    if(buf > arr[i]) buf = arr[i];
  }

  return buf;
}


/* ============================================================ */
/*                     create axerage value                     */
/*                                                              */
/*  argument   arr  : input array                               */
/*             init :                                           */
/*             size :                                           */
/*  return     buf  :                                           */
/* ============================================================ */

double array_average(double* arr, int init, int size)
{
  int i;
  double buf = 0;

  // n = fin - init; + 1
  // buf = arr[0];
  for (i = init; i < init + size; i++)
  {
    if(arr[i] == 0)
    {
      printf("worning : array has arr[%d]=%lf value\n", i, arr[i]);
      exit (1);
    }

    buf += arr[i];
    // if(buf > arr[i]) buf = arr[i];
  }

  return buf / size;
}


/* ============================================================ */
/*                     create axerage value                     */
/*                                                              */
/*  argument   arr  : input array                               */
/*             init :                                           */
/*             size :                                           */
/*  return     buf  :                                           */
/* ============================================================ */

double array_standerd(double* arr, int init, int size, double average)
{
  int i;
  double buf = 0;

  // n = fin - init; + 1
  // buf = arr[0];
  for (i = init; i < init + size; i++)
  {
    if(arr[i] == 0)
    {
      printf("worning : array has arr[%d]=%lf value\n", i, arr[i]);
      exit (1);
    }
    buf += pow(arr[i] - average, 2);
    // if(buf > arr[i]) buf = arr[i];
  }
  buf = sqrt(buf / size);
  // buf /= size;

  return buf;
}



/* ============================================================ */
/*                  create random number array                  */
/*                                                              */
/* ============================================================ */

void create_rand_array(void)
{
  int i;

  for (i = 0; i < DEEP_DAY; i++)
  {
    // rand_array[i] = rand() % (PAST - 200) + 1;
    rand_array[i] = rand() % (400) + 1;
    // printf("%d ", rand_array[i]);
  }
}


/* ============================================================ */
/*                 fast Fourier transform (fft)                 */
/*  This function carry out the fft.                            */
/*                                                              */
/*  argument   input  : original data array                     */
/*             output : output data array after fft             */
/*             FNDATA : the number of data                      */
/*             FNPRE  : the number of future data               */
/* ============================================================ */

void fft(double _Complex* input, double _Complex* output, int FNDATA, int FNPRE)
{
  int m, n;
  double val_real, val_imag;

  for (m = 0; m < FNDATA; m++)
  {
    val_real = 0;
    val_imag = 0;
    for (n = 0; n < FNDATA; n++)
    {
      val_real += creal(input[n]) * cos(2 * n * m * M_PI / FNDATA);
      val_imag -= creal(input[n]) * sin(2 * n * m * M_PI / FNDATA);
    }
    output[m] = val_real + I * val_imag;
  }
}


/* ============================================================ */
/*            inverse fast Fourier transform (ifft)             */
/*  This function carry out the ifft.                           */
/*                                                              */
/*  argument   input  : original data array                     */
/*             output : output data array after ifft            */
/*             FNDATA : the number of data                      */
/*             FNPRE  : the number of future data               */
/* ============================================================ */

void ifft(double _Complex* input, double _Complex* output, int FNDATA, int FNPRE)
{
  int m, n;
  double val_real, val_imag;

  for (m = 0; m < FNDATA + FNPRE; m++)
  {
    val_real = 0;
    val_imag = 0;
    // for (n = 0; n < FNDATA; n++)
    for (n = 0; n < FNDATA; n++)
    {
      val_real += creal(input[n]) * cos(2 * n * m * M_PI / FNDATA) - cimag(input[n]) * sin(2 * n * m * M_PI / FNDATA);
      val_imag += creal(input[n]) * sin(2 * n * m * M_PI / FNDATA) + cimag(input[n]) * cos(2 * n * m * M_PI / FNDATA);
    }
    output[m] = (val_real + I * val_imag) / FNDATA;
  }
}



/* ============================================================ */
/*                residual sum of squares : RSS                 */
/*                                                              */
/*  argument   x :                                              */
/*  return     h :                                              */
/* ============================================================ */

double RSS(double* dataA, double* dataB, int SNDATA)
{
  int n;
  double rss = 0;

  for (n = 0; n < SNDATA; n++) rss += pow(dataA[n] - dataB[n], 2);

  return rss;
}


/* ============================================================ */
/*                       Sigmoid function                       */
/*  One point value by spline function                          */
/*                                                              */
/*  argument   x :                                              */
/*  return     h :                                              */
/* ============================================================ */

double Sigmoid(double x)
{
  double h;

  h = pow(1 + exp(-x), -1);

  return h;
}
