/*================= defolt seting =================

  This file is read by all code.

  FUNCTION
  pick_base : Create basic unit parameters on all codes in input redshift.
  rtoz      : Convert distance to redshift
  buffer    : Skip any lines when loading a file.

==================================================*/
#ifndef _COMMON_H_
#define _COMMON_H_


extern char fx_data_ori[256];
extern char fx_data_pre[256];
extern char fx_deep_res[256];

extern int rand_array[1000];

// extern int N_WEIGHT;
// extern int WEIGHT_DAY;
extern int DEEP_DAY;
extern double WEIGHT[50];


void buffer(int k, FILE *s);
int file_number(char *file_name, int ignore_line);
double array_Max(double* arr, int init, int size);
double array_Min(double* arr, int init, int size);
double array_average(double* arr, int init, int size);
double array_standerd(double* arr, int init, int size, double average);
void create_rand_array(void);
void fft(double _Complex* input, double _Complex* output, int FNDATA, int FNPRE);
void ifft(double _Complex* input, double _Complex* output, int FNDATA, int FNPRE);
double RSS(double* dataA, double* dataB, int SNDATA);
double Sigmoid(double x);

#endif //_COMMON_H_
