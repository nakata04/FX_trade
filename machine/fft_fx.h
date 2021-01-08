/*================= defolt seting =================

  This file is read by all code.

  FUNCTION
  pick_base : Create basic unit parameters on all codes in input redshift.
  rtoz      : Convert distance to redshift
  buffer    : Skip any lines when loading a file.

==================================================*/
#ifndef _FFT_FX_H_
#define _FFT_FX_H_


extern int RUN;
extern int NDATA;
extern int NEW_ML_CONF;
extern int RESTART;
extern int DATA_TYPE;
extern int DATA_TYPE;
extern int TARGET_TIME;
extern int PAST;
extern int ML_TRIAL;
extern int ML_BIN;
extern int TRIAL;
extern int PIPS_THRESHOLD;
extern int FUTURE;
extern int TIME_ADJUST;
extern int LEARNING_BIN;
extern int MACHINE_DAY;
extern int N_WEIGHT;
extern int SAMPLE_LINE;

extern double LEARNING_RATE;
extern double ML_THRESHOLD;
extern double RELIABILITY;

extern char histrical_data[512];
extern char fft_data[512];
extern char fft_accuracy_data[512];
extern char param_conf[512];

#define TEST_DAY 1000

void load_parameter();
void data_set(char** timec, double _Complex* data_input, double* Mclose_log, double* data_dif, double* data_buf, double* data_Max, double* data_Min);
void fft_fftw3();
void data_set_accuracy(double* data_dif, double ***data_machine, double* Mclose_log, int MD, int SAMPLE_DAYS);
// void data_set_accuracy(double _Complex* data_input, double* Mclose_log, double* data_dif, double* data_buf, double* data_Max, double* data_Min, int SAMPLE_DAYS);
// void machine_array(double *data_buf, double* data_Max, double* data_Min, double pre_ifft, double ***data_machine, int N_LINE);
void fft_accuracy();

void machine_learning(double *data_dif, double ***data_machine, double **wconf, double *INTERSEPT);
double M_f(double *data_dif, double ***data_machine, double **wconf, double *INTERSEPT);
void M_df(double *data_dif, double ***data_machine, double **wconf, double *INTERSEPT);

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])


#endif //_FFT_FX_H_
