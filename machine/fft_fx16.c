/*================= defolt seting =================

  INCL = common.h figure.h
  MAIN = read_snapshot_data.o common.o figure9.o main_figure3.o
  ALL  = figure

  Analyze data for the last 'PAST' days.
  https://www.central-tanshifx.com/market/finder/popn-csv-download.html
  https://sec.himawari-group.co.jp/report/chart/

==================================================*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <dirent.h>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_fft_complex.h>
#include <complex.h> /* This liblary has to include before "fftw3.h" */
// #include <fftw3.h>
#include "common.h"
#include "fft_fx.h"


int RUN;
int NDATA;
int NEW_ML_CONF;
int RESTART;
int BIN;
int DATA_TYPE;
int TARGET_TIME;
int PAST;
int TRIAL;
int ML_TRIAL;
int ML_BIN;
int PIPS_THRESHOLD;
int FUTURE;
int TIME_ADJUST;
int LEARNING_BIN;
int MACHINE_DAY;
int N_WEIGHT;
int SAMPLE_LINE;

double LEARNING_RATE;
double ML_THRESHOLD;
double RELIABILITY;

char histrical_data[512];
char fft_data[512];
char fft_accuracy_data[512];
char param_conf[512];


void load_parameter()
{
  char FX_PAIR[128], DATA_BIN[128], DATA_FORMAT[128];
  char base[64], data_dir[128], log_dir[128];
  struct stat buf;

  /*===== read parameter.txt =====*/
  FILE *fp;
  char line[128];
  char name[128];
  char dummy[128];
  char value[128];

  fp = fopen("./parameter.txt", "r");
  while(fgets(line, 256, fp) != NULL)
  {
    /*===== Skip lines starting with "#" or line breaks =====*/
    if((line[0] == '#') || (!strcmp(line, "\n")));
    else{
      sscanf(line, "%s %s %s\n", name, dummy, value);
      if(!strcmp(name, "FX_PAIR"))
        strcpy(FX_PAIR, value);
      else if(!strcmp(name, "DATA_BIN"))
        strcpy(DATA_BIN, value);
      else if(!strcmp(name, "DATA_FORMAT"))
        strcpy(DATA_FORMAT, value);
      else if(!strcmp(name, "base_dir"))
        strcpy(base, value);
      else if(!strcmp(name, "data_dir"))
        strcpy(data_dir, value);
      else if(!strcmp(name, "log_dir"))
        strcpy(log_dir, value);
      else if(!strcmp(name, "RUN"))
        RUN            = atoi(value);
      else if(!strcmp(name, "NEW_ML_CONF"))
        NEW_ML_CONF    = atoi(value);
      else if(!strcmp(name, "RESTART"))
        RESTART        = atoi(value);
      else if(!strcmp(name, "BIN"))
        BIN            = atoi(value);
      else if(!strcmp(name, "PAST"))
        PAST           = atoi(value);
      else if(!strcmp(name, "TARGET_TIME"))
        TARGET_TIME    = atoi(value);
      else if(!strcmp(name, "TRIAL"))
        TRIAL          = atoi(value);
      else if(!strcmp(name, "PIPS_THRESHOLD"))
        PIPS_THRESHOLD = atoi(value);
      else if(!strcmp(name, "FUTURE"))
        FUTURE         = atoi(value);
      else if(!strcmp(name, "TIME_ADJUST"))
        TIME_ADJUST    = atoi(value);
      else if(!strcmp(name, "N_WEIGHT"))
        N_WEIGHT       = atoi(value);
      else if(!strcmp(name, "LEARNING_BIN"))
        LEARNING_BIN   = atoi(value);
      else if(!strcmp(name, "MACHINE_DAY"))
        MACHINE_DAY    = atoi(value);
      else if(!strcmp(name, "SAMPLE_LINE"))
        SAMPLE_LINE    = atoi(value);
      else if(!strcmp(name, "ML_TRIAL"))
        ML_TRIAL       = atoi(value);
      else if(!strcmp(name, "ML_BIN"))
        ML_BIN         = atoi(value);
      else if(!strcmp(name, "LEARNING_RATE"))
        LEARNING_RATE  = atof(value);
      else if(!strcmp(name, "ML_THRESHOLD"))
        ML_THRESHOLD   = atof(value);
      else if(!strcmp(name, "RELIABILITY"))
        RELIABILITY    = atof(value);
    }
  }
  fclose(fp);
  // printf("%lf\n", LEARNING_RATE);

  sprintf(fft_data, "%s/%s_fft.txt", data_dir, FX_PAIR);
  sprintf(histrical_data, "%s/%s_%s.%s", base, FX_PAIR, DATA_BIN, DATA_FORMAT);
  sprintf(fft_accuracy_data, "%s/%s_acc.txt", data_dir, FX_PAIR);
  sprintf(param_conf, "%s/conf_%s_B%dN%dL%dM%d.txt", data_dir, FX_PAIR, BIN, N_WEIGHT, LEARNING_BIN, MACHINE_DAY);
  if(mkdir(data_dir, 0777) == 0) stat(data_dir, &buf);

  if(!strcmp(DATA_BIN, "M1"))
  {
    DATA_TYPE = 1;
    ML_BIN    = 127;
  }
  else if(!strcmp(DATA_BIN, "H1"))
  {
    DATA_TYPE = 60;
    ML_BIN    = 31;
  }
  else if(!strcmp(DATA_BIN, "D"))
  {
    DATA_TYPE = 1440;
    ML_BIN    = 3;
  }
  else
  {
    printf("DATA TYPE error : %s\n", DATA_BIN);
    exit (1);
  }
}


void data_set(char** timec, double _Complex* data_input, double* Mclose_log, double* data_dif, double* data_buf, double* data_Max, double* data_Min)
{
  int i, j, count_bin, count_data;
  char datac[32];
  double Mopen, Mhigh, Mlow, Mclose, Mclose_buf, Mclose_ave, Mclose_ave_buf, average25;

  FILE *fp;

  fp = fopen(histrical_data, "r");
  if(fp == NULL)
  {
    printf("Could not open file : %s (  )\n", histrical_data);
    exit (1);
  }

  buffer(1, fp);
  if(NDATA - ((PAST + LEARNING_BIN + 25) * BIN) < 0)
  {
    printf("setting error\n");
    exit (1);
  }
  buffer(NDATA - ((PAST + LEARNING_BIN + 25) * BIN), fp);

  fscanf(fp, "\n%[^,],%lf,%lf,%lf,%lf", datac, &Mopen, &Mhigh, &Mlow, &Mclose);
  Mclose_log[0] = Mclose;
  Mclose_log[1] = Mclose;
  Mclose_buf    = Mclose;
  Mclose_ave_buf = Mclose;

  double *data_average_buf;
  data_average_buf = (double*)malloc(sizeof(double) * 25);

  count_bin  = 1;
  count_data = 0;
  Mclose_ave = 0;
  for (i = 0; i < (PAST + LEARNING_BIN + MACHINE_DAY + 25) * BIN; i++) // worning : check number
  {
    if(i % BIN == BIN - 1)
    {
      fscanf(fp, "\n%[^,],%lf,%lf,%lf,%lf", datac, &Mopen, &Mhigh, &Mlow, &Mclose);

      if(count_bin < 25-1)
      {
        data_average_buf[count_bin] = Mclose;
      }
      count_bin++;

      if(count_bin >= 25-1)
      {
        if(count_bin == 25-1)
        {
          data_average_buf[count_bin] = Mclose;
          Mclose_log[0]  = Mclose;
          Mclose_log[1]  = Mclose;
          Mclose_buf     = Mclose;
          average25 = array_average(data_average_buf, 0, 25);
          data_buf[count_data] = Mclose - average25;
          data_Max[count_data] = Mhigh  - average25;
          data_Min[count_data] = Mlow   - average25;
          count_data++;
        }
        else
        {
          for (j = 0; j < 25 - 1; j++)
          {
            data_average_buf[j] = data_average_buf[j+1];
          }
          data_average_buf[25-1] = Mclose;
          average25 = array_average(data_average_buf, 0, 25);
          data_dif[count_data] = Mclose - Mclose_buf;
          data_buf[count_data] = Mclose - average25;
          data_Max[count_data] = Mhigh  - average25;
          data_Min[count_data] = Mlow   - average25;
          Mclose_buf = Mclose;

          if(i < PAST * BIN)
          {
            data_input[count_data] = ((Mclose_ave + Mclose) / (BIN)) - Mclose_ave_buf;
            Mclose_ave_buf = (Mclose_ave + Mclose) / (BIN);
            Mclose_ave = 0;
          }

          count_data++;
        }
      }
    }
    else
    {
      // fscanf(fp, "\n%[^,],%*lf,%*lf,%*lf,%lf", datac, &Mclose);
      fscanf(fp, "\n%[^,],%lf,%lf,%lf,%lf", datac, &Mopen, &Mhigh, &Mlow, &Mclose);
      Mclose_ave += Mclose;
    }
  }

  fclose(fp);
}


void fft_fftw3()
{
  int i, j, wd, nw, md;
  double wx, Mclose_log[2], pips=0;
  int accuracy=0, accuracy_tot=0, accuracy_threshold=0, accuracy_threshold_tot=0;
  FILE *fp2;

  double *data_dif;
  data_dif   = (double*)malloc(sizeof(double) * (FUTURE));
  double *E_wb;
  E_wb   = (double*)malloc(sizeof(double) * (FUTURE));

  double ***data_machine = (double***)malloc(sizeof(double**) * LEARNING_BIN);
  if(data_machine == NULL){
    printf("d3arrayメモリが確保できません\n");
    exit(1);
  }
  data_machine[0] = (double**)malloc(sizeof(double*) * LEARNING_BIN * 13);
  if(data_machine[0] == NULL){
    printf("d3array[0]メモリが確保できません\n");
    exit(1);
  }
  data_machine[0][0] = (double*)malloc(sizeof(double) * LEARNING_BIN * 13 * FUTURE);
  if(data_machine[0][0] == NULL){
    printf("d3array[0][0]=%d %d %dメモリが確保できません\n", LEARNING_BIN, 13, FUTURE);
    exit(1);
  }
  for(i = 0; i < LEARNING_BIN; i++)
  {
    data_machine[i] = data_machine[0] + i * 13;
    for (j = 0; j < 13; j++)
    {
      data_machine[i][j] = data_machine[0][0] + i * (13 * FUTURE) + j * FUTURE;
    }
  }

  /*===== initialized parameters =====*/
  double **wconf;
  wconf = malloc(sizeof(double *) * LEARNING_BIN);
  for(i = 0; i < LEARNING_BIN; i++)
  {
    wconf[i] = malloc(sizeof(double) * N_WEIGHT);
  }
  double *INTERSEPT;
  INTERSEPT = malloc(sizeof(double *) * 2);
  INTERSEPT[0] = 1;

  // fp2 = fopen("./param_conf.txt", "r");
  fp2 = fopen(param_conf, "r");
  fscanf(fp2, "%lf", &INTERSEPT[0]);
  for (i = 0; i < LEARNING_BIN; i++)
  {
    /*=====        1   2   3   4   5   6   7   8   9  10  11  12  13 =====*/
    fscanf(fp2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      &wconf[i][0],  &wconf[i][1],  &wconf[i][2], &wconf[i][3], &wconf[i][4],
      &wconf[i][5],  &wconf[i][6],  &wconf[i][7], &wconf[i][8], &wconf[i][9],
      &wconf[i][10], &wconf[i][11], &wconf[i][12]);
  }
  fclose(fp2);

  double *Y_out;
  Y_out = (double*)malloc(sizeof(double) * FUTURE);


  for (md = 0; md < FUTURE; md++)
  {
    data_set_accuracy(data_dif, data_machine, Mclose_log, md, FUTURE-md);
    wx = 0;
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        wx += data_machine[wd][nw][md] * wconf[wd][nw];
      }
    }
    E_wb[md] = fabs(data_dif[md] - (wx + INTERSEPT[0])) / (N_WEIGHT * LEARNING_BIN);
  }

  for (md = 0; md < FUTURE; md++)
  {
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < 13; nw++)
      {
        if(nw == 12 && wd == LEARNING_BIN - 1)
        {
          printf("%e\n\n", data_machine[wd][nw][md]);
        }
        else if(nw == 12)
        {
          printf("%e\n", data_machine[wd][nw][md]);
        }
        else
        {
          printf("%e ", data_machine[wd][nw][md]);
        }
      }
    }
  }

  for (md = 0; md < FUTURE; md++)
  {
    Y_out[md] = 0;
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        Y_out[md] += data_machine[wd][nw][md] * wconf[wd][nw];
      }
    }
    Y_out[md] += INTERSEPT[0];

    Mclose_log[0] += data_dif[md];
    printf("%d %e dif : %e pips : %e reliability : %e\n", md+1, Mclose_log[0], data_dif[md], Y_out[md]*100, E_wb[md]);
  }

  free(data_dif);
  free(wconf);
  free(INTERSEPT);
  free(data_machine);
  free(Y_out);
  free(E_wb);

}


void data_set_accuracy(double* data_dif, double ***data_machine, double* Mclose_log, int MD, int SAMPLE_DAYS)
{
  int i, lb, nw, count_data=0;
  char datac[32];
  double Mopen, Mhigh, Mlow, Mclose, Mhigh_buf=0, Mlow_buf=1.0e8, Mclose_ave, average25, standard25;

  FILE *fp;

  double *data_buf, *data_Max, *data_Min;
  data_buf   = (double*)malloc(sizeof(double) * (PAST + LEARNING_BIN + 1));
  data_Max   = (double*)malloc(sizeof(double) * (PAST + LEARNING_BIN + 1));
  data_Min   = (double*)malloc(sizeof(double) * (PAST + LEARNING_BIN + 1));

  double average_sum25=0;

  fp = fopen(histrical_data, "r");
  if(fp == NULL)
  {
    printf("Could not open file : %s (  )\n", histrical_data);
    exit (1);
  }

  buffer(1, fp);
  if(NDATA - ((PAST + LEARNING_BIN + 1) * BIN + (SAMPLE_DAYS + 1)) < 0)
  {
    printf("setting error\n");
    exit (1);
  }
  if(RUN == 1) buffer(NDATA - ((PAST + LEARNING_BIN) * BIN + SAMPLE_DAYS), fp);
  else buffer(NDATA - ((PAST + LEARNING_BIN + 1) * BIN + (SAMPLE_DAYS + 1)), fp);

  fscanf(fp, "\n%[^,],%lf,%lf,%lf,%lf", datac, &Mopen, &Mhigh, &Mlow, &Mclose);

  for (i = 0; i < (PAST + LEARNING_BIN + 1) * BIN; i++)
  {
    if(i % BIN == BIN - 1)
    {
      fscanf(fp, "\n%[^,],%lf,%lf,%lf,%lf", datac, &Mopen, &Mhigh, &Mlow, &Mclose);
      if(Mhigh > Mhigh_buf) Mhigh_buf = Mhigh;
      if(Mlow  < Mlow_buf)  Mlow_buf  = Mlow;

      data_buf[count_data] = Mclose;
      data_Max[count_data] = Mhigh_buf;
      data_Min[count_data] = Mlow_buf;

      Mhigh_buf = 0;
      Mlow_buf  = 1.0e8;

      count_data++;
    }
    else
    {
      fscanf(fp, "\n%[^,],%lf,%lf,%lf,%lf", datac, &Mopen, &Mhigh, &Mlow, &Mclose);
      if(Mhigh > Mhigh_buf) Mhigh_buf = Mhigh;
      if(Mlow  < Mlow_buf)  Mlow_buf  = Mlow;
      Mclose_ave += Mclose;
    }
  }
  fclose(fp);

  for (lb = 0; lb < LEARNING_BIN; lb++)
  {
    average25 = array_average(data_buf, PAST + lb - (25+1), 25);

    nw = 0;
    data_machine[lb][nw][MD]   = data_buf[PAST+lb];
    nw += 1;
    /*===== mean average =====*/
    data_machine[lb][nw][MD]   = array_average(data_buf, PAST+lb-5, 5);
    data_machine[lb][nw+1][MD] = array_average(data_buf, PAST+lb-25, 25);
    data_machine[lb][nw+2][MD] = array_average(data_buf, PAST+lb-75, 75);
    data_machine[lb][nw+3][MD] = array_average(data_buf, PAST+lb-200, 200);
    nw += 4;
    /*===== trans =====*/
    data_machine[lb][nw][MD]   = (array_Max(data_Max, PAST+lb-9, 9)      + array_Min(data_Min, PAST+lb-9, 9)) / 2;
    data_machine[lb][nw+1][MD] = (array_Max(data_Max, PAST+lb-26, 26)    + array_Min(data_Min, PAST+lb-26, 26)) / 2;
    data_machine[lb][nw+2][MD] = (array_Max(data_Max, PAST+lb-9-26, 9)   + array_Min(data_Min, PAST+lb-9-26, 9)
                               +  array_Max(data_Max, PAST+lb-26-26, 26) + array_Min(data_Min, PAST+lb-26-26, 26)) / 4;
    data_machine[lb][nw+3][MD] = (array_Max(data_Max, PAST+lb-52-26, 52) + array_Min(data_Min, PAST+lb-52-26, 52)) / 2;
    nw += 4;
    /*===== standard25 =====*/
    standard25                 = array_standerd(data_buf, PAST+lb-25, 25, average25);
    data_machine[lb][nw][MD]   = average25 + 1.0 * standard25;
    data_machine[lb][nw+1][MD] = average25 - 1.0 * standard25;
    data_machine[lb][nw+2][MD] = average25 + 2.0 * standard25;
    data_machine[lb][nw+3][MD] = average25 - 2.0 * standard25;
    nw += 4;

    if(lb == 0 && MD == 0)
    {
      printf("%d parameter\n", nw);
    }

    for (nw = 0; nw < 13; nw++)
    {
      average_sum25 += data_machine[lb][nw][MD];
    }
  }

  average_sum25 /= N_WEIGHT * LEARNING_BIN;

  for (lb = 0; lb < LEARNING_BIN; lb++)
  {
    for (nw = 0; nw < N_WEIGHT; nw++)
    {
      data_machine[lb][nw][MD] -= average_sum25;
    }
  }

  if(MD == 0)
  {
    Mclose_log[0] = data_buf[PAST+LEARNING_BIN-1];
    // printf("ml : %lf\n", Mclose_log[0]);
  }
  data_dif[MD]  = data_buf[PAST+LEARNING_BIN] - data_buf[PAST+LEARNING_BIN-1];
  // printf("%d data : %lf %lf dif : %lf\n", PAST+LEARNING_BIN, data_buf[PAST+LEARNING_BIN], data_buf[PAST+LEARNING_BIN-1], data_dif[MD]);

  free(data_buf);
  free(data_Max);
  free(data_Min);
}


void fft_accuracy()
{
  int i, j, wd, nw, md;
  double wx, Mclose_log[2], pips=0;
  int accuracy=0, accuracy_tot=0, accuracy_threshold=0, accuracy_threshold_tot=0;
  FILE *fp2, *fp3;
  double *data_dif;
  data_dif   = (double*)malloc(sizeof(double) * (MACHINE_DAY));
  double *E_wb;
  E_wb   = (double*)malloc(sizeof(double) * (MACHINE_DAY));

  double ***data_machine = (double***)malloc(sizeof(double**) * LEARNING_BIN);
  if(data_machine == NULL){
    printf("d3arrayメモリが確保できません\n");
    exit(1);
  }
  data_machine[0] = (double**)malloc(sizeof(double*) * LEARNING_BIN * 13);
  if(data_machine[0] == NULL){
    printf("d3array[0]メモリが確保できません\n");
    exit(1);
  }
  data_machine[0][0] = (double*)malloc(sizeof(double) * LEARNING_BIN * 13 * MACHINE_DAY);
  if(data_machine[0][0] == NULL){
    printf("d3array[0][0]=%d %d %dメモリが確保できません\n", LEARNING_BIN, 13, MACHINE_DAY);
    exit(1);
  }
  for(i = 0; i < LEARNING_BIN; i++)
  {
    data_machine[i] = data_machine[0] + i * 13;
    for (j = 0; j < 13; j++)
    {
      data_machine[i][j] = data_machine[0][0] + i * (13 * MACHINE_DAY) + j * MACHINE_DAY;
    }
  }

  /*===== initialized parameters =====*/
  double **wconf;
  wconf = malloc(sizeof(double *) * LEARNING_BIN);
  for(i = 0; i < LEARNING_BIN; i++)
  {
    wconf[i] = malloc(sizeof(double) * 13);
  }
  double *INTERSEPT;
  INTERSEPT = malloc(sizeof(double *) * 2);
  INTERSEPT[0] = 1;

  // wconf = (double *)malloc(nbin * 13 * sizeof(double));
  if(NEW_ML_CONF == 0 || RESTART == 1 || RUN == 1)
  {
    fp2 = fopen(param_conf, "r");
    fscanf(fp2, "%lf", &INTERSEPT[0]);
    for (i = 0; i < LEARNING_BIN; i++)
    {
      /*=====        1   2   3   4   5   6   7   8   9  10  11  12  13 =====*/
      fscanf(fp2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &wconf[i][0],  &wconf[i][1],  &wconf[i][2], &wconf[i][3], &wconf[i][4],
        &wconf[i][5],  &wconf[i][6],  &wconf[i][7], &wconf[i][8], &wconf[i][9],
        &wconf[i][10], &wconf[i][11], &wconf[i][12]);
    }
    fclose(fp2);
  }
  else
  {
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        wconf[wd][nw] = 0;
      }
    }
  }

  double *Y_out;
  Y_out = (double*)malloc(sizeof(double) * MACHINE_DAY);

  // int p;
  if(NEW_ML_CONF == 1 && RUN == 0)
  {
    for (md = 0; md < MACHINE_DAY; md++)
    {
      data_set_accuracy(data_dif, data_machine, Mclose_log, md, ML_BIN*md+1);
    }

    printf("start machine learning\n");
    machine_learning(data_dif, data_machine, wconf, INTERSEPT);

    for (md = 0; md < MACHINE_DAY; md++)
    {
      Y_out[md] = 0;
      for (wd = 0; wd < LEARNING_BIN; wd++)
      {
        for (nw = 0; nw < N_WEIGHT; nw++)
        {
          Y_out[md] += data_machine[wd][nw][md] * wconf[wd][nw];
        }
      }
      Y_out[md] += INTERSEPT[0];
    }
  }


  for (md = 0; md < MACHINE_DAY; md++)
  {
    data_set_accuracy(data_dif, data_machine, Mclose_log, md, ML_BIN*(MACHINE_DAY-(md+1))+SAMPLE_LINE);
    wx = 0;
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        wx += data_machine[wd][nw][md] * wconf[wd][nw];
      }
    }
    E_wb[md] = fabs(data_dif[md] - (wx + INTERSEPT[0])) / (N_WEIGHT * LEARNING_BIN);
  }

  md = MACHINE_DAY - 1;
  for (wd = 0; wd < LEARNING_BIN; wd++)
  {
    for (nw = 0; nw < 13; nw++)
    {
      if(nw == 12 && wd == LEARNING_BIN - 1)
      {
        printf("%e\n\n", data_machine[wd][nw][md]);
      }
      else if(nw == 12)
      {
        printf("%e\n", data_machine[wd][nw][md]);
      }
      else
      {
        printf("%e ", data_machine[wd][nw][md]);
      }
    }
  }

  fp3 = fopen(fft_accuracy_data, "w");
  if(fp3 == NULL)
  {
    printf("Could not open file : %s (  )\n", fft_accuracy_data);
    exit (1);
  }

  for (md = 0; md < MACHINE_DAY; md++)
  {
    Y_out[md] = 0;
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        Y_out[md] += data_machine[wd][nw][md] * wconf[wd][nw];
      }
    }
    Y_out[md] += INTERSEPT[0];

    if(Y_out[md] > PIPS_THRESHOLD * 0.01)
    {
      if(E_wb[md] < RELIABILITY)
      {
        accuracy_threshold_tot++;
        if(data_dif[md] > 0)
        {
          accuracy_threshold++;
        }
      }
      if(data_dif[md] > 0)
      {
        accuracy++;
      }
      accuracy_tot++;
      pips += (data_dif[md]) * 100;
    }
    else if(Y_out[md] < -PIPS_THRESHOLD * 0.01)
    {
      if(E_wb[md] < RELIABILITY)
      {
        accuracy_threshold_tot++;
        if(data_dif[md] < 0)
        {
          accuracy_threshold++;
        }
      }
      if(data_dif[md] < 0)
      {
        accuracy++;
      }
      accuracy_tot++;
      pips -= (data_dif[md]) * 100;
    }
    else
    {
      if(E_wb[md] < RELIABILITY)
      {
        accuracy_threshold_tot++;
      }
    }
    Mclose_log[0] += data_dif[md];
    fprintf(fp3, "%d %e %e %e %e\n", md+1, Mclose_log[0], pips, data_dif[md], Y_out[md]);

  }

  fclose(fp3);

  printf("accuracy : %lf (%d/%d)\n", (double)accuracy/accuracy_tot*100, accuracy, accuracy_tot);
  printf("Degree of reliability : %lf (%d/%d)\n",
    (double)accuracy_threshold/accuracy_threshold_tot*100, accuracy_threshold, accuracy_threshold_tot);

  free(data_dif);
  free(wconf);
  free(INTERSEPT);
  free(data_machine);
  free(Y_out);
  free(E_wb);
}


/* ============================================================ */
/*                       machine learning                       */
/*   data[input BIN][N parameter]                               */
/* ============================================================ */

void machine_learning(double* data_dif, double ***data_machine, double **wconf, double *INTERSEPT)
{
  int i, wd, nw, iter=0;
  double E_wb, E_wb_buf=1.0e8;
  FILE *fp;

  printf("%5d      ", iter);
  for (i = 0; i < N_WEIGHT; i++) printf("%.5f ", wconf[0][i]);
  printf("%.5f\n", INTERSEPT[0]);

  for (iter = 0; iter < ML_TRIAL; iter++)
  {
    E_wb = M_f(data_dif, data_machine, wconf, INTERSEPT);
    if(E_wb < ML_THRESHOLD)
    {
      printf("machine learning is completed!\n");
      printf("accuracy E_wb = %e\n", E_wb);
      break;
    }
    if(E_wb > 1.0e10)
    {
      printf("error : 'LEARNING_RATE' too learge!\n");
      exit(1);
    }

    M_df(data_dif, data_machine, wconf, INTERSEPT);

    if(iter % 1000 == 0)
    {
      printf("%5d %e ", iter+1, E_wb);
      for (i = 0; i < N_WEIGHT; i++) printf("%.5f ", wconf[0][i]);
      printf("%.5f\n", INTERSEPT[0]);
      if(E_wb > E_wb_buf)
      {
        printf("worning : 'LEARNING_RATE' too learge!\n");
        break;
      }
      E_wb_buf = E_wb;
    }
  }

  fp = fopen(param_conf, "w");
  fprintf(fp, "%e\n", INTERSEPT[0]);
  for (wd = 0; wd < LEARNING_BIN; wd++)
  {
    for (nw = 0; nw < 13; nw++)
    {
      if(nw == 12)
      {
        printf("%e\n", wconf[wd][nw]);
        fprintf(fp, "%e\n", wconf[wd][nw]);
      }
      else
      {
        printf("%e ", wconf[wd][nw]);
        fprintf(fp, "%e ", wconf[wd][nw]);
      }
    }
  }
  fclose(fp);

}

double M_f(double *data_dif, double ***data_machine, double **wconf, double *INTERSEPT)
{
  int nw, wd, md;
  double wx, E_wb = 0;

  for (md = 0; md < MACHINE_DAY; md++)
  {
    wx = 0;
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        wx += data_machine[wd][nw][md] * wconf[wd][nw];
      }
    }
    E_wb += pow(data_dif[md] - (wx + INTERSEPT[0]), 2);
  }

  return E_wb;
}

void M_df(double *data_dif, double ***data_machine, double **wconf, double *INTERSEPT)
{
  int i, nw, wd, nw2, wd2, md;
  double wx, dwx, db;

  double **wconf_buf;
  wconf_buf = malloc(sizeof(double *) * LEARNING_BIN);
  for(i = 0; i < LEARNING_BIN; i++)
  {
    wconf_buf[i] = malloc(sizeof(double) * 13);
  }

  for (wd = 0; wd < LEARNING_BIN; wd++)
  {
    for (nw = 0; nw < N_WEIGHT; nw++)
    {
      dwx = 0;
      for (md = 0; md < MACHINE_DAY; md++)
      {
        wx = 0;
        for (wd2 = 0; wd2 < LEARNING_BIN; wd2++)
        {
          for (nw2 = 0; nw2 < N_WEIGHT; nw2++)
          {
            wx += data_machine[wd2][nw2][md] * wconf[wd2][nw2];
          }
        }
        dwx += -2 * data_machine[wd][nw][md] * (data_dif[md] - (wx + INTERSEPT[0]));
      }
      wconf_buf[wd][nw] = wconf[wd][nw] - LEARNING_RATE * dwx;
    }
  }

  db = 0;
  for (md = 0; md < MACHINE_DAY; md++)
  {
    wx = 0;
    for (wd = 0; wd < LEARNING_BIN; wd++)
    {
      for (nw = 0; nw < N_WEIGHT; nw++)
      {
        wx += data_machine[wd][nw][md] * wconf[wd][nw];
      }
    }
    db += -2 * (data_dif[md] - (wx + INTERSEPT[0]));
  }

  /*===== Exchange data =====*/
  for (wd = 0; wd < LEARNING_BIN; wd++)
  {
    for (nw = 0; nw < N_WEIGHT; nw++)
    {
      wconf[wd][nw] = wconf_buf[wd][nw];
    }
  }
  INTERSEPT[0] -= LEARNING_RATE * db;

  free(wconf_buf);
}
