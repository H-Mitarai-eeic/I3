#include <assert.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define Fs 44100.0
#define KeyNum 88
#define A0_FREQ 27.5;
#define Default_Threshold 100000
#define PeakNum 5 //ACFから一度に取り出すピーク数の上限
#define ACF_Cal_Times 8  //ACFの計算回数

typedef short sample_t;
void sortarrays(double *peaks, int *index, int N);
void die(char * s) {
  perror(s); 
  exit(1);
}
/* fd から 必ず n バイト読み, bufへ書く.
   n バイト未満でEOFに達したら, 残りは0で埋める.
   fd から読み出されたバイト数を返す */
ssize_t read_n(int fd, ssize_t n, void * buf) {
  ssize_t re = 0;
  while (re < n) {
    ssize_t r = read(fd, buf + re, n - re);
    if (r == -1) die("read");
    if (r == 0) break;
    re += r;
  }
  memset(buf + re, 0, n - re);
  return re;
}
/* 標本(整数)を複素数へ変換 */
void sample_to_complex(sample_t * s, 
		       complex double * X, 
		       long n) {
  long i;
  for (i = 0; i < n; i++) X[i] = s[i];
}

/* 複素数を標本(整数)へ変換. 虚数部分は無視 */
void complex_to_sample(complex double * X, 
		       sample_t * s, 
		       long n) {
  long i;
  for (i = 0; i < n; i++) {
    s[i] = creal(X[i]);
  }
}

/* 高速(逆)フーリエ変換;
   w は1のn乗根.
   フーリエ変換の場合   偏角 -2 pi / n
   逆フーリエ変換の場合 偏角  2 pi / n
   xが入力でyが出力.
   xも破壊される
 */
void fft_r(complex double * x, 
	   complex double * y, 
	   long n, 
	   complex double w) {
  if (n == 1) { y[0] = x[0]; }
  else {
    complex double W = 1.0; 
    long i;
    for (i = 0; i < n/2; i++) {
      y[i]     =     (x[i] + x[i+n/2]); /* 偶数行 */
      y[i+n/2] = W * (x[i] - x[i+n/2]); /* 奇数行 */
      W *= w;
    }
    fft_r(y,     x,     n/2, w * w);
    fft_r(y+n/2, x+n/2, n/2, w * w);
    for (i = 0; i < n/2; i++) {
      y[2*i]   = x[i];
      y[2*i+1] = x[i+n/2];
    }
  }
}

void fft(complex double * x, 
	 complex double * y, 
	 long n) {
  long i;
  double arg = 2.0 * M_PI / n;
  complex double w = cos(arg) - 1.0j * sin(arg);
  fft_r(x, y, n, w);
  for (i = 0; i < n; i++) y[i] /= n;
}

void ifft(complex double * y, 
	  complex double * x, 
	  long n) {
  double arg = 2.0 * M_PI / n;
  complex double w = cos(arg) + 1.0j * sin(arg);
  fft_r(y, x, n, w);
}

int pow2check(long N) {
  long n = N;
  while (n > 1) {
    if (n % 2) return 0;
    n = n / 2;
  }
  return 1;
}

void print_complex(FILE * wp, 
		   complex double * Y, long n) {
  long i;
  for (i = 0; i < n; i++) {
    fprintf(wp, "%ld %f %f %f %f\n", 
	    i, 
	    creal(Y[i]), cimag(Y[i]),
	    cabs(Y[i]), atan2(cimag(Y[i]), creal(Y[i])));
  }
}
void print_PSD(FILE * wp, complex double * PSD, long n, double T) {
  long i;
  for (i = 0; i < n; i++) {
    //fprintf(wp, "%f %f\n", (double)i / T, PSD[i]);
    fprintf(wp, "%ld %f %f\n", i, creal(PSD[i]), cimag(PSD[i]));
  }
}

void print_ACF(FILE * wA, double *ACF_re, long n) {
  long i;
  for (i = 0; i < n; i++) {
    fprintf(wA, "%ld %f\n", i, ACF_re[i]);
  }
}
void CALC_PSD(complex double * x, complex double * PSD, long n, double T){
  long i;
  for(i = 0; i < n; i++){
    PSD[i] = (x[i] * conj(x[i])) / T;
  }
}
double max_peak(double *signal, long n){
  double max_peak = 0;
  int i;
  int index;
  for(i = 1; i < n / 2 - 1; i++ ){
    if(signal[i - 1] <= signal[i] && signal[i] > signal[i + 1]){
      if (max_peak < signal[i]){
        max_peak = signal[i];
        index = i;
      }
    }
  }
  return index;
}
int multi_peak(double *signal, int nt, int *peak_index, int NUM){  //0 ~ ntまでの領域で NUM個大きい順にpeak インデックスを取得し, peak_indexに書き込む(書き込まれる順番は大きい順ではない) 返り値は書き込んだピーク数
  double *signal_MAX = calloc(sizeof(double), NUM);
  //double signal_low_lim = 0;
  int i, j;
  j = 0;
  for(i = 1; i <= nt; i++){
      if(signal[i - 1] <= signal[i] && signal[i] > signal[i + 1]){
        if(j == 0){
          peak_index[0] = i;
          signal_MAX[0] = signal[i];
          j++;
        }
        else if(j < NUM){
          peak_index[j] = i;
          signal_MAX[j] = signal[i];
          j++;
          if(j == NUM){
            //sort
            sortarrays(signal_MAX, peak_index, NUM);
          }
        }
        else if(j == NUM){
          if (signal_MAX[NUM - 1] < signal[i]){
            peak_index[NUM - 1] = i;
            signal_MAX[NUM - 1] = signal[i];
            sortarrays(signal_MAX, peak_index, NUM);
          }
        }
      }
  }
  sortarrays(signal_MAX, peak_index, j);
  free(signal_MAX);
  return j;
}
void complex_to_re(complex double * X,double * s, long n) { //実部のみ取り出す
  long i;
  for (i = 0; i < n; i++) {
    s[i] = creal(X[i]);
  }
}
void clear_array(char *array, int n){ //配列要素を0にする
  int i;
  for (i = 0; i < n; i++){
    array[i] = 0;
  }
}
void store2key(double f, char *key, int num_key){
  int i;
  double A0 = A0_FREQ;
  double Center_freq = A0;

  for(i = 0; i < num_key; i++){
    if(Center_freq * pow(2, -1.0/24) <= f && f < Center_freq * pow(2, 1.0/24)){
      key[i] = 1;
    }
    Center_freq = Center_freq * pow(2, 1.0/12.0);
  }
}
int check_PSD(double *PSD, int nf, double threshold){ //PSD[nf - 1] + PSD[nf] + PSD[nf + 1] が閾値以上なら1, 以下なら0 を返す
  /*if(PSD[nf - 1] + PSD[nf] + PSD[nf + 1] >= threshold){
    return 1;
  }
  else{
    return 0;
  }*/
  
  if(PSD[nf] >= threshold){
    return 1;
  }
  else{
    return 0;
  }
}
void update_PSD(complex double *PSD, int nf, long n){
  int i, j;
  complex double PSD_peak;
  PSD_peak = PSD[nf];
  //PSD_peak = PSD[nf - 1] + PSD[nf] + PSD[nf + 1];
  /*
  for(i = nf; i < n/2 - 1; i += nf){
    for (j = -1; j <= 1; j++){
      if(creal(PSD[i + j]) <= creal(PSD_peak)){
        PSD[i + j] = 0;
        PSD[n - 1 - i - j] = 0;
      }
      else{
        PSD[i + j] -= PSD_peak;
        PSD[n - 1 - i - j] -= PSD_peak;
      }
    }
  }*/
  for(i = nf; i < n/2 - 1; i += nf){
      if(creal(PSD[i]) <= creal(PSD_peak)){
        PSD[i] = 0;
        PSD[n - 1 - i] = 0;
      }
      else{
        PSD[i] -= PSD_peak;
        PSD[n - 1 - i] -= PSD_peak;
      }
    }
  return;
}
int rounding(double A){  //四捨五入
  double dec;
  dec = A - (int)A;
  if(dec >= 0.5){
    return (int)A + 1;
  }
  else
  return (int)A;
}
void sortarrays(double *peaks, int *index, int N){  //peakの値のよって２つの配列を降順にソート
  int i, j;
  double temp_p;
  int temp_i;
  for(i = 1; i <= N - 1; i++){
    for(j = N - 1; j >= i; j--){
        if(peaks[j] > peaks[j - 1]){
            temp_p = peaks[j - 1];
            temp_i = index[j - 1];
            peaks[j - 1] = peaks[j];
            index[j - 1] = index[j];
            peaks[j] = temp_p;
            index[j] = temp_i;
        }
    }
  }
}
void print_keyboard(char *key, int N){
    //int N = 52;
    int i;
    printf("\x1b[30m\x1b[47m|");
    for(i = 0; i < N; i++){
        if(i % 12 == 1 || i % 12 == 4 || i % 12 == 6 || i % 12 == 9 || i % 12 == 11){
            if(key[i] == 1){
                printf("\x1b[41m \x1b[47m");
            }
            else{
                printf("\x1b[40m \x1b[47m");
            }
        }
        else if (i % 12 == 2 || i % 12 == 7){
            printf(" |");
        }
        else{
            printf(" ");
        }
    }
    printf("\n\x1b[47m|");
    for (i = 0; i < N; i++){
        if(i % 12 == 1 || i % 12 == 4 || i % 12 == 6 || i % 12 == 9 || i % 12 == 11){
            printf("|");
        }
        else if (i % 12 == 2 || i % 12 == 7){
            if (key[i] == 1){
                printf("\x1b[41m_\x1b[47m|");            
            }
            else{
                printf("_|");
            }
        }
        else{
            if (key[i] == 1){
                printf("\x1b[41m_\x1b[47m");            
            }
            else{
                printf("_");
            }
        }
    }
    printf("|\x1b[49m\x1b[39m\n");
    fflush(stdout);
}
void copy_complex(complex double *Original, complex double *Copy, long n){
  long i;
  for (i = 0; i < n; i++){
    Copy[i] = Original[i];
  }
  return;
}
int main(int argc, char ** argv) {
  (void)argc;
  long n = atol(argv[1]);
  double threshold;
  if (argc < 3){
    threshold = Default_Threshold;
  }
  else{
    threshold = atof(argv[2]);
  }

  if (!pow2check(n)) {
    fprintf(stderr, "error : n (%ld) not a power of two\n", n);
    exit(1);
  }
  /*
  FILE * wp = fopen("./output/fft.dat", "wb");
  if (wp == NULL) die("fopen");
  FILE * wA = fopen("./output/ACF.dat", "wb");
  if (wA == NULL) die("fopen");
  FILE * wP = fopen("./output/PSD_re.dat", "wb");
  if (wP == NULL) die("fopen");
  */
  sample_t * buf = calloc(sizeof(sample_t), n);
  complex double * X = calloc(sizeof(complex double), n);
  complex double * Y = calloc(sizeof(complex double), n);
  complex double * PSD_comp = calloc(sizeof(complex double), n);
  complex double * PSD_comp_cp = calloc(sizeof(complex double), n);
  complex double * ACF_comp = calloc(sizeof(complex double), n);
  double *ACF_re = calloc(sizeof(double), n);
  double *PSD_re = calloc(sizeof(double), n);

  double f0, T;
  int nt, nf, Np; //Np ピーク数
  int nt_list[PeakNum] = {0};
  int i, j;  //counter

  char key[KeyNum] = {0};

  T = (double)n / Fs;
  print_keyboard(key, KeyNum);
  printf("\e[%dA", 1);
  while (1) {
    /* 標準入力からn個標本を読む */
    ssize_t m = read_n(0, n * sizeof(sample_t), buf);
    if (m == 0) {
      printf("faild to read data\n");
      fflush(stdout);
      break;
    }
    /* 複素数の配列に変換 */
    sample_to_complex(buf, X, n);
    /* FFT -> Y */

    fft(X, Y, n);

    CALC_PSD(Y, PSD_comp, n, T);
    
    //print_PSD(wP, PSD_comp, n, T);

    //print_complex(wp, Y, n);
    //fprintf(wp, "----------------\n");

    for(j = 0; j < ACF_Cal_Times; j++){
      //print_PSD(wP, PSD_comp, n, T);
      /* IFFT で　ACFを求める */
      copy_complex(PSD_comp, PSD_comp_cp, n);
      ifft(PSD_comp_cp, ACF_comp, n);
      //f0 = f_at_PSD_max(PSD, n, (double)n / Fs);

      /* ACF, PSDの実部を取り出す */
      complex_to_re(ACF_comp, ACF_re, n);
      complex_to_re(PSD_comp, PSD_re, n);

      /* AFCのピークを求め、基本周波数を求める */
      nt = max_peak(ACF_re, n);
      Np = multi_peak(ACF_re, nt, nt_list, PeakNum);
      //printf("### NP = %d #######\n", Np);
      fflush(stdout);
      for (i = 0; i < Np; i++){
        f0 = Fs / nt_list[i];
        ///printf("f0 = %f", f0);
        //fflush(stdout);
        nf = rounding((double)n / nt_list[i]);
        //printf("nt, = %d, nf = %d, f0 = %f\n",nt_list[i], nf, f0);
        //fflush(stdout);


        //鍵盤の何番目に値するかを求める
        if(check_PSD(PSD_re, nf, threshold)){
              store2key(f0 ,key, KeyNum);
              update_PSD(PSD_comp, nf, n);
              break;
              //printf("stored\n");
              //fflush(stdout);
        }
      }
    }
    //出力
    //printf("%d, %f\n", n0, f0);
    printf("\e[%dA", 2);
    print_keyboard(key, KeyNum);
    fflush(stdout);
    clear_array(key, KeyNum);
    /* ACFデータを書き込み */
    ///print_ACF(wA, ACF_re, n);
  }
  //fclose(wp);
  //fclose(wA);
  //fclose(wP);
  return 0;
}
