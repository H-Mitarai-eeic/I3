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

typedef short sample_t;

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

double f_at_PSD_max(double * PSD, long n, double T){
  double max = 0;
  double f = 0;
  long i;
  for (i = 0; i < n / 2; i++) {
    if (max < PSD[i]){
      max = PSD[i];
      f = (double)i / T;
    }
  }
  return f;
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
int check_PSD(double *PSD, int nf, double threshold){ //PSD[nf] が閾値以上なら1, 以下なら0 を返す
  if(PSD[nf] >= threshold){
    return 1;
  }
  else{
    return 0; 
  }
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
  FILE * wp = fopen("./output/fft.dat", "wb");
  if (wp == NULL) die("fopen");
  FILE * wA = fopen("./output/ACF.dat", "wb");
  if (wA == NULL) die("fopen");
  FILE * wP = fopen("./output/PSD_re.dat", "wb");
  if (wP == NULL) die("fopen");

  sample_t * buf = calloc(sizeof(sample_t), n);
  complex double * X = calloc(sizeof(complex double), n);
  complex double * Y = calloc(sizeof(complex double), n);
  complex double * PSD_comp = calloc(sizeof(complex double), n);
  complex double * PSD_comp_cp = calloc(sizeof(complex double), n);
  complex double * ACF_comp = calloc(sizeof(complex double), n);
  double *ACF_re = calloc(sizeof(double), n);
  double *PSD_re = calloc(sizeof(double), n);

  double f0, T;
  int nt;
  int nf;

  char key[KeyNum] = {0};

  T = (double)n / Fs;

  print_keyboard(key, KeyNum);
  printf("\e[%dA", 1);
  while (1) {
    /* 標準入力からn個標本を読む */
    ssize_t m = read_n(0, n * sizeof(sample_t), buf);
    if (m == 0) break;
    /* 複素数の配列に変換 */
    sample_to_complex(buf, X, n);
    /* FFT -> Y */

    fft(X, Y, n);

    CALC_PSD(Y, PSD_comp, n, T);
    copy_complex(PSD_comp, PSD_comp_cp, n);
    print_PSD(wP, PSD_comp, n, T);

    print_complex(wp, Y, n);
    fprintf(wp, "----------------\n");

    /* IFFT で　ACFを求める */
    ifft(PSD_comp_cp, ACF_comp, n);
    //f0 = f_at_PSD_max(PSD, n, (double)n / Fs);

    /* ACF, PSDの実部を取り出す */
    complex_to_re(ACF_comp, ACF_re, n);
    complex_to_re(PSD_comp, PSD_re, n);

    /* AFCのピークを求め、基本周波数を求める */
    nt = max_peak(ACF_re, n);
    f0 = Fs / nt;

    nf = rounding((double)n / nt);

    //鍵盤の何番目に値するかを求める
    clear_array(key, KeyNum);
    if(check_PSD(PSD_re, nf, threshold)){
          store2key(f0 ,key, KeyNum);
    }
    //出力
    //printf("%d, %f\n", nt, f0);
    printf("\e[%dA", 2);
    print_keyboard(key, KeyNum);
    fflush(stdout);
    /* ACFデータを書き込み */
    print_ACF(wA, ACF_re, n);
  }
  fclose(wp);
  fclose(wA);
  fclose(wP);
  return 0;
}
