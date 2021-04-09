#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#define Pi       3.14159265358979323846264338328
#define Pih      1.57079632679489661923132169164 // Pi/2
#define Piq      0.78539816339744830961566084582 // Pi/4
#define Pi2      9.86960440108935861883449099988 // Pi^2
#define Third    0.33333333333333333333333333333 // 1/3
#define TwoThird 0.66666666666666666666666666667 // 2/3
#define Sqrt_2   1.41421356237309504880168872421 // Sqrt[2]

#define TINY 1.0e-08
#define BIG  1.0e+10
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*
#define NMAX 2000

typedef struct DCOMPLEX {double r,i;} dcomplex;

void nrerror(char error_text[]);
int      *ivector( int nl,  int nh);
double   *dvector( int nl,  int nh);
long double *ldvector( int nl,  int nh);
int     **imatrix( int nrl, int nrh, int ncl, int nch);
double  **dmatrix( int nrl, int nrh, int ncl, int nch);
long double **ldmatrix( int nrl, int nrh, int ncl, int nch);
double ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
void free_ivector(int       *v, int nl,  int nh);
void free_dvector(double    *v, int nl,  int nh);
void free_ldvector(long double *v, int nl,  int nh);
void free_imatrix(int      **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double   **m, int nrl, int nrh, int ncl, int nch);
void free_ldmatrix(long double **m, int nrl, int nrh, int ncl, int nch);
void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
								int ndl, int ndh);
void fill0_dvector(double *X, int nl,  int nh);
void fill0_ivector(int *X, int nl,  int nh);
void fill0_dmatrix(double **X, int nrl, int nrh, int ncl, int nch);
void fill0_imatrix(int **X, int nrl, int nrh, int ncl, int nch);
void copy_dvector(double  *aout, double  *ain, int nl, int nh);
void copy_dmatrix(double **aout, double **ain, int nrl, int nrh, int ncl, int nch);
double norm1(double *v, int nl,  int nh);
double norm2(double *v, int nl,  int nh);
double scalarproduct(double *v, double *w, int nl,  int nh);

int minimum2(int i,int j);
int minimum3(int i,int j,int k);
int maximum2(int i,int j);
double dmaximum2(double a, double b);
int maximum3(int i,int j,int k);
int pow_int(int mantisse,int exponent);
double sinch(double x);
double sinc(double x);
double Sqrt(double x);
double sqr(double x);
long double sqrl(long double x);



dcomplex Cadd(dcomplex a, dcomplex b);
dcomplex Csub(dcomplex a, dcomplex b);
dcomplex Cmul(dcomplex a, dcomplex b);
dcomplex RCmul(double x, dcomplex a);
dcomplex Cdiv(dcomplex a, dcomplex b);
dcomplex Complex(double re, double im);
dcomplex Conjg(dcomplex z);
double Cabs(dcomplex z);

dcomplex Csqrt(dcomplex z);
dcomplex Cexp(dcomplex z);
dcomplex Clog(dcomplex z);
dcomplex Csin(dcomplex z);
dcomplex Ccos(dcomplex z);
dcomplex Ctan(dcomplex z);
dcomplex Ccot(dcomplex z);
dcomplex Csinh(dcomplex z);
dcomplex Ccosh(dcomplex z);
dcomplex Ctanh(dcomplex z);
dcomplex Ccoth(dcomplex z);

void banmul(double **a, int n, int m1, int m2, double x[], double b[]);
void bandec(double **a, int n, int m1, int m2, double **al, int indx[], double *d);
void banbks(double **a, int n, int m1, int m2, double **al, int indx[], double b[]);
void tridag(double a[], double b[], double c[], double r[], double u[],
   unsigned long n);
void ludcmp(double **a, int N, int *indx, double *d, int FLAG);
void lubksb(double **a, int N, int *indx, double *b, int FLAG);
void mprove(double **a, double **alud, int n, int indx[], double b[], double x[]);
void mprove_band(double **a, double **au, int n, int m1, int m2, double **al, int indx[], double b[], double x[]);

double plgndr(int l, int m, double x);
