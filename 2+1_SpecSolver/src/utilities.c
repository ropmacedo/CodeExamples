#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "utilities.h"


// -----------------------------------------------------------------------------------
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
// -----------------------------------------------------------------------------------
int *ivector(int nl, int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	if(nl>nh) return 0;
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
// -----------------------------------------------------------------------------------
double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	if(nl>nh) return 0;
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}
// -----------------------------------------------------------------------------------
long double *ldvector(int nl, int nh)
/* allocate a long double vector with subscript range v[nl..nh] */
{
	if(nl>nh) return 0;
	long double *v;

	v=(long double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long double)));
	if (!v) nrerror("allocation failure in ldvector()");
	return v-nl+NR_END;
}
// -----------------------------------------------------------------------------------
int **imatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	if(nrl>nrh || ncl>nch) return 0;
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
// -----------------------------------------------------------------------------------
double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	if(nrl>nrh || ncl>nch) return 0;
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
// -----------------------------------------------------------------------------------
long double**ldmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a long double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	if(nrl>nrh || ncl>nch) return 0;
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	long double**m;

	/* allocate pointers to rows */
	m=(long double**) malloc((size_t)((nrow+NR_END)*sizeof(long double*)));
	if (!m) nrerror("allocation failure 1 in ldmatrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(long double*) malloc((size_t)((nrow*ncol+NR_END)*sizeof(long double)));
	if (!m[nrl]) nrerror("allocation failure 2 in ldmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
// -----------------------------------------------------------------------------------
double ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}
// -----------------------------------------------------------------------------------
void free_ivector(int *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
// -----------------------------------------------------------------------------------
void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
// -----------------------------------------------------------------------------------
void free_ldvector(long double *v, int nl, int nh)
/* free a long double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
// -----------------------------------------------------------------------------------
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
// -----------------------------------------------------------------------------------
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
// -----------------------------------------------------------------------------------
void free_ldmatrix(long double **m, int nrl, int nrh, int ncl, int nch)
/* free a long double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
// -----------------------------------------------------------------------------------
void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch,
	int ndl, int ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
// -----------------------------------------------------------------------------------
void fill0_dvector(double *X, int nl,  int nh)
{
	int n;
	for (n=nl; n<=nh; n++) X[n]=0.;
}
// -----------------------------------------------------------------------------------
void fill0_ivector(int *X, int nl,  int nh)
{
	int n;	
	for (n=nl; n<=nh; n++) X[n]=0;
}
// -----------------------------------------------------------------------------------
void fill0_dmatrix(double **X, int nrl, int nrh, int ncl, int nch)
{
	int n, m;
	for (m=nrl; m<=nrh; m++)
		for (n=ncl; n<=nch; n++) X[m][n]=0.;
}
// -----------------------------------------------------------------------------------
void fill0_imatrix(int **X, int nrl, int nrh, int ncl, int nch)
{
	int n, m;	
	for (m=nrl; m<=nrh; m++)
		for (n=ncl; n<=nch; n++) X[m][n]=0;
}
// -----------------------------------------------------------------------------------
void copy_dvector(double *aout, double *ain, int nl,  int nh)
{
	int n;
	for (n=nl; n<=nh; n++) aout[n]=ain[n];
}
// -----------------------------------------------------------------------------------
void copy_dmatrix(double **aout, double **ain, int nrl, int nrh, int ncl, int nch)
{
	int n, m;
	for (m=nrl; m<=nrh; m++)
		for (n=ncl; n<=nch; n++) aout[m][n]=ain[m][n];
}
// -----------------------------------------------------------------------------------
double norm1(double *v, int nl,  int nh)
{
	int i;
	double result=-1;
	
	for (i=nl; i<=nh; i++)
		if (fabs(v[i]) > result) result = fabs(v[i]);
	
	return result;
}
// -----------------------------------------------------------------------------------
double norm2(double *v, int nl,  int nh)
{
	int i;
	double result=0;
	
	for (i=nl; i<=nh; i++)
		result += sqr(v[i]);
	
	return sqrt(result);
}
// -----------------------------------------------------------------------------------
double scalarproduct(double *v, double *w, int nl,  int nh)
{
	int i;
	double result=0;
	
	for (i=nl; i<=nh; i++)
		result += v[i]*w[i];
	
	return result;
}
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
int minimum2(int i,int j)
{	
	int result=i;
	if (j<result)	result = j;
	return result;
}
// -----------------------------------------------------------------------------------
int minimum3(int i,int j,int k)
{
	int result=i;
	if (j<result)	result = j;
	if (k<result)	result = k;
	return result;
}
// -----------------------------------------------------------------------------------
int maximum2(int i,int j)
{	
	int result=i;
	if (j>result)	result = j;
	return result;
}
// -----------------------------------------------------------------------------------
double dmaximum2(double a, double b)
{	
	double result=a;
	if (b>result) result = b;
	return result;
}
// -----------------------------------------------------------------------------------
int maximum3(int i,int j,int k)
{
	int result=i;
	if (j>result)	result = j;
	if (k>result)	result = k;
	return result;
}
// -----------------------------------------------------------------------------------
int pow_int(int mantisse,int exponent)
{
	int i, result =1;
	
	for (i=1; i<=exponent; i++)
		result *= mantisse;
		
	return result;
}
// -----------------------------------------------------------------------------------
double sinch(double x)
{
	double result = 1.;
	
	if (fabs(x) > TINY) result = sinh(x)/x;

	return result;
}
// -----------------------------------------------------------------------------------
double sinc(double x)
{
	double result = 1.;
	
	if (fabs(x) > TINY) result = sin(x)/x;

	return result;
}
// -----------------------------------------------------------------------------------
double Sqrt(double x)
{
	return sqrt(fabs(x));
}
// -----------------------------------------------------------------------------------
double sqr(double x)
{
	return x*x;
}
// -----------------------------------------------------------------------------------
long double sqrl(long double x)
{
	return x*x;
}
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
dcomplex Cadd(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Csub(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Cmul(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex RCmul(double x, dcomplex a)
{
	dcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Cdiv(dcomplex a, dcomplex b)
{
	dcomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Complex(double re, double im)
{
	dcomplex c;
	c.r=re;
	c.i=im;
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Conjg(dcomplex z)
{
	dcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}
// -----------------------------------------------------------------------------------
double Cabs(dcomplex z)
{
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}
// -----------------------------------------------------------------------------------
dcomplex Csqrt(dcomplex z)
{
	dcomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}
// -----------------------------------------------------------------------------------
dcomplex Cexp(dcomplex z)
{
	dcomplex c;
	double exp_r=exp(z.r);

	c.r=exp_r*cos(z.i);
	c.i=exp_r*sin(z.i);
	
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Clog(dcomplex z)
{   
	dcomplex c;
	
	c.r = 0.5*log(z.r*z.r+z.i*z.i);
	c.i = atan2(z.i,z.r);
	
	return c;
}	
// -----------------------------------------------------------------------------------
dcomplex Csin(dcomplex z)
{
	dcomplex c;
	
	c.r= sin(z.r)*cosh(z.i);
	c.i= cos(z.r)*sinh(z.i);
	
	return c;
}// -----------------------------------------------------------------------------------
dcomplex Ccos(dcomplex z)
{
	dcomplex c;
	
	c.r= cos(z.r)*cosh(z.i);
	c.i=-sin(z.r)*sinh(z.i);
	
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Ctan(dcomplex z)
{
	return Cdiv(Csin(z), Ccos(z));
}
// -----------------------------------------------------------------------------------
dcomplex Ccot(dcomplex z)
{
	return Cdiv(Ccos(z), Csin(z));
}
// -----------------------------------------------------------------------------------
dcomplex Csinh(dcomplex z)
{
	dcomplex c;
	
	c.r= sinh(z.r)*cos(z.i);
	c.i= cosh(z.r)*sin(z.i);
	
	return c;
}// -----------------------------------------------------------------------------------
dcomplex Ccosh(dcomplex z)
{
	dcomplex c;
	
	c.r= cosh(z.r)*cos(z.i);
	c.i= sinh(z.r)*sin(z.i);
	
	return c;
}
// -----------------------------------------------------------------------------------
dcomplex Ctanh(dcomplex z)
{
	return Cdiv(Csinh(z), Ccosh(z));
}
// -----------------------------------------------------------------------------------
dcomplex Ccoth(dcomplex z)
{
	return Cdiv(Ccosh(z), Csinh(z));
}
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
void banmul(double **a, int n, int m1, int m2, double x[], double b[])
{
	//Small alteration on banbks from Numerical Recipes: Here ist x and b offset 0, while
    //the matrices a, au and al are offset 1
	int i,j,k,k_o, i_o, tmploop;

	for (i=1;i<=n;i++) {
		i_o = i-1;
		k=i-m1-1;
		k_o = k-1;
		tmploop = minimum2(m1+m2+1,n-k);

		b[i_o]=0.0;
		for (j = maximum2(1,1-k); j<=tmploop;j++) b[i_o] += a[i][j]*x[j+k_o];
	}
}
//-----------------------------------------------------------------------------
void bandec(double **a, int n, int m1, int m2, double **al, int indx[], double *d)
{
	unsigned long i,j,k,l;
	int mm;
	double temp;

	mm=m1+m2+1;
	l=m1;
	for (i=1;i<=m1;i++) {
		for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];
		l--;
		for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
	}
	*d=1.0;
	l=m1;
	for (k=1;k<=n;k++) {
		temp=a[k][1];
		i=k;
		if (l < n) l++;
		for (j=k+1;j<=l;j++) {
			if (fabs(a[j][1]) > fabs(temp)) {
				temp=a[j][1];
				i=j;
			}
		}
		indx[k]=i;
		if (temp == 0.0) a[k][1]=TINY;
		if (i != k) {
			*d = -(*d);
			for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
		}
		for (i=k+1;i<=l;i++) {
			temp=a[i][1]/a[k][1];
			al[k][i-k]=temp;
			for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-temp*a[k][j];
			a[i][mm]=0.0;
		}
	}
}
//-----------------------------------------------------------------------------
void banbks(double **a, int n, int m1, int m2, double **al, int indx[], double b[])
{
  //Small alteration on banbks from Numerical Recipes: Here ist b offset 0, while
  //the matrices a and al are offset 1
	unsigned long i,k,l;
	int mm;
	double temp;

	mm=m1+m2+1;
	l=m1;
	for (k=1;k<=n;k++) {
		i=indx[k];
		if (i != k) SWAP(b[k-1],b[i-1])
		if (l < n) l++;
		for (i=k+1;i<=l;i++) b[i-1] -= al[k][i-k]*b[k-1];
	}
	l=1;
	for (i=n;i>=1;i--) {
		temp=b[i-1];
		for (k=2;k<=l;k++) temp -= a[i][k]*b[k+i-1-1];
		b[i-1]=temp/a[i][1];
		if (l < mm) l++;
	}
}
// -----------------------------------------------------------------------------------
void tridag(double a[], double b[], double c[], double r[], double u[],
	unsigned long n)
{
	unsigned long j;
	double bet,*gam;

	gam=dvector(1,n);
	if (b[1] == 0.0) nrerror("Error 1 in tridag");
	u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0)	nrerror("Error 2 in tridag");
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
	free_dvector(gam,1,n);
}
// -----------------------------------------------------------------------------------
void ludcmp(double **a, int N, int *indx, double *d,int FLAG)
{
  /*Small change in the NR routine: Introduce FLAG to
  caracterise offset of input quantites*/


	int i,imax=0,j,k;
	double big,dum,sum,temp;
	double *vv;
	
switch (FLAG)
{
  
  case 0:
  {
    // Version of 'ludcmp' of the numerical recipes for
    // matrices a[0:n-1][0:n-1]


int n=N+1;
	vv=dvector(0,n-1);
	*d=1.0;

	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0){
			free_dvector(vv,0,n-1); 
			printf("ludcmp: Row i=%d is identical 0\n",i);
			nrerror("Singular matrix in routine ludcmp");
		}
		vv[i]=1.0/big;
	}

	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,0,n-1);
	

  break;
  }
  
  case 1:
  {	// Version of 'ludcmp' of the numerical recipes for
	// matrices a[1:n][1:n] and vectors b[1:n]
int n=N;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0){
			free_dvector(vv,1,n);
			printf("ludcmp_1: Row i=%d is identical 0\n",i);
			nrerror("Singular matrix in routine ludcmp_1");
		}
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
    break;
  }
  
  default:
  {
    printf("\nError in routine ludcmp: vector and matrices are not offset 0 or 1...\n");
    exit(1);
  }
 
}
 
}
// -----------------------------------------------------------------------------------

/* (C) Copr. 1986-92 Numerical Recipes Software V,3. */

void lubksb(double **a, int N, int *indx, double *b, int FLAG)
{
  /*Small change in the NR routine: Introduce FLAG to
  caracterise offset of input quantites*/


	int i,ii=0,ip,j;
	double sum;

	
	switch (FLAG)
	{
	  
	  case 0:
	  {
	   
	    // Version of 'lubksp' of the numerical recipes for
	    // matrices a[0:n-1][0:n-1]
	int n=N+1;

		for (i=0;i<n;i++) {
			ip=indx[i];
			sum=b[ip];
			b[ip]=b[i];
			if (ii)
				for (j=ii-1;j<=i-1;j++) sum -= a[i][j]*b[j];
			else if (sum) ii=i+1;
	/*Modification when going from offset 1 to offset 0: 
	Condition if(ii) checks if ii is diferent from 0. In the first iteration, this is false, so it used set ii=i.
	But here ii would continue to be zero, and it has to change. So I added ii=i+1. However, ii is also used in the
	loop in j, which in turns has to start from 0 at the first time it is called. Therefor adding '-1' in the j range.*/
			b[i]=sum;
		}
		for (i=n-1;i>=0;i--) {
			sum=b[i];
			for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
			b[i]=sum/a[i][i];
		}

	  break;
	  }
	  
	  case 1:
	  {	// Version of 'lubksb' of the numerical recipes for
		// matrices a[1:n][1:n] and vectors b[1:n]
		int n=N;
		
			for (i=1;i<=n;i++) {
			ip=indx[i];
			sum=b[ip];
			b[ip]=b[i];
			if (ii)
				for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
			else if (sum) ii=i;
			b[i]=sum;
		}
		for (i=n;i>=1;i--) {
			sum=b[i];
			for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
			b[i]=sum/a[i][i];
		}

	    break;
	  }
	  
	  default:
	  {
	    printf("\nError in routine lubksp: vector and matrices are not offset 0 or 1...\n");
	    exit(1);
	  }
	 
	}
}
//-------------------------------------------------------------------
void mprove(double **a, double **alud, int n, int indx[], double b[], double x[])
{
	int j,i;
	double sdp;
	double *r;

	r=dvector(1,n);
	for (i=1;i<=n;i++) {
		sdp = -b[i];
		for (j=1;j<=n;j++) sdp += a[i][j]*x[j];
		r[i]=sdp;
	}
	lubksb(alud,n,indx,r, 1);
	for (i=1;i<=n;i++) x[i] -= r[i];
	free_dvector(r,1,n);
}
//-------------------------------------------------------------------
void mprove_band(double **a, double **au, int n, int m1, int m2, double **al, int indx[], double b[], double x[])
// (double **a, double **alud, int n, int indx[], double b[], double x[])
{
	int i;
	double *r;

	r=dvector(1,n);

	banmul(a, n, m1, m2, x, r);
	for (i=1;i<=n;i++)
		r[i] -= b[i];

	banbks(au, n, m1, m2, al, indx, r);
	
	for (i=1;i<=n;i++) x[i] -= r[i];
	free_dvector(r,1,n);
}
//------------------------------------------------------------------
double plgndr(int l, int m, double x)
{
	void nrerror(char error_text[]);
	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;

	if (m < 0 || m > l || fabs(x) > 1.0)
		nrerror("Bad arguments in routine plgndr");
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=m+2;ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}
