#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "utilities.h"
#include "spectral_utilities.h"

//---------------------------------------------------------------
// ********* Clenshaw Algorithms ********************************
//---------------------------------------------------------------
double Clenshaw_cos(double *alpha, int N, double phi)
{
	double bk = 0, bkp1 = 0, bkp2, cos_phi = cos(phi), c2 = 2.*cos_phi;
	int k;
	
	for (k=N; k>=0; k--){
		bkp2 = bkp1;
		bkp1 = bk;
		bk   = alpha[k] + c2*bkp1 - bkp2;
	}
	return bk - bkp1*cos_phi - 0.5*alpha[0];
}
//---------------------------------------------------------------
double Clenshaw_sin(double *beta, int N, double phi)
{
	double bk = 0, bkp1 = 0, bkp2, cos_phi = cos(phi), sin_phi = sin(phi), c2 = 2.*cos_phi;
	int k;
	
	for (k=N; k>=1; k--){
		bkp2 = bkp1;
		bkp1 = bk;
		bk   = beta[k] + c2*bkp1 - bkp2;
	}
	return bk*sin_phi;
}
//---------------------------------------------------------------
double Clenshaw_Fourier(double *alpha, double *beta, int N, double phi)
{
	double bck = 0, bckp1 = 0, bckp2, bsk = 0, bskp1 = 0, bskp2, 
		cos_phi = cos(phi), sin_phi = sin(phi), c2 = 2.*cos_phi;
	int k;
	
	for (k=N; k>=0; k--){
		bckp2 = bckp1;
		bckp1 = bck;
		bck   = alpha[k] + c2*bckp1 - bckp2;
		if(k>0){
			bskp2 = bskp1;
			bskp1 = bsk;
			bsk   = beta[k]  + c2*bskp1 - bskp2;
		}
	}
	return bck - bckp1*cos_phi - 0.5*alpha[0] + bsk*sin_phi;
}	
//---------------------------------------------------------------
double Clenshaw_Chebyshev(double *c, int N, double x)
{
	double bk = 0, bkp1 = 0, bkp2 = 0, x2 = 2*x;
	int k;
	
	for (k=N; k>=0; k--){
		bkp2 = bkp1;
		bkp1 = bk;
		bk   = c[k] + x2*bkp1 - bkp2;
	}
	return 0.5*(bk - bkp2);
}
//---------------------------------------------------------------
// ********* Fourier-Coefficients  ****************************
//---------------------------------------------------------------
void Fourier_Coefficients(double *psi, double *alpha, double *beta, int N)
{
	int m;
	double N2=2*N, fac1 = 2./(N2+1), delta_phi = Pi*fac1, aux = 0.5*psi[0]*fac1;
	
	for(m=0; m<=N; m++){ // phi_m = 2*Pi*m/(2N+1)
		alpha[m] = Clenshaw_cos(psi, N2, delta_phi*m)*fac1 + aux;
		beta[m]  = Clenshaw_sin(psi, N2, delta_phi*m)*fac1;
	}
}
//---------------------------------------------------------------
// ********* Chebyshev-Coefficients  ****************************
//---------------------------------------------------------------
void Chebyshev_Coefficients_Radau_RHS(double *psi, double *c, int N)
{
	int m;
	double fac1 = 1./(2*N+1), delta_phi = 2.*Pi*fac1, fac3 = 4.*fac1;
	
	for(m=0; m<=N; m++)
		c[m] = Clenshaw_cos(psi, N, delta_phi*m)*fac3; // phi_m = 2*Pi*m/(2N+1)
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Radau_LHS(double *psi, double *c, int N)
{
	int m, sign = 1;
	double fac1 = 1./(2*N+1), delta_phi = 2.*Pi*fac1, fac3 = 4.*fac1;
	
	for(m=0; m<=N; m++){
		c[m] = Clenshaw_cos(psi, N, delta_phi*m)*fac3*sign; // tilde_phi_m = 2*Pi*m/(2N+1)
		sign = -sign;
	}		
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Gauss(double *psi, double *c, int N)
{
	int m;
	double fac1 = 1./(N+1), delta_phi = Pi*fac1, fac3 = 2.*fac1, aux = 0.5*psi[0];
	
	for(m=0; m<=N; m++){
		double 
			tilde_phi_m = delta_phi*m,  // tilde_phi_m = Pi*m/(N+1)
			arg         = 0.5*tilde_phi_m,
			fac_cos     = fac3*cos(arg),
			fac_sin     = fac3*sin(arg);
			
		c[m] =    fac_cos*(Clenshaw_cos(psi, N, tilde_phi_m) + aux)
		       -  fac_sin* Clenshaw_sin(psi, N, tilde_phi_m);
	}
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Lobatto(double *psi, double *c, int N)
{
	int m, sign = 1;
	double fac1 = 1./N, delta_phi = Pi*fac1, fac3 = 2.*fac1, aux = 0.5*psi[N];
	
	for(m=0; m<=N; m++){
		c[m] = fac3*(Clenshaw_cos(psi, N, delta_phi*m) - aux*sign); // phi_m = Pi*m/N
		sign = -sign;
	}
	c[N] *= 0.5;
}
//---------------------------------------------------------------------------------------
void Chebyshev_Coefficients(double *psi, double *c, int N, char *grid)
{
    if(strcmp( grid,"Radau_RHS") ==0)
      Chebyshev_Coefficients_Radau_RHS(psi, c, N);
    else if(strcmp( grid,"Radau_LHS")==0)
      Chebyshev_Coefficients_Radau_LHS(psi, c, N);
    else if(strcmp( grid,"Gauss")==0)
      Chebyshev_Coefficients_Gauss(psi, c, N);
    else if(strcmp( grid,"Lobatto")==0)
      Chebyshev_Coefficients_Lobatto(psi, c, N);
    else{
      printf("Error in Chebyshev_Coefficients: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto\n grid was: %s\n", grid);
      exit(1);
    }
}
//---------------------------------------------------------------
// ********* Fourier-Collocation values from Coefficients ****
//---------------------------------------------------------------
void Fourier_Collocations(double *psi, double *alpha, double *beta, int N)
{
	int j, N2=2*N;
	double h = 2.*Pi/(N2+1); // h: Stepsize in phi
	
	for(j=0; j<=N2; j++)
		psi[j] = Clenshaw_Fourier(alpha, beta, N, h*j);
}
//---------------------------------------------------------------
// ********* Chebyshev- Collocation values from Coefficients ****
//---------------------------------------------------------------
void Chebyshev_Collocations_Radau_RHS(double *psi, double *c, int N)
{
	int j;
	double h = 2.*Pi/(2*N+1); // h: Stepsize in phi=arccos(x)
	for(j=0; j<=N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, cos(h*j));
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Radau_LHS(double *psi, double *c, int N)
{
	int j;
	double h = 2.*Pi/(2*N+1); // h: Stepsize in phi=arccos(x)
	for(j=0; j<=N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, -cos(h*j));
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Gauss(double *psi, double *c, int N)
{
	int j;
	double h = Pi/(N+1), hh = 0.5*h; // h: Stepsize in phi=arccos(x)
	for(j=0; j<=N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, cos(h*j+hh));
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Lobatto(double *psi, double *c, int N)
{
	int j;
	double h = Pi/N; // h: Stepsize in phi=arccos(x)
	for(j=0; j<=N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, cos(h*j));
}
//---------------------------------------------------------------------------------------
void Chebyshev_Collocations(double *psi, double *c, int N, char *grid)
{
    if(strcmp( grid,"Radau_RHS") ==0)
      Chebyshev_Collocations_Radau_RHS(psi, c, N);
    else if(strcmp( grid,"Radau_LHS")==0)
      Chebyshev_Collocations_Radau_LHS(psi, c, N);
    else if(strcmp( grid,"Gauss")==0)
      Chebyshev_Collocations_Gauss(psi, c, N);
    else if(strcmp( grid,"Lobatto")==0)
      Chebyshev_Collocations_Lobatto(psi, c, N);
    else{
      printf("Error in Chebyshev_Collocations: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto\n grid was: %s\n", grid);
      exit(1);
    }
}
//---------------------------------------------------------------
// ********* Fourier-Coefficients of the derivative ***********
//---------------------------------------------------------------
void Fourier_Coefficients_Derivative(double *alpha, double *beta, double *d_alpha, double *d_beta, int N)
{
	int k;
	
	d_alpha[0] = 0.; 
	for(k=1; k<=N; k++){
		d_alpha[k] =  beta[k]*k;
		d_beta[k]  = -alpha[k]*k;
	}
}
//---------------------------------------------------------------
// ********* Chebyshev-Coefficients of the derivative ***********
//---------------------------------------------------------------
void Chebyshev_Coefficients_Derivative(double *c, double *dc, int N)
{
	int k;
	
	dc[N]   = 0.;
	dc[N-1] = 2.*N*c[N];
	for(k=N-1; k>=1; k--)
		dc[k-1] = 2.*k*c[k] + dc[k+1];
}
//---------------------------------------------------------------
// ********* Chebyshev-Coefficients of the integral *************
//---------------------------------------------------------------
void Chebyshev_Coefficients_Integral(double *c, double *C, int N)
{ // undetermined value C[0] (put to 0); C is indexed 0..N+1
	int k;
	
	for(k=1;k<=N;k++)
		C[k] = 0.5*(c[k-1]-c[k+1])/k;
	
	C[N+1] = 0.5*c[N]/(N+1);
	C[0]   = 0.;
}
//---------------------------------------------------------------
// ********* Chebyshev- definite integrals **********************
//---------------------------------------------------------------
double Chebyshev_Definite_Integral_Coefficients(double *c, int N)
{ 
	int k,  kmax = 2*floor(0.5*N);
	double sum = 0.;
	
	for(k=kmax; k>=2; k -= 2)
		sum += c[k]/(-1.+k*k);
	
	return c[0] - 2.*sum;
}
//---------------------------------------------------------------
void Chebyshev_Integration_Vector_Gauss(int N, double *I)
{
	int j;
	double *c, h = Pi/(N+1), hh = 0.5*h, aux = 4./(N+1);
	
	c = dvector(0, N);
	
	for(j=0; j<=N; j++){
		if(j % 2 == 0) c[j] = aux/(1.-j*j);
		else           c[j] = 0.;
	}
	
	for(j=0; j<=N; j++)
		I[j] = Clenshaw_Chebyshev(c, N, cos(h*j+hh));
	
	free_dvector(c, 0, N);
}
//---------------------------------------------------------------
double Chebyshev_Definite_Integral_Collocations(double *psi, double *I, int N)
{
	int j;
	double s=0.;
	
	for(j=0; j<=N; j++)
		s += I[j]*psi[j];
	
	return s;
}
//---------------------------------------------------------------
double Chebyshev_Integration_VectorComponent_Lobatto(int i, int N)
{
	double G, th_i = Pi*i/N;;
	int k, NN=floor(N/2);
	
	G=1.;
	for(k=1; k<=NN; k++){
		int delta = 2*k==N? 1:0;
		double aux = i*(i-N)==0? 1.:cos(2*k*th_i);
	
		G-=(2.-delta)/(4.*k*k-1.)*aux;
	}
	
	G*=i*(i-N)==0? 1./N:2./N;
	
	return G;
	
}
//---------------------------------------------------------------
double ChebyshevLobatto_Interpolation_MatrixComponents(int alpha, int N_alpha, int i, int N){
// Given a Lobatto-grid {x_i} with i = 0 ... N, x_i \in [-1,1]
// Output interpolation at the point {x_alpha} alpha = 0... N_alpha  x_alpha \in[-1,1], so that 
// f_alpha =  I_alpha^i f_i 


double I_alpha_i, th_alpha, th_i;

th_alpha = Pi*alpha/N_alpha;
th_i = Pi*i/N;



int j;

I_alpha_i = 1.;

for(j=1;j<=N; j++){
		int delta = j==N? 1:0;
		double aux = i==N? pow(-1.,j):cos(j*th_i);
	
		I_alpha_i+= (2.-delta)*cos(j*th_alpha)*aux;
	}
	
I_alpha_i*= i*(i-N)==0? 0.5/N:1./N;


return I_alpha_i;
}
// ---------------------------------------------------------------
// void ChebyshevLobatto_Interpolation_Matrix(int N_alpha, int N, double **I){
// 
// int alpha, i;
// for(alpha=0; alpha<=N_alpha; alpha++)
// 	for(i=0; i<=N; i++)
// 		I[alpha][I]=ChebyshevLobatto_Interpolation_MatrixComponents(alpha, N_alpha, i, N);
// 
// 
// return;
// }

//---------------------------------------------------------------
double Chebyshev_Definite_Integral_Collocations_Matrix(int i, int j, int N){

int alpha, N_alpha = 2*N;

double H_i_j;


H_i_j = 0.;
for(alpha = 0; alpha<=N_alpha; alpha++){
	double I_alpha_i =  ChebyshevLobatto_Interpolation_MatrixComponents(alpha, N_alpha, i, N),
	       I_alpha_j =  ChebyshevLobatto_Interpolation_MatrixComponents(alpha, N_alpha, j, N),
	       G_alpha = Chebyshev_Integration_VectorComponent_Lobatto(alpha, N_alpha);;
			   
	H_i_j += I_alpha_i*G_alpha*I_alpha_j;
}


return H_i_j;
}
//---------------------------------------------------------------
double Chebyshev_Definite_Integral_Collocations_TwoFunc(double *f, double *g, double **H, int N){

	int i,j;
	double I=0.;
	for(i=0; i<=N; i++)
		for(j=0;j<=N; j++){
		I+=f[i]*H[i][j]*g[j];
		}

return I;
}


//---------------------------------------------------------------
// ********* Chebyshev- Differentiation matrices ****************
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrix_Radau_RHS(int N, double **D1)
{
	long double h = 2.*Pi/(2*N+1), m1m = -1.; // h: Stepsize in phi=arccosl(x)
	int m, j;
	
	for(m=0; m<=N; m++){
		long double argm = h*m, argm2 = 0.5*argm,
			root_1pxm = Sqrt_2*cosl(argm2),                          // root_1pxm = Sqrt[1+xm]
			root_1mxm = Sqrt_2*sinl(argm2), xmm1 = -sqrl(root_1mxm), // root_1mxm = Sqrt[1-xm]
			m1j = -1.;  
		m1m *= -1;  // m1m = (-1)^m
		D1[m][m]=0.;
		for(j=0; j<=N; j++){
			long double argj = h*j, argj2 = 0.5*argj,
				root_1pxj = Sqrt_2*cosl(argj2),                          // root_1pxj = Sqrt[1+xj],
				root_1mxj = Sqrt_2*sinl(argj2), xjm1 = -sqrl(root_1mxj), // root_1mxj = Sqrt[1-xj], xjm1 = xj-1
				// xmmxj = xm-xj = cos(argm)-cos(argj)= -2.*sin(argm2+argj2)*sin(argm2-argj2)
				xmmxj     = -2.*sinl(argm2+argj2)*sinl(argm2-argj2);     
			m1j *= -1;   // m1j = (-1)^j
			if(m == 0 && j != 0)           D1[m][j] = -m1j*Sqrt_2*root_1pxj/xjm1;
			if(m != 0 && j == 0)           D1[m][j] = m1m/(Sqrt_2*xmm1*root_1pxm);
			if(m != 0 && j != 0 && m != j) D1[m][j] = m1j*m1m*root_1pxj/(xmmxj*root_1pxm);
			if(m != j) D1[m][m] -= D1[m][j];
		}
	}
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrix_Radau_LHS(int N, double **D1)
{
	int m, j;
	Chebyshev_Differentiation_Matrix_Radau_RHS(N, D1);
	for(m=0; m<=N; m++)
		for(j=0; j<=N; j++)
			D1[m][j] *= -1.;
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrix_Gauss(int N, double **D1)
{
	long double h = Pi/(N+1), hh = 0.5*h, m1m = -1.; // h: Stepsize in phi=arccosl(x)
	int m, j;
	
	for(m=0; m<=N; m++){
		long double argm = h*m+hh, argm2 = 0.5*argm, // arg2m = argj/m
			root_1pxm = Sqrt_2*cosl(argm2),            // root_1pxm = Sqrt[1+xm]
			root_1mxm = Sqrt_2*sinl(argm2),            // root_1mxm = Sqrt[1-xm]
			m1j = -1.;  
		m1m *= -1;  // m1m = (-1)^m
		D1[m][m]=0.;
		for(j=0; j<=N; j++){
			long double argj = h*j+hh, argj2 = 0.5*argj,           // arg2j = argj/2
				root_1pxj = Sqrt_2*cosl(argj2),                      // root_1pxj = Sqrt[1+xj]
				root_1mxj = Sqrt_2*sinl(argj2),                      // root_1mxj = Sqrt[1-xj]
				// xmmxj = xm-xj = cos(argm)-cos(argj)= -2.*sin(argm2+argj2)*sin(argm2-argj2)
				xmmxj     = -2.*sinl(argm2+argj2)*sinl(argm2-argj2); 
			m1j *= -1;   // m1j = (-1)^j
			if(m != j){
				long double big_root = root_1pxj*root_1mxj/(root_1pxm*root_1mxm);
				D1[m][j] =  m1m*m1j*big_root/xmmxj;
			}
// 			if(m == j) D1[m][j] = -0.5*xm/xmm1_2;
			if(m != j) D1[m][m] -= D1[m][j];
		}
	}
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrix_Lobatto(int N, double **D1)
{
	long double h = Pi/N, m1m = -1.; // h: Stepsize in phi=arccosl(x)
	int m, j;
	
	for(m=0; m<=N; m++){
		long double argm = h*m, argm2 = 0.5*argm, m1j = -1., ka_m, ka_j;
		if(m == 0 || m == N) ka_m = 2.;
		else                 ka_m = 1.;
		m1m *= -1;  // m1m = (-1)^m
		D1[m][m]=0.;
		for(j=0; j<=N; j++){
			long double argj = h*j, argj2 = 0.5*argj, 
			// xmmxj = xm-xj = cos(argm)-cos(argj)= -2.*sin(argm2+argj2)*sin(argm2-argj2)
				xmmxj     = -2.*sinl(argm2+argj2)*sinl(argm2-argj2); 
			if(j == 0 || j == N) ka_j = 2.;
			else                 ka_j = 1.;
			m1j *= -1;   // m1j = (-1)^j
			if(m != j)             D1[m][j] = ka_m*m1m*m1j/(ka_j*xmmxj);
			if(m != j) D1[m][m] -= D1[m][j];
		}
	}
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrix(int N, double **D1, char *grid){
    if(strcmp( grid,"Radau_RHS") ==0)
      Chebyshev_Differentiation_Matrix_Radau_RHS(N, D1);
    else if(strcmp( grid,"Radau_LHS")==0)
      Chebyshev_Differentiation_Matrix_Radau_LHS(N, D1);
    else if(strcmp( grid,"Gauss")==0)
      Chebyshev_Differentiation_Matrix_Gauss(N, D1);
    else if(strcmp( grid,"Lobatto")==0)
      Chebyshev_Differentiation_Matrix_Lobatto(N, D1);
    else{
      printf("Error in Chebyshev_Differentiation_Matrix: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto\n grid was: %s\n", grid);
      exit(1);
    }

}
//---------------------------------------------------------------
void Chebyshev_Collocations_Derivatives(double *psi, double *d1psi, double **D1, int N)
{
	int m,j;
	long double s1, s2;
	
	for(m=0; m<=N; m++){
		s1=s2=0.;
		for(j=0; j<=m; j++)
			s1 += D1[m][j]*psi[j];
		for(j=N; j>m; j--)
			s2 += D1[m][j]*psi[j];
		d1psi[m] = s1+s2;
	}
}
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
// ********* 2D Routines ****************
//---------------------------------------------------------------
//---------------------------------------------------------------------------------------
void Chebyshev_Coefficients_2D(double **X, double **c2D, int N1, char *grid1, int N2, char *grid2)
{
	int i, j, N = maximum2(N1, N2);
	double *Y, **c_tilde;
	
	Y       = dvector(0, N);
	c_tilde = dmatrix(0, N2, 0, N1);

	for(j=0; j<=N2; j++){
		for(i=0; i<=N1; i++) Y[i] = X[i][j];
		
		Chebyshev_Coefficients(Y, c_tilde[j], N1, grid1);
	}
	for(i=0; i<=N1; i++){
		for(j=0; j<=N2; j++) Y[j] = c_tilde[j][i];
		
		Chebyshev_Coefficients(Y, c2D[i], N2, grid2);
	}
	
	free_dvector(Y, 0, N);
	free_dmatrix(c_tilde, 0, N2, 0, N1);
}
//---------------------------------------------------------------
void ChebyshevFourier_Coefficients_2D(double **X, double **alpha2D, double **beta2D, int N1, char *grid1, int N2)
{
	int i, j, N = maximum2(N1, 2*N2), n2=2*N2+1;
	double *Y, **c_tilde;
	
	Y       = dvector(0, N);
	c_tilde = dmatrix(0, 2*N2, 0, N1);

	for(j=0; j<n2; j++){
		for(i=0; i<=N1; i++) Y[i] = X[i][j];
		Chebyshev_Coefficients(Y, c_tilde[j], N1, grid1);
		
	}
	for(i=0; i<=N1; i++){
		for(j=0; j<n2; j++) Y[j] = c_tilde[j][i];
		Fourier_Coefficients(Y, alpha2D[i], beta2D[i], N2);
	}
	
	free_dvector(Y, 0, N);
	free_dmatrix(c_tilde, 0, 2*N2, 0, N1);
}
//---------------------------------------------------------------
double Clenshaw_Chebyshev_2D(double **c2D, int N1, int N2, double x1, double x2)
{
	 int i;
 	double *c_tilde, result;
 	
 	c_tilde = dvector(0, N1);
	
	for(i=0; i<=N1; i++)
 		c_tilde[i] = Clenshaw_Chebyshev(c2D[i], N2, x2);
 		
 	result = Clenshaw_Chebyshev(c_tilde, N1, x1);
 	
 	free_dvector(c_tilde, 0, N1);
 	
 	return result;
}
//---------------------------------------------------------------
double Clenshaw_ChebyshevFourier(double **alpha2D, double **beta2D, int N1, int N2, double x1, double x2)
{
	int i;
	double *c_tilde, result;
	
	c_tilde = dvector(0, N1);
	
	for(i=0; i<=N1; i++)
		c_tilde[i] = Clenshaw_Fourier(alpha2D[i], beta2D[i], N2, x2);
		
	result = Clenshaw_Chebyshev(c_tilde, N1, x1);
	
	free_dvector(c_tilde, 0, N1);
	
	return result;
}
//---------------------------------------------------------------
//---------------------------------------------------------------
// ********* 3D Routines ****************
//---------------------------------------------------------------
//---------------------------------------------------------------------------------------
void Chebyshev_Coefficients_3D(double ***X, double ***c3D, int N0, char *grid0, int N1, char *grid1, int N2, char *grid2)
{
	int i0, i1, i2, N = maximum2(N1, N2);
	double *Z, **Y, ***c_tilde;
	
	N = maximum2(N0, N);
	 
	Y  = dmatrix(0, N0, 0, N1);
	c_tilde = d3tensor(0, N2, 0, N0, 0, N1);
	
	Z = dvector(0, N2);
	
	for(i2=0; i2<=N2; i2++){
		for(i0=0; i0<=N0; i0++){
			for(i1=0; i1<=N1; i1++){
				Y[i0][i1]=X[i0][i1][i2];
			}					
		}
		Chebyshev_Coefficients_2D(Y, c_tilde[i2], N0, grid0, N1, grid1);	
	}
	
	for(i0=0; i0<=N0; i0++){
		for(i1=0; i1<=N1; i1++){
			for(i2=0; i2<=N2; i2++){
				Z[i2] = c_tilde[i2][i0][i1];
			}
		Chebyshev_Coefficients(Z, c3D[i0][i1], N2, grid2);
		}
	}
	
	free_dmatrix(Y, 0, N0, 0, N1);
	free_d3tensor(c_tilde, 0, N2, 0, N0, 0, N1);
	free_dvector(Z, 0 , N2);
		

// 	int i0, i1, i2, N = maximum2(N1, N2);
// 	double *Y, ***c_tilde0, ***c_tilde1;
// 	
// 	N = maximum2(N0, N);
// 	 
// 	Y       = dvector(0, N);
// 	c_tilde0 = d3tensor(0, N2, 0, N1, 0, N0);
// 	c_tilde1 = d3tensor(0, N2, 0, N0, 0, N1);
// 
// 	for(i2=0; i2<=N2; i2++){
// 		for(i1=0; i1<=N1; i1++){
// 			for(i0=0; i0<=N0; i0++){
// 				Y[i0] = X[i2][i1][i0];
// 			}
// 			Chebyshev_Coefficients(Y, c_tilde0[i2][i1], N0, grid0);
// 		}
// 		
// 		for(i0=0; i0<=N0; i0++){
// 			for(i1=0; i1<=N1; i1++){
// 				Y[i1] = c_tilde0[i2][i1][i0];
// 			}
// 			Chebyshev_Coefficients(Y, c_tilde1[i2][i0], N1, grid1);
// 		}			
// 	}
// 	
// 	for(i0=0; i0<=N0; i0++){
// 		for(i1=0; i1<=N1; i1++){
// 			for(i2=0; i2<=N2; i2++){
// 				Y[i2] = c_tilde1[i2][i0][i1];				
// 			}
// 			Chebyshev_Coefficients(Y, c3D[i0][i1], N2, grid2);
// 		}
// 	}
// 		
// 
// 	
// 	free_dvector(Y, 0, N);
// 	free_d3tensor(c_tilde0, 0, N2, 0, N1, 0, N0);
// 	free_d3tensor(c_tilde1, 0, N2, 0, N0, 0, N1);
}
//---------------------------------------------------------------
double Clenshaw_Chebyshev_3D(double ***c3D, int N0, int N1, int N2, double x0, double x1, double x2)
{
	int i0;
 	double *c_tilde, result;
 	
 	c_tilde = dvector(0, N0);
	
	for(i0=0; i0<=N0; i0++){
 		c_tilde[i0] = Clenshaw_Chebyshev_2D(c3D[i0], N1, N2, x1, x2);
 	}
 		
 	result = Clenshaw_Chebyshev(c_tilde, N0, x0);
 	
 	free_dvector(c_tilde, 0, N0);
 	
 	return result;
}
