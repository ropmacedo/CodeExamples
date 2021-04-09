#include "2+1_Free_Boundaries.h"

//-----------------------------------------------------------------------------
void Jacobian(parameters par, double *X, double **J)
{
	int j, l, ntotal = par.ntotal, Ntotal=ntotal-1;
	double *DX, *JDX;

	DX  = dvector(0, Ntotal);
	JDX = dvector(0, Ntotal);
	
	fill0_dvector(DX,   0, Ntotal);
	fill0_dvector(JDX,  0, Ntotal);

	for(j=0; j <= Ntotal; j++){
		DX[j] = 1.;
		J_times_DX(par, X, DX, JDX);
		for(l=0; l <= Ntotal; l++)
			J[l][j] = JDX[l];
		DX[j] = 0.;
	}
	
	free_dvector(DX,  0, Ntotal);
	free_dvector(JDX, 0, Ntotal);
}
// -------------------------------------------------------------------------------
void DF_of_X(parameters par, double *X, double **J)
{
	int j, k, ntotal = par.ntotal, Ntotal=ntotal-1;
	double *Xp, *Fp, *Xm, *Fm, eps = 5.e-06;
	
	Xp = dvector(0, Ntotal);
	Fp = dvector(0, Ntotal);
	Xm = dvector(0, Ntotal);
	Fm = dvector(0, Ntotal);
	
	for(j=0; j<= Ntotal; j++)
		Xp[j] = Xm[j] = X[j];

	for(j=0; j<= Ntotal; j++){
		Xp[j] += eps; 
		Xm[j] -= eps; 
		F_of_X(par, Xp, Fp);
		F_of_X(par, Xm, Fm);
		for(k=0; k<= Ntotal; k++)
			J[k][j] = 0.5*(Fp[k]-Fm[k])/eps;
		Xp[j] = Xm[j] = X[j];
	}
	
	free_dvector(Xp, 0, Ntotal);
	free_dvector(Xm, 0, Ntotal);
	free_dvector(Fm, 0, Ntotal);
	free_dvector(Fp, 0, Ntotal);
}
// -------------------------------------------------------------------------------
int newton_direct(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int ntotal = par.ntotal, Ntotal=ntotal-1, *indx, iter=0, j, FLAG=-1;
	double *F, *DX, *J[ntotal], d, norm;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);
	indx  = ivector(0, Ntotal);
	for(j=0; j<ntotal; j++)	J[j]=dvector(0,Ntotal);
	// J     = dmatrix(0, Ntotal, 0, Ntotal);
		
	F_of_X(par, X, F);
	// PrintVector(stdout, F, 0, Ntotal, 1.e-20);

	norm = norm1(F, 0, Ntotal);	
	copy_dvector(DX, F, 0, Ntotal); 

	double tol=Newton_tol;//norm*Newton_tol;

	if (Newton_verb == 1){
		printf(" Newton Raphson Method Main Equation: LU decomposition\n");
		printf(" Initial Residual: \t |F| = %e (%e)\n", norm, tol);
		printf(" ------------------------------------------------------------------\n");
	}

//  	while(iter < Newton_itmin || (norm > tol  && iter < Newton_itmax)){
	while(iter < Newton_itmin || (FLAG < 0  && iter < Newton_itmax)){
		if(norm < tol) FLAG++;
		iter += 1;

		Jacobian(par, X, J);		
// 		DF_of_X( par, X, J);

		ludcmp(J, Ntotal, indx, &d, 0);
		lubksb(J, Ntotal, indx, DX, 0);
		
		for (j = 0; j < ntotal; j++) 
			X[j] -= DX[j];

		F_of_X(par, X, F);
		norm = norm2(F, 0, Ntotal);
		copy_dvector(DX, F, 0, Ntotal); 
		
		if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
			printf("\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n");
			exit(1);
		}
		if (Newton_verb == 1){

			printf(" Newton: iter = %3d \t |F| = %e (%3.2e), %d \n", iter, norm, tol, FLAG);
		}	
	}
	if(norm > tol)
		printf(
		" Newton Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);

	free_dvector(F,    0, Ntotal);
	free_dvector(DX,   0, Ntotal);
	free_ivector(indx, 0, Ntotal);
	// free_dmatrix(J,    0, Ntotal, 0, Ntotal);
	for(j=0; j<ntotal; j++) free_dvector(J[j],   0, Ntotal);


	return iter;
}
//-----------------------------------------------------------------------------
