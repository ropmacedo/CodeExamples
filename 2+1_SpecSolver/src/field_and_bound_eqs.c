#include "2+1_Free_Boundaries.h"


void F_of_XY(parameters par, double x0, double x1, double x2, derivs W, derivs W1D, double *F){
	
	double t=func_t_of_x0(par, x0), r=func_r_of_x1(par, x1), x=func_x_of_x2(par, x2), kappa=par.kappa;
	
	 
	
	derivs W_trx, W1D_t;
	allocate_derivs(&W_trx, 0, nFields-1);
  	if(nFields1D!=0) allocate_derivs(&W1D_t, 0, nFields1D-1);
	
	scale_derivs_x0x1x2_to_trx(par, W, W1D, W_trx, W1D_t, nFields, nFields1D);
	
	double S, L1_0, L1_r, L1_x, L2_0, L2_r, L2_rr, L2_x, L2_xx, L2_rx, Source , l=1.*par.ell;

	get_LinearEquationCoefficients(kappa, t, r, x, &S, &L1_0, &L1_r, &L1_x,&L2_0, &L2_r, &L2_x, &L2_rr, &L2_rx, &L2_xx);
	
	Source = get_Source_f2(par, t, r, x,  W1D_t )* plgndr(l, 0, x);
	
	
	double f2 = W_trx.X[0],  f2_r = W_trx.X1[0],  f2_x = W_trx.X2[0],	f2_t = W_trx.Y[0],   f2_rx = W_trx.X11[0], f2_rr = W_trx.X11[0], f2_xx = W_trx.X22[0],
		  df2 = W_trx.X[1], df2_r = W_trx.X1[1], df2_x = W_trx.X2[1],	df2_t = W_trx.Y[1];
	
		
	F[0] = df2-f2_t;	
	F[1] = S*df2_t  + L1_r*df2_r + L1_x*df2_x + L1_0*df2
		 + L2_rr*f2_rr+ L2_r*f2_r + L2_x*f2_x + L2_xx*f2_xx + L2_rx*f2_rx + L2_0*f2 
		 + Source;  

	free_derivs(&W_trx,      0, nFields-1);
	if(nFields1D!=0) free_derivs(&W1D_t,    0, nFields1D-1);
	
	return;
}
//----------------------------------------------------------------------------------------------
void JDX_of_XY(parameters par, double x0, double x1, double x2, derivs W, derivs W1D, derivs DW, derivs DW1D, double *JDX){
	
	double t=func_t_of_x0(par, x0), r=func_r_of_x1(par, x1), x=func_x_of_x2(par, x2), kappa=par.kappa;
	

	
	derivs DW_trx, DW1D_t;
	allocate_derivs(&DW_trx, 0, nFields-1);
  	if(nFields1D!=0) allocate_derivs(&DW1D_t, 0, nFields1D-1);
	
	scale_derivs_x0x1x2_to_trx(par, DW, DW1D, DW_trx, DW1D_t, nFields, nFields1D);
	
	double S, L1_0, L1_r, L1_x, L2_0, L2_r, L2_rr, L2_x, L2_xx, L2_rx, DSource, l=1.*par.ell;

	get_LinearEquationCoefficients(kappa, t, r, x, &S, &L1_0, &L1_r, &L1_x,&L2_0, &L2_r, &L2_x, &L2_rr, &L2_rx, &L2_xx);
	
	par.alpha=0;
	DSource = get_Source_f2(par, t, r, x,  DW1D_t )* plgndr(l, 0, x);
	
	double Df2 = DW_trx.X[0],  Df2_r = DW_trx.X1[0],  Df2_x = DW_trx.X2[0],	Df2_t = DW_trx.Y[0],   Df2_rx = DW_trx.X11[0], Df2_rr = DW_trx.X11[0], Df2_xx = DW_trx.X22[0],
		  dDf2 = DW_trx.X[1], dDf2_r = DW_trx.X1[1], dDf2_x = DW_trx.X2[1],	dDf2_t = DW_trx.Y[1];
	
		
	JDX[0] = dDf2-Df2_t;	
	JDX[1] = S*dDf2_t  + L1_r*dDf2_r + L1_x*dDf2_x + L1_0*dDf2
		 + L2_rr*Df2_rr+ L2_r*Df2_r + L2_x*Df2_x + L2_xx*Df2_xx + L2_rx*Df2_rx + L2_0*Df2 
		 + DSource; 
		   

	free_derivs(&DW_trx,      0, nFields-1);
	if(nFields1D!=0) free_derivs(&DW1D_t,    0, nFields1D-1);
	
	return;
}
//----------------------------------------------------------------------------------
void FBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, double *F){

	double t=func_t_of_x0(par, x0), t2=sqr(t);
	
	derivs Wslice_trx, W1D_t;
	allocate_derivs(&Wslice_trx, 0, par.nslice-1);
  	if(nFields1D!=0) allocate_derivs(&W1D_t, 0, nFields1D-1);
  	
  	scale_derivs_x0x1x2_to_trx(par, Wslice, W1D, Wslice_trx, W1D_t, par.nslice, nFields1D);
	
	double f1 = W1D_t.X[0], f1_t = W1D_t.Y[0],
		  df1 = W1D_t.X[1], df1_t = W1D_t.Y[1],
		  Source_f1, l=1.*par.ell;
		
	Source_f1 = get_Source_f1(par, t);
		   
	F[0] = df1-f1_t;
	if(l==0 && t==1.){
			F[1] = 2*df1_t + 2.*df1 + Source_f1;
	}
	else{
			F[1] = -(1.-t2)*df1_t - 2.*(1.-t)*df1 - l*(l+1)*f1 + Source_f1;		
	}

	
	free_derivs(&Wslice_trx, 0, par.nslice-1);
	if(nFields1D!=0) free_derivs(&W1D_t,    0, nFields1D-1);
	return;
}
// -------------------------------------------------------------------------------
void JDXBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, derivs DWslice, derivs DW1D, double *JDX){

	double t=func_t_of_x0(par, x0), t2=sqr(t);
	
	derivs DWslice_trx, DW1D_t;
	allocate_derivs(&DWslice_trx, 0, par.nslice-1);
  	if(nFields1D!=0) allocate_derivs(&DW1D_t, 0, nFields1D-1);
  	
  	scale_derivs_x0x1x2_to_trx(par, DWslice, DW1D, DWslice_trx, DW1D_t, par.nslice, nFields1D);
	
	double Df1 = DW1D_t.X[0], Df1_t = DW1D_t.Y[0],
		  dDf1 = DW1D_t.X[1], dDf1_t = DW1D_t.Y[1],
		  l=1.*par.ell;
		
			   
	JDX[0] = dDf1-Df1_t;
	if(l==0 && t==1){
			JDX[1] = 2*dDf1_t + 2.*dDf1 ;
	}
	else{
			JDX[1] = -(1.-t2)*dDf1_t - 2.*(1.-t)*dDf1 - l*(l+1)*Df1;		
	}
	
	free_derivs(&DWslice_trx, 0, par.nslice-1);
	if(nFields1D!=0) free_derivs(&DW1D_t,    0, nFields1D-1);
	
	return;
}
