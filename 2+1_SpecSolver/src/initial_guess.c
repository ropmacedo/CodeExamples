#include "2+1_Free_Boundaries.h"

void get_initial_guess(parameters par, double *X){

	printf("Calculating initial guess....\n");fflush(0);

	if(par.EXACT_SOL_CHECK == 1){
		get_ExactSolution(par, X);
		return;
	}

	if(par.SOLVER_METHOD == 0)
		fill0_dvector(X,0, par.ntotal-1);
	else if(par.SOLVER_METHOD == 1){
		SDIRK_Initial_Guess(par, X);
	}
	else if(par.SOLVER_METHOD == 2){
		SDIRK_Initial_Guess(par, X);
	}
	else{
		fprintf(stderr, "Error in get_initial_guess: Method = %d not known\n", par.SOLVER_METHOD );
	}
	printf("...Done\n");fflush(0);
	
	return;
}
//------------------------------------------------------------------
void SDIRK_Initial_Guess(parameters par, double *X)
{
	int j, iField, i_slice, nslice_total = par.nslice_total, Nslice_total=nslice_total-1, 
		N0=par.N0, i1, N1=par.N1, i2, N2=par.N2;
	double *X_previous = dvector(0,Nslice_total), *X_next = dvector(0,Nslice_total);
	
	for(i_slice=0; i_slice < nslice_total; i_slice++){
		X_previous[i_slice] = par.ID.X[i_slice];
	}
	
	
	 for(j=0; j<=N0; j++){
		//Order of vector in time direction is contrary to time flow
		int i0_previous=idx_order_time_flow(j-1, N0),
			i0_next=idx_order_time_flow(j, N0);

		double  x0_previous = j==0? -1.:par.grid_points_x0[i0_previous],
				x0_next=par.grid_points_x0[i0_next], 
				h = x0_next-x0_previous;

		if (Newton_SDIRK_verb == 1){ 
			if (par.SOLVER_METHOD == 1){printf("SDIRK Initial Guess: LU Newton Raphson Method\n");fflush(0);}
			if (par.SOLVER_METHOD == 2){printf("SDIRK Initial Guess: Finite Difference BiCGStab Newton Raphson Method\n");fflush(0);}
			printf("Solving for i0=%d, (x0=%lf)\n", i0_next, x0_next );fflush(0);
		}
		SDIRK_Step_IG(par, x0_previous, h, X_previous, X_next); 

		for(iField=0; iField<nFields; iField++){
			for(i1=0; i1<=N1; i1++){
				for(i2=0;i2<=N2;i2++){
					int Idx_Slice = Index_Slice(par, iField, i1, i2),
						Idx = Index(iField, i0_next, N0, i1, N1, i2, N2);
					X[Idx] = (X_next[Idx_Slice] - par.ID.X[Idx_Slice])/(1.+x0_next);

				}
			}
		}
		for(iField=0; iField<nFields1D; iField++){
			int Idx_Slice = Index_Slice_Fields1D(iField, par.nslice),
				Idx = Index_Fields1D(iField, par.ngrid, i0_next, N0);
				
			X[Idx] = (X_next[Idx_Slice] - par.ID.X[Idx_Slice])/(1.+x0_next);
		}

		copy_dvector(X_previous, X_next, 0, Nslice_total);


	 }
	
	free_dvector(X_previous,0,Nslice_total);
	free_dvector(X_next,0,Nslice_total);
	return;
}
// -------------------------------------------------------------------------------
void SDIRK_Step_IG(parameters par, double x0, double h, double *X_previous, double *X_next)
{
	int i_RK, s = par.SDIRK_s, i_slice, Nslice_total=par.nslice_total-1, m;
	double **k  = dmatrix(1, s, 0, Nslice_total), *Mj=dvector(0, Nslice_total);

	par.gah = par.gamma*h;

	fill0_dvector(X_next, 0, Nslice_total);	
	copy_dvector(Mj, X_previous, 0, Nslice_total);

	
	for(i_RK=1; i_RK<=s;i_RK++){
		double x0_RK = x0 + h*par.SDIRK_C[i_RK];
		if (Newton_SDIRK_verb == 1) printf("\n\t Runge Kutta Step = %d\n", i_RK);


		if (par.SOLVER_METHOD == 1)
			Newton_SDIRK_IG(par, x0_RK, Mj, k[i_RK]);
		else if (par.SOLVER_METHOD == 2)
			Newton_SDIRK_IG_bicgstab(par, x0_RK, Mj, k[i_RK]);


		fill0_dvector(Mj, 0, Nslice_total);

		for(i_slice=0; i_slice<= Nslice_total; i_slice++){
			for(m=1; m<=i_RK; m++){
				if(i_RK<s) Mj[i_slice]     += par.SDIRK_A[i_RK+1][m]*k[m][i_slice];
				else       X_next[i_slice] += par.SDIRK_B[m]     *k[m][i_slice];
			}
			if(i_RK<s)   Mj[i_slice]   = X_previous[i_slice] + h*Mj[i_slice];
			else      X_next[i_slice]  = X_previous[i_slice] + h*X_next[i_slice];
		}
	}

	
	free_dvector(Mj, 0, Nslice_total);
	free_dmatrix(k,  1, s, 0, Nslice_total);

	return;
}
// -------------------------------------------------------------------------------
int Newton_SDIRK_IG(parameters par, double x0, double *Mj, double *kj)
{
	int ntotal = par.nslice_total, Ntotal = ntotal-1, *indx, iter=0, j;
	double *F, *Dkj, **J, d, norm;
	
	F     = dvector(0, Ntotal);
	Dkj   = dvector(0, Ntotal);
	indx  = ivector(0, Ntotal);
	J     = dmatrix(0, Ntotal, 0, Ntotal);	
	
	fill0_dvector(kj, 0, Ntotal);
	
	F_of_kj_IG(par, x0, Mj, kj, F);
	
	norm = norm2(F, 0, Ntotal);
	copy_dvector(Dkj, F, 0, Ntotal); 
	
	if (Newton_SDIRK_verb == 1){
		printf("\t Initial Residual: \t |F| = %e \n", norm);
		printf("\t -------------------------------------------------------------------------\n");
	}
 	while(iter < Newton_SDIRK_itmin || (norm > Newton_SDIRK_tol  && iter < Newton_SDIRK_itmax)){
		iter += 1;

		/*************ATTENTION: DF_of_X_SDIRK_IG NOT WORKING***********/
		// DF_of_X_SDIRK_IG( par, x0, Mj, kj, /*Wslice, W1D,*/ J); 


		Jacobian_SDIRK_IG(par, x0, Mj, kj, J);
		
		//******************************************
		// ELLIPTIC LU SOLVER FOR INITIAL GUESS
		ludcmp(J, Ntotal, indx, &d, 0);
		lubksb(J, Ntotal, indx, Dkj, 0);
		//******************************************


		for (j = 0; j < ntotal; j++) 
			kj[j] -= Dkj[j];

		F_of_kj_IG(par, x0, Mj, kj, /*Wslice, W1D,*/ F);
		norm = norm2(F, 0, Ntotal);
		copy_dvector(Dkj, F, 0, Ntotal); 
		if (Newton_SDIRK_verb == 1){
// 			printf(" Newton_SDIRK: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e \n",
			if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
				printf("\n\t No Convergence of the Newton_SDIRK Raphson Method. Now exiting to system.\n\n");
// 				exit(1);
			}
			printf("\t Newton_SDIRK: iter = %3d \t |F| = %e (%e)\n", iter, norm, Newton_SDIRK_tol);
		}	
	}
	if(norm > Newton_SDIRK_tol)
		printf(
		"\t Newton_SDIRK Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);

	free_dvector(F,      0, Ntotal);
	free_dvector(Dkj,    0, Ntotal);
	free_ivector(indx,   0, Ntotal);
	free_dmatrix(J,      0, Ntotal, 0, Ntotal);


	return iter;
}
//-----------------------------------------------------------------------------
void Jacobian_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, double **J)
{
	int j, l, ntotal = par.nslice_total, Ntotal = ntotal - 1;
	double *Dkj, *JDkj;
	
	Dkj  = dvector(0, Ntotal);
	JDkj = dvector(0, Ntotal);
	
	fill0_dvector(Dkj,   0, Ntotal);
	fill0_dvector(JDkj,  0, Ntotal);

	for(j=0; j <= Ntotal; j++){
		Dkj[j] = 1.;
		J_times_Dkj_IG(par, x0, Mj, kj, Dkj, JDkj);
		for(l=0; l <= Ntotal; l++)
			J[l][j] = JDkj[l];
		Dkj[j] = 0.;
	}
	
	free_dvector(Dkj,  0, Ntotal);
	free_dvector(JDkj, 0, Ntotal);
}
// -------------------------------------------------------------------------------
void DF_of_X_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, double **J)
{
	/*************ATTENTION: DF_of_X_SDIRK_IG NOT WORKING***********/
	int j, k, ntotal = par.nslice_total, Ntotal = ntotal - 1;
	double *kjp, *Fp, *kjm, *Fm, eps = 5.e-06;
	
	kjp = dvector(0, Ntotal);
	Fp  = dvector(0, Ntotal);
	kjm = dvector(0, Ntotal);
	Fm  = dvector(0, Ntotal);
	
	for(j=0; j<= Ntotal; j++)
		kjp[j] = kjm[j] = kj[j];

	for(j=0; j<= Ntotal; j++){
		kjp[j] += eps; 
		kjm[j] -= eps; 
		F_of_kj_IG(par, x0, Mj, kjp, Fp);
		F_of_kj_IG(par, x0, Mj, kjm, Fp);
		for(k=0; k<= Ntotal; k++)
			J[k][j] = 0.5*(Fp[k]-Fm[k])/eps;
		kjp[j] = kjm[j] = kj[j];
	}
	
	free_dvector(kjp, 0, Ntotal);
	free_dvector(kjm, 0, Ntotal);
	free_dvector(Fm,  0, Ntotal);
	free_dvector(Fp,  0, Ntotal);
}
// -------------------------------------------------------------------------------
void get_Ws_from_M_and_k(parameters par, double *Mj, double *kj, derivs Wslice, derivs W1D){
	int islice, nslice=par.nslice, iField;

	for(islice=0; islice<nslice; islice ++){
		Wslice.Y[islice] = kj[islice];
		Wslice.X[islice] = Mj[islice] + par.gah*kj[islice];
	}

	for(iField=0; iField < nFields1D; iField++){
		int i_slice = Index_Slice_Fields1D(iField, nslice);
		W1D.Y[iField] = kj[i_slice];
		W1D.X[iField] = Mj[i_slice] + par.gah*kj[i_slice];
	}

	return;
}
// -------------------------------------------------------------------------------
void get_W2D_at_Grid(parameters par, int i1, int i2, derivs Wslice, derivs W){
	int iField;

	for(iField=0; iField < nFields; iField++){
				int indx_slice = Index_Slice(par, iField, i1, i2);
				
				W.X[iField]    = Wslice.X[indx_slice];  
				W.X1[iField]   = Wslice.X1[indx_slice]; 
				W.X11[iField]  = Wslice.X11[indx_slice];
				W.X2[iField]   = Wslice.X2[indx_slice]; 
				W.X12[iField]  = Wslice.X12[indx_slice];
				W.X22[iField]  = Wslice.X22[indx_slice];

				W.Y[iField]    = Wslice.Y[indx_slice];  
				W.Y1[iField]   = Wslice.Y1[indx_slice]; 
				W.Y11[iField]  = Wslice.Y11[indx_slice];
				W.Y2[iField]   = Wslice.Y2[indx_slice]; 
				W.Y12[iField]  = Wslice.Y12[indx_slice];
				W.Y22[iField]  = Wslice.Y22[indx_slice];
	}

	return;
}
// -------------------------------------------------------------------------------
void F_of_kj_IG(parameters par, double x0, double *Mj, double *kj, double *F)
{
	
	int nslice=par.nslice, iField, i1, i2, N1=par.N1, N2=par.N2, Naux = maximum2(nFields, nFields1D)-1;
	double *Eq  = dvector(0, Naux);
	
	derivs W, Wslice, W1D;	

	allocate_derivs(&W, 0, nFields-1);

	if(nFields1D>0) allocate_derivs(&W1D,    0, nFields1D-1);
	allocate_derivs(&Wslice, 0, par.nslice-1);

	get_Ws_from_M_and_k(par, Mj, kj, Wslice, W1D);
	Derivatives_Slice(par, Wslice);
	
	for(i1=0; i1<=N1; i1++){
		double x1 = par.grid_points_x1[i1];
		for(i2=0; i2 <=N2 ; i2++){
			double x2 = par.grid_points_x2[i2];	

			get_W2D_at_Grid(par, i1, i2, Wslice, W);		
			F_of_XY(par, x0, x1, x2, W, W1D, Eq);
		
			for(iField=0; iField < nFields; iField++){
				int i_slice = Index_Slice(par, iField, i1, i2);
				F[i_slice] = Eq[iField];
			}
		}
	}
	FBound_of_XY(par, x0, Wslice, W1D, Eq);
	for(iField=0; iField < nFields1D; iField++){
		int i_slice = Index_Slice_Fields1D(iField, nslice);
		F[i_slice] = Eq[iField];
	}
	
	free_dvector(Eq,   0, Naux);
	free_derivs(&W,      0, nFields-1);

	if(nFields1D>0) free_derivs(&W1D,    0, nFields1D-1);
	free_derivs(&Wslice, 0, par.nslice-1);
}
// -------------------------------------------------------------------------------
void J_times_Dkj_IG(parameters par, double x0, double *Mj, double *kj, double *Dkj, double *JDqj)
{	
	int nslice=par.nslice, nslice_total = par.nslice_total, iField, i1, i2, N1=par.N1, N2=par.N2, Naux = maximum2(nFields, nFields1D)-1;
	double *DEq, *Null;
	derivs W, DW, DW1D, DWslice, Wslice, W1D;
	
	DEq  = dvector(0, Naux);
	
	allocate_derivs(&W,       0, nFields-1);
	if(nFields1D>0) allocate_derivs(&W1D,    0, nFields1D-1);
	allocate_derivs(&Wslice, 0, par.nslice-1);

	allocate_derivs(&DW,      0, nFields-1);
	if(nFields1D>0) allocate_derivs(&DW1D,    0, nFields1D-1);
	allocate_derivs(&DWslice, 0, nslice-1);

	

	Null = dvector(0, nslice_total-1);
  	fill0_dvector(Null, 0, nslice_total-1);

	get_Ws_from_M_and_k(par, Mj, kj, Wslice, W1D);
	Derivatives_Slice(par, Wslice);

	get_Ws_from_M_and_k(par, Null, Dkj, DWslice, DW1D);
	Derivatives_Slice(par, DWslice);
	
	for(i1=0; i1<=N1; i1++){
		double x1 = par.grid_points_x1[i1];
		for(i2=0; i2 <=N2 ; i2++){
			double x2 = par.grid_points_x2[i2];
			get_W2D_at_Grid(par, i1, i2, Wslice, W);
			get_W2D_at_Grid(par, i1, i2, DWslice, DW);	

			JDX_of_XY(par, x0, x1, x2, W, W1D, DW, DW1D, DEq); // Jaux contains normal derivatives at inner boundaries
			for(iField=0; iField < nFields; iField++){
				int i_slice = Index_Slice(par, iField, i1, i2);
				JDqj[i_slice] = DEq[iField];
			}
		}
	}
	JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, DEq);
	for(iField=0; iField < nFields1D; iField++){
		int i_slice = Index_Slice_Fields1D(iField, nslice);	
		JDqj[i_slice] = DEq[iField];
	}
	
	free_dvector(Null, 0, nslice_total-1);
	free_dvector(DEq,    0, Naux);
	free_derivs(&W,       0, nFields-1);
	free_derivs(&DW,      0, nFields-1);
	if(nFields1D>0) free_derivs(&DW1D,    0, nFields1D-1);
	free_derivs(&DWslice, 0, nslice-1);

	if(nFields1D>0) free_derivs(&W1D,    0, nFields1D-1);
	free_derivs(&Wslice, 0, par.nslice-1);
}
