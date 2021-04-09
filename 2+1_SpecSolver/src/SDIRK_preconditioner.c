#include "2+1_Free_Boundaries.h"

int idx_order_time_flow(int j, int N0){
	return N0-j;
}
// -------------------------------------------------------------------------------
void get_chebF(parameters par, double *F, double **chebF)
{
	
	int i0, N0=par.N0 , i1, N1=par.N1, i2, N2=par.N2, iField;
	double *p=dvector(0,N0), *cp=dvector(0,N0);

	
	for(iField=0; iField<nFields; iField++){
		for(i1=0; i1<=N1; i1++){
			for(i2=0;i2<=N2;i2++){
				for(i0=0; i0<=N0; i0++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					p[i0] = F[Idx];
				}
				
				Chebyshev_Coefficients(p, cp, N0, par.grid_x0);

				int Idx_slice=Index_Slice(par, iField, i1, i2);
				copy_dvector(chebF[Idx_slice], cp, 0, N0);
			}
		}
	}

	for(iField=0; iField<nFields1D; iField++){
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_Fields1D(iField, par.ngrid, i0, N0);
			p[i0] = F[Idx];
		}
		Chebyshev_Coefficients(p, cp, N0, par.grid_x0);

		int Idx_slice=Index_Slice_Fields1D(iField, par.nslice);
		copy_dvector(chebF[Idx_slice], cp, 0, N0);
	}


	free_dvector(p, 0, N0);
	free_dvector(cp, 0, N0);  
}
// -------------------------------------------------------------------------------
void PreCond(parameters par, cheb_derivs chebV, double *F,  double *X)
{
	int j, iField, N0=par.N0, i1, N1=par.N1, i2, N2=par.N2, Nslice_total=par.nslice_total-1;
	double *X_previous = dvector(0,Nslice_total), *X_next = dvector(0,Nslice_total), **chebF=dmatrix(0,Nslice_total,0,N0);

	fill0_dvector(X_previous,   0, Nslice_total);
	fill0_dvector(X_next,   0, Nslice_total);
	get_chebF(par, F, chebF);

	for(j=0; j<=N0; j++){
		//Order of vector in time direction is contrary to time flow
		int i0_previous=idx_order_time_flow(j-1, N0),
			i0_next=idx_order_time_flow(j, N0); 

		double  x0_previous = j==0? -1.:par.grid_points_x0[i0_previous],
				x0_next=par.grid_points_x0[i0_next], 
				h = x0_next-x0_previous;

		if (Newton_LinSolve_SDIRK_verb == 1){ 
			if (par.SOLVER_METHOD == 1){printf("\t\t\t SDIRK Linear Solver: LU decomposition\n");fflush(0);}
			if (par.SOLVER_METHOD == 2){printf("\t\t\t SDIRK Linear Solver: Finite Difference BiCGStab\n");fflush(0);}
			printf("\t\t\t Solving for i0=%d, (x0=%lf)\n", i0_next, x0_next );fflush(0);
		}
		SDIRK_Step(par, chebV, chebF, x0_previous, h, X_previous, X_next);

		for(iField=0; iField<nFields; iField++){
			for(i1=0; i1<=N1; i1++){
				for(i2=0;i2<=N2;i2++){
					int Idx_Slice = Index_Slice(par, iField, i1, i2),
						Idx = Index(iField, i0_next, N0, i1, N1, i2, N2);
					X[Idx] = X_next[Idx_Slice]/(1.+x0_next);

				}
			}
		}
		for(iField=0; iField<nFields1D; iField++){
			int Idx_Slice = Index_Slice_Fields1D(iField, par.nslice),
				Idx = Index_Fields1D(iField, par.ngrid, i0_next, N0);
				
			X[Idx] = X_next[Idx_Slice]/(1.+x0_next);
		}

		copy_dvector(X_previous, X_next, 0, Nslice_total);
		// pause();
	}

	free_dvector(X_previous,0,Nslice_total);
	free_dvector(X_next,0,Nslice_total);
	free_dmatrix(chebF, 0, Nslice_total, 0, N0);

	return;
}
//-------------------------------------------------------------------------
void SDIRK_Step(parameters par, cheb_derivs chebV, double **chebF, double x0, double h, double *X_previous, double *X_next){

	int i_RK, s = par.SDIRK_s, i_slice, Nslice_total=par.nslice_total-1, m;

	double **k  = dmatrix(1, s, 0, Nslice_total), *Mj=dvector(0, Nslice_total);

	par.gah = par.gamma*h;

	fill0_dvector(X_next, 0, Nslice_total);	
	copy_dvector(Mj, X_previous, 0, Nslice_total);

	for(i_RK=1; i_RK<=s;i_RK++){
	
		double x0_RK = x0 + h*par.SDIRK_C[i_RK];
		if (Newton_LinSolve_SDIRK_verb == 1) printf("\n\t\t\t\t Runge Kutta Step = %d\n", i_RK);

		if(par.SOLVER_METHOD == 1){
			LinSolve_SDIRK(par, chebV, chebF, x0_RK, Mj, k[i_RK]);
		}
		else if(par.SOLVER_METHOD == 2){
			LinSolve_SDIRK_bicgstab(par, chebV, chebF, x0_RK, Mj, k[i_RK]);
		}

		// LinSolve_SDIRK(par, chebV, chebF, x0_RK, Mj, k[i_RK]);
		// Newton_LinSolve_SDIRK(par, chebV, chebF, x0_RK, Mj, k[i_RK]);
		// Newton_LinSolve_SDIRK_bicgstab(par, chebV, chebF, x0_RK, Mj, k[i_RK]);
		


		fill0_dvector(Mj, 0, Nslice_total);
		for(i_slice=0; i_slice<= Nslice_total; i_slice++){
			for(m=1; m<=i_RK; m++){
				if(i_RK<s) Mj[i_slice]     += par.SDIRK_A[i_RK+1][m]*k[m][i_slice];
				else       X_next[i_slice] += par.SDIRK_B[m]     *k[m][i_slice];
			}
			if(i_RK<s)   Mj[i_slice]      = X_previous[i_slice] + h*Mj[i_slice];
			else      X_next[i_slice]  = X_previous[i_slice] + h*X_next[i_slice];
		}
	}
	free_dvector(Mj, 0, Nslice_total);
	free_dmatrix(k,  1, s, 0, Nslice_total);
	return;
}
//----------------------------------------------------------------
void Jacobian_LinSolve_SDIRK(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double  **J){
	int j, l, nslice_total = par.nslice_total, Nslice_total=nslice_total-1;

	double *Dk  = dvector(0, Nslice_total);
	double *JDk = dvector(0, Nslice_total);
	
	fill0_dvector(Dk,  0, Nslice_total);
	fill0_dvector(JDk,  0, Nslice_total);

	
	for(j=0; j <= Nslice_total; j++){
		Dk[j] = 1.;
		F_of_kj(par, chebV, chebF, x0, Mj, Dk, JDk);
		for(l=0; l <= Nslice_total; l++)
			J[l][j] = JDk[l]; 
		Dk[j] = 0.;
	}

	free_dvector(Dk,  0, Nslice_total);
	free_dvector(JDk, 0, Nslice_total);
	return;
}
//----------------------------------------------------------------
int Newton_LinSolve_SDIRK(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj)
{	
	int	nslice_total=par.nslice_total, Nslice_total = nslice_total-1, *indx, iter=0, j;
	double *F, *Dk, **J, d, norm;

	F       = dvector(0, Nslice_total);   
	Dk      = dvector(0, Nslice_total);          
	indx    = ivector(0, Nslice_total);
	J       = dmatrix(0, Nslice_total, 0, Nslice_total); 

	fill0_dvector(kj, 0, Nslice_total);
	F_of_kj(par, chebV, chebF, x0, Mj, kj, F);

	norm = norm2(F, 0, Nslice_total);
	copy_dvector(Dk, F, 0, Nslice_total);

	if (Newton_LinSolve_SDIRK_verb == 1){
		printf(" \tLinear Solver SDIRK: Initial Residual: \t |F| = %e \n", norm);
		printf(" \t------------------------------------------------------------------\n");
		// exit(-1);
	}

	double tol=Newton_LinSolve_SDIRK_tol*norm;
	while(iter < Newton_LinSolve_SDIRK_itmin || (norm > tol  && iter < Newton_LinSolve_SDIRK_itmax)){
		iter += 1;
		
		Jacobian_LinSolve_SDIRK(par, chebV, chebF, x0, Mj, kj, J);

		
		ludcmp(J, Nslice_total, indx, &d, 0);
		lubksb(J, Nslice_total, indx, Dk, 0);

		for (j = 0; j < nslice_total; j++) kj[j] -= Dk[j];

		F_of_kj(par, chebV, chebF, x0, Mj, kj, F);
		norm = norm2(F, 0, Nslice_total);
		copy_dvector(Dk, F, 0, Nslice_total);

		if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
			printf("\n\t No Convergence of the Linear Solver SDIRK. Now exiting to system.\n\n");
			exit(1);
		}
		if (Newton_LinSolve_SDIRK_verb == 1){
			printf("\t Linear Solver SDIRK: iter = %3d \t |F| = %e (%3.2e)\n", iter, norm, tol);
		}
	}
	if(norm > tol)
		printf(
		" Linear Solver SDIRK failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);


	free_dvector(F,    0, Nslice_total);
	free_dvector(Dk,   0, Nslice_total);
	free_ivector(indx, 0, Nslice_total);
	free_dmatrix(J,    0, Nslice_total, 0, Nslice_total);

	return iter;
}
//-------------------------------------------------------------------------
void LinSolve_SDIRK(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj){
	int	nslice_total=par.nslice_total, Nslice_total = nslice_total-1, 
		*indx, j, l;
	double *F0, **K, d, *F0_m_Kkj, norm=0., tol=0.;

	F0       = dvector(0, Nslice_total);             // F(qj) = F0 - K * qj
	F0_m_Kkj = dvector(0, Nslice_total);
	indx     = ivector(0, Nslice_total);
	K        = dmatrix(0, Nslice_total, 0, Nslice_total); 


	fill0_dvector(kj, 0, Nslice_total);
	F_of_kj(par, chebV, chebF, x0, Mj, kj, F0);

	if (Newton_LinSolve_SDIRK_verb == 1){
		norm = norm2(F0, 0, Nslice_total);
		tol=norm*Newton_LinSolve_SDIRK_tol;
		printf("\t\t\t\t Initial Residual: \t |F| = %e (%e)\n", norm,tol );
		printf("\t\t\t\t  ------------------------------------------------------------------\n");
	}


	for(j=0; j <= Nslice_total; j++){
		kj[j] = 1.;
		F_of_kj(par, chebV, chebF, x0, Mj, kj, F0_m_Kkj);
		for(l=0; l <= Nslice_total; l++)
			K[l][j] = F0[l] - F0_m_Kkj[l]; 
		kj[j] = 0.;
	}
	
	//********ELIPTICAL LU Solver*******
	ludcmp(K, Nslice_total, indx, &d, 0);
	copy_dvector(kj, F0, 0, Nslice_total);
	lubksb(K, Nslice_total, indx, kj, 0);
	//---***************************----
	

	if (Newton_LinSolve_SDIRK_verb == 1){
		F_of_kj(par, chebV, chebF, x0, Mj, kj, F0);
		norm = norm2(F0, 0, Nslice_total);
		printf("\t\t\t\t LU linear solver: |F| = %e (%3.2e) \n", norm, tol);
	}

	free_dvector(F0,     0, Nslice_total);
	free_ivector(indx,   0, Nslice_total);
	free_dmatrix(K,      0, Nslice_total, 0, Nslice_total);
}
// -------------------------------------------------------------------------------
void Get_Wslice_W1D_Fslice(parameters par, double x0, double **chebF, cheb_derivs chebV, derivs Wslice, derivs W1D, double *Fslice)
{
	int i_slice, iField, nslice=par.nslice, N0=par.N0;
	double x0p1 = x0+1.;
	
	for(i_slice=0; i_slice < nslice; i_slice ++){
		double 
			VX    = Clenshaw_Chebyshev(chebV.X[i_slice],   N0, x0),
			VX1   = Clenshaw_Chebyshev(chebV.X1[i_slice],  N0, x0),
			VX11  = Clenshaw_Chebyshev(chebV.X11[i_slice], N0, x0),
			VX2   = Clenshaw_Chebyshev(chebV.X2[i_slice],  N0, x0),
			VX12  = Clenshaw_Chebyshev(chebV.X12[i_slice], N0, x0),
			VX22  = Clenshaw_Chebyshev(chebV.X22[i_slice], N0, x0),
			VY    = Clenshaw_Chebyshev(chebV.Y[i_slice],   N0, x0),
			VY1   = Clenshaw_Chebyshev(chebV.Y1[i_slice],  N0, x0),
			VY11  = Clenshaw_Chebyshev(chebV.Y11[i_slice], N0, x0),
			VY2   = Clenshaw_Chebyshev(chebV.Y2[i_slice],  N0, x0),
			VY12  = Clenshaw_Chebyshev(chebV.Y12[i_slice], N0, x0),
			VY22  = Clenshaw_Chebyshev(chebV.Y22[i_slice], N0, x0);

		Wslice.X[i_slice]    = par.ID.X[i_slice]   + x0p1*VX;
		Wslice.X1[i_slice]   = par.ID.X1[i_slice]  + x0p1*VX1;
		Wslice.X11[i_slice]  = par.ID.X11[i_slice] + x0p1*VX11;
		Wslice.X1[i_slice]   = par.ID.X2[i_slice]  + x0p1*VX2;
		Wslice.X12[i_slice]  = par.ID.X12[i_slice] + x0p1*VX12;
		Wslice.X22[i_slice]  = par.ID.X22[i_slice] + x0p1*VX22;
		Wslice.Y[i_slice]    = VX                  + x0p1*VY;
		Wslice.Y1[i_slice]   = VX1                 + x0p1*VY1;
		Wslice.Y11[i_slice]  = VX11                + x0p1*VY11;
		Wslice.Y2[i_slice]   = VX2                 + x0p1*VY2;
		Wslice.Y12[i_slice]  = VX12                + x0p1*VY12;
		Wslice.Y22[i_slice]  = VX22                + x0p1*VY22;
		
		Fslice[i_slice] = Clenshaw_Chebyshev(chebF[i_slice], N0, x0);
		
	}
	for(iField=0; iField < nFields1D; iField++){
		i_slice = Index_Slice_Fields1D(iField, nslice);
		double 
			VX   = Clenshaw_Chebyshev(chebV.X[i_slice],   N0, x0),
			VY   = Clenshaw_Chebyshev(chebV.Y[i_slice],   N0, x0);
			
		W1D.X[iField]    = par.ID.X[i_slice] + x0p1*VX;
		W1D.Y[iField]    = VX                + x0p1*VY;

		Fslice[i_slice] = Clenshaw_Chebyshev(chebF[i_slice], N0, x0);
	}
}
//-------------------------------------------------------------------------
void F_of_kj(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F){

	int nslice=par.nslice, Nslice = nslice-1, Nslice_total=par.nslice_total -1, iField, i1, i2, N1=par.N1, N2=par.N2, Naux = maximum2(nFields, nFields1D)-1;
	double *Faux   = dvector(0, Naux), *Fslice = dvector(0, Nslice_total);

	derivs DW, W1D, Wslice, W, DWslice, DW1D;
	allocate_derivs(&DW, 0, nFields-1);
	allocate_derivs(&W,  0, nFields-1);
	if(nFields1D>0) allocate_derivs(&W1D,    0, nFields1D-1);
	allocate_derivs(&Wslice, 0, Nslice);

	allocate_derivs(&DWslice, 0, Nslice);
	if(nFields1D>0) allocate_derivs(&DW1D,    0, nFields1D-1);


	Get_Wslice_W1D_Fslice(par, x0, chebF, chebV, Wslice, W1D, Fslice);

	get_Ws_from_M_and_k(par, Mj, kj, DWslice, DW1D);
	Derivatives_Slice(par, DWslice);

	
	for(i1=0; i1<=N1; i1++){
		double x1 = par.grid_points_x1[i1];
		for(i2=0; i2 <=N2 ; i2++){
			double x2 = par.grid_points_x2[i2];	

			get_W2D_at_Grid(par, i1, i2, Wslice, W);
			get_W2D_at_Grid(par, i1, i2, DWslice, DW);

			JDX_of_XY(par, x0, x1, x2, W, W1D, DW, DW1D, Faux);
			for(iField=0; iField < nFields; iField++){
				int i_slice = Index_Slice(par, iField, i1, i2);
				F[i_slice] = Faux[iField] - Fslice[i_slice];
			}	
		}		
	}
	JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, Faux);
	for(iField=0; iField < nFields1D; iField++){
		int i_slice = Index_Slice_Fields1D(iField, nslice);
		F[i_slice] = Faux[iField] - Fslice[i_slice];
	}
	
	
	free_dvector(Faux,   0, Naux);
 	free_dvector(Fslice, 0, Nslice_total);

	free_derivs(&DW,     0, nFields-1);
	free_derivs(&W,      0, nFields-1);
	if(nFields1D>0) free_derivs(&W1D,    0, nFields1D-1); 
	free_derivs(&Wslice, 0, Nslice); 
	
	free_derivs(&DWslice, 0, Nslice); 
	if(nFields1D>0) free_derivs(&DW1D,    0, nFields1D-1); 

	return;
}

