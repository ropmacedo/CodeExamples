#include "2+1_Free_Boundaries.h"

//----------------------------------------------------------------
int LinSolve_SDIRK_bicgstab(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj)
{	
	int	nslice_total=par.nslice_total, Nslice_total = nslice_total-1, iter=0, j, output_bicgstab=0;
	double *F, *Dkj, norm;

	F       = dvector(0, Nslice_total);   
	Dkj      = dvector(0, Nslice_total);          

	fill0_dvector(kj, 0, Nslice_total);
	F_of_kj(par, chebV, chebF, x0, Mj, kj, F);

	norm = norm2(F, 0, Nslice_total);

	fill0_dvector(Dkj, 0, Nslice_total);

	double tol=Newton_LinSolve_SDIRK_tol*norm;


	if (Newton_LinSolve_SDIRK_verb == 1){
		printf("\t\t\t\t Initial Residual: \t |F| = %e (%e)\n", norm, tol);
		printf("\t\t\t\t ------------------------------------------------------------------\n");
	}
	while(iter < Newton_LinSolve_SDIRK_itmin || (norm > tol  && iter < Newton_LinSolve_SDIRK_itmax)){
		double normres=-1.0;
		iter += 1;

					

		par.SDIRK_bicgstab_findiff_tol = SDIRK_bicgstab_findiff_decr*norm;
		output_bicgstab = SDIRK_bicgstab_findiff(par, chebV, chebF, x0, Mj, kj, F, &normres);
		

		for (j = 0; j < nslice_total; j++) kj[j] += Dkj[j];

		F_of_kj(par, chebV, chebF, x0, Mj, kj, F);
		norm = norm2(F, 0, Nslice_total);
		copy_dvector(Dkj, F, 0, Nslice_total);

		if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
			printf("\n\t\t\t\t No Convergence of the Linear Solver SDIRK. Now exiting to system.\n\n");
			exit(1);
		}
		if (Newton_LinSolve_SDIRK_verb == 1){
			printf("\t\t\t\t Finite Differece BiCGStab Linear Solver: iter = %3d \t Number of bicgstab-steps = %3d (res=%e) \t |F| = %e (%e) \n",
					iter, output_bicgstab, normres, norm, tol);
		}
	}
	if(norm > tol)
		printf(
		"\t\t Linear Solver SDIRK failed to converge to prescribed tolerance (norm = %3.15e, tol = %3.15e). Now move on to the next sequence element.\n", norm, tol);


	free_dvector(F,    0, Nslice_total);
	free_dvector(Dkj,   0, Nslice_total);
	
	return iter;
}
// -------------------------------------------------------------------------------
int SDIRK_bicgstab_findiff(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, /*double *Dkj,*/ double *F, double *normres){
	int ntotal=par.nslice_total, Ntotal = ntotal-1, iter, n;
	double alpha = 0., beta = 0., rho = 0., rho1 = 1., rhotol = 1e-50, omega = 1.0e-50, omegatol = 1e-50,
		*p, *ph, *q, *r, *rt, *s, *sh, *t;

	

	r  = dvector(0, Ntotal), rt = dvector(0, Ntotal), p  = dvector(0, Ntotal);
	// 	compute initial residual rt = r = p = F - J*DX
	// F_of_kj(par, chebV, chebF, x0, Mj, Dkj, r);
	// J_times_Dkj_IG(par, x0, Mj, kj, Dkj, r);
	fill0_dvector(r,0,Ntotal);
	for (n = 0; n < ntotal; n++) rt[n] = r[n] = p[n] = F[n] /*- r[n]*/;


	*normres = norm2(r, 0, Ntotal);
	
	

	if (*normres <= par.SDIRK_bicgstab_findiff_tol){
		free_dvector(r,  0, Ntotal);
		free_dvector(rt,  0, Ntotal);
		free_dvector(p,  0, Ntotal); 
		return 0;
	}

	ph = dvector(0, Ntotal);
	q  = dvector(0, Ntotal);    
	s  = dvector(0, Ntotal);    
	sh = dvector(0, Ntotal);
	t  = dvector(0, Ntotal);	
	
	JFD_Components JFD;	
	get_SDIRK_BandMatrix(par, chebV, chebF, x0, Mj, kj, F, &JFD);

	double resprecond;

	for (iter = 0; iter < SDIRK_bicgstab_findiff_itmax; iter++) {
		rho = scalarproduct(rt, r, 0, Ntotal);
		if (fabs(rho) < rhotol)	  break;
	
		
		if (iter > 0){ // compute direction vector p
			beta = (rho/rho1)*(alpha/omega);
			for (n = 0; n < ntotal; n++) 
				p[n] = r[n] + beta * (p[n] - omega * q[n]);
		}
		// compute direction adjusting vector ph and scalar alpha :
		
		PreCond_FinDiff(par, JFD, p, ph, &resprecond); //J_FD*ph=p (solve for ph)
		
		F_of_kj(par, chebV, chebF, x0, Mj, ph, q);
		for(n=0; n<ntotal; n++) q[n] = F[n]-q[n]; // q = J*ph
		
		alpha = rho/scalarproduct(rt, q, 0, Ntotal);
		
		for (n = 0; n < ntotal; n++) s[n] = r[n] - alpha * q[n];
		
		// early check of tolerance:
		*normres = norm2(s, 0, Ntotal);
	
		
		if (*normres <= par.SDIRK_bicgstab_findiff_tol){
			for (n = 0; n < ntotal; n++)
				kj[n] += alpha * ph[n];
				// Dkj[n] += alpha * ph[n];
			break;
		}
		// compute stabilizer vector sh and scalar omega:

		PreCond_FinDiff(par, JFD, s, sh, &resprecond);//J_FD*sh=s (solve for sh)
		
		F_of_kj(par, chebV, chebF, x0, Mj, sh, t);
		for(n=0; n<ntotal; n++) t[n] = F[n]-t[n]; // t = J*sh
		

		omega = scalarproduct(t, s, 0, Ntotal) / scalarproduct(t, t, 0, Ntotal);

		
		// compute new solution approximation:
		for (n = 0; n < ntotal; n++) {
			// if(n<ntotal) Dkj[n] += alpha * ph[n] + omega * sh[n];
			if(n<ntotal) kj[n] += alpha * ph[n] + omega * sh[n];
			r[n]   = s[n] - omega * t[n];
		}
		
		rho1 = rho;
		// are we done? 
		*normres = norm2(r, 0, Ntotal);
		if (SDIRK_bicgstab_findiff_verb == 1) {printf("\t\t\t\t\t Finite difference BiCGStab: iter = %2d  norm = %6.1e (%6.1e)\r", iter+1, *normres, par.SDIRK_bicgstab_findiff_tol);fflush(0); }
		if (*normres <= par.SDIRK_bicgstab_findiff_tol)  break;
		
		if (fabs(omega) < omegatol)  break;
	}
	

	free_bandMatrix(par,&JFD);
	

	free_dvector(p,  0, ntotal-1); 		free_dvector(ph, 0, ntotal-1); 
	free_dvector(q,  0, ntotal-1); 		free_dvector(r,  0, ntotal-1);
	free_dvector(s,  0, ntotal-1); 		free_dvector(sh, 0, ntotal-1); 
	free_dvector(rt, 0, ntotal-1); 		free_dvector(t,  0, ntotal-1);
		
	/* iteration failed */
	if (iter > SDIRK_bicgstab_findiff_itmax) return  -1;
	
	/* breakdown */
	if (fabs(rho) < rhotol) return -10;
	if (fabs(omega) < omegatol) return -11;
	
	/* success! */
	if (SDIRK_bicgstab_findiff_verb == 1) {printf("\t\t\t\t\t Finite difference BiCGStab: iter = %2d  norm = %6.1e (%6.1e)\n", iter+1, *normres, par.SDIRK_bicgstab_findiff_tol);fflush(0); }
	return ++iter;
}
//----------------------------------------------------------
void get_SDIRK_BandMatrix(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F, JFD_Components *JFD){

    int n=par.nslice, N = n-1, NFields1D = nFields1D-1,       
        maxcol = nFields*STENCILSIZE;
    
//  allocating the components of JFD
    (*JFD).J_UL         = dmatrix(0, N, 0, maxcol-1);
    if(nFields1D>0){
      (*JFD).J_UR         = dmatrix(0, N, 0, NFields1D);
      (*JFD).J_LL         = dmatrix(0, NFields1D, 0, N);
      (*JFD).J_LR         = dmatrix(0, NFields1D, 0, NFields1D);
    }
    
    (*JFD).cols_J_UL    = imatrix(0, N, 0, maxcol-1);
    (*JFD).ncols_J_UL   = ivector(0, N);
    
    fill0_ivector((*JFD).ncols_J_UL,  0, N+1);
    
    Get_SDIRK_JFD_Matrix(par, chebV, chebF, x0, Mj, kj, F, JFD);

    
    
    (*JFD).K_UL         = dmatrix(1, N+1,  1, (*JFD).m1+(*JFD).m2+1);
    (*JFD).Ku_UL         = dmatrix(1, N+1,  1, (*JFD).m1+(*JFD).m2+1);
    fill0_dmatrix((*JFD).K_UL,    1, N+1,  1, (*JFD).m1+(*JFD).m2+1);
    fill0_dmatrix((*JFD).Ku_UL,    1, N+1,  1, (*JFD).m1+(*JFD).m2+1);
    
    (*JFD).Kl_UL        = dmatrix(1, N+1,    1, (*JFD).m1);
    (*JFD).iK_UL        = ivector(1, N+1);
//     
    if(nFields1D>0){
      (*JFD).K_UR        = dmatrix(1, N+1,    1, nFields1D);
      (*JFD).K_LR        = dmatrix(1, nFields1D,    1, nFields1D); 
      (*JFD).Klu_LR        = dmatrix(1, nFields1D,    1, nFields1D);      
      (*JFD).indx_LR        = ivector(1, nFields1D); 
    }

    
    Get_JFD_Components(par, JFD);
    return;
}//----------------------------------------------------------------------------------------------------------
void Get_SDIRK_JFD_Matrix(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F, JFD_Components *JFD){ // Calculating JFD
	 int N1=par.N1, n1=par.n1, 
	     N2=par.N2, n2=par.n2,
	     ntotal=par.nslice_total, Ntotal = ntotal-1, 
	     n=par.nslice, N = n-1, NFields1D = nFields1D-1,	    
	     I, i, j, row, column, mcol, ii, jj, J, i0, i1, j0, j1;
	
	
	// double *EqPar=dvector(0,nFields1D);
	double *Fslice, *DX, *JDX, *Eq=dvector(0, nFields-1), *EqPar=dvector(0,nFields1D);
	derivs Wslice, W1D, DWslice, DW1D, W_atGrid, DW_atGrid;
	
	Fslice = dvector(0, Ntotal);
	allocate_derivs(&Wslice, 0, N);
	if(nFields1D>0) allocate_derivs(&W1D,    0, NFields1D);

	Get_Wslice_W1D_Fslice(par, x0, chebF, chebV, Wslice, W1D, Fslice);
	
	allocate_derivs(&DWslice, 0, N);
	if(nFields1D>0) allocate_derivs(&DW1D,    0, NFields1D);
 	fill0_derivs(DWslice, 0, N);	
	fill0_derivs(DW1D,0, NFields1D);
	
	allocate_derivs(&W_atGrid,  0, nFields-1);
	allocate_derivs(&DW_atGrid, 0, nFields-1);

	DX=dvector(0, Ntotal);
	fill0_dvector(DX,  0, Ntotal);


	JDX = dvector(0, Ntotal);
	fill0_dvector(JDX,  0, Ntotal);

	(*JFD).m1 = (*JFD).m2 = 0;



		for(column = 0; column < n; column ++){
			get_indices_from_Index_Slice(par, column, &I, &i, &j);
			double x1=par.grid_points_x1[i],
				   x2=par.grid_points_x2[j];

			DX[column]=1.;
			get_Ws_from_M_and_k(par, Mj, DX, DWslice, DW1D);

			int d = FD_ORDER/2 + 1;
			i0 = maximum2(0,i-d);
			i1 = minimum2(i+d,n1-1);
			j0 = maximum2(0,j-d);
			j1 = minimum2(j+d,n2-1);  

			for(ii=i0; ii<=i1; ii++) Get_DerivativesFinDif_grid(par, I, ii,  j, DWslice);
			for(jj=j0; jj<=j1; jj++) Get_DerivativesFinDif_grid(par, I,  i, jj, DWslice); 

			for(jj=j0; jj<=j1; jj++){
				get_W2D_at_Grid(par, i,  jj, Wslice, W_atGrid);
				get_W2D_at_Grid(par, i,  jj, DWslice, DW_atGrid);
				JDX_of_XY(par, x0, x1, x2, W_atGrid, W1D, DW_atGrid, DW1D, Eq);


				for(J=0;J<nFields; J++){
					row = Index_Slice(par, J, i, jj);
					JDX[row]=Eq[J] - Fslice[row];

					if (fabs(JDX[row]) > TINY){
						mcol = (*JFD).ncols_J_UL[row];
						(*JFD).cols_J_UL[row][mcol] =  column;
						(*JFD).J_UL[row][mcol]      =  JDX[row];
						(*JFD).ncols_J_UL[row]     +=  1;

						if(column > row) // determining the numbers m1 and m2
							(*JFD).m2=maximum2((*JFD).m2, column-row);
						else
							(*JFD).m1=maximum2((*JFD).m1, row-column);
					}
					JDX[row] =  0.;
				}
			}

			for(ii=i0; ii<=i1; ii++){
				get_W2D_at_Grid(par, ii, j, Wslice, W_atGrid);
				get_W2D_at_Grid(par, ii, j, DWslice, DW_atGrid);

				JDX_of_XY(par, x0, x1, x2, W_atGrid, W1D, DW_atGrid, DW1D, Eq);

				for(J=0;J<nFields; J++){
					row = Index_Slice(par,J,ii,j);
					JDX[row]=Eq[J]- Fslice[row];

					if (fabs(JDX[row]) > TINY){
						mcol = (*JFD).ncols_J_UL[row];
						(*JFD).cols_J_UL[row][mcol] =  column;
						(*JFD).J_UL[row][mcol]      =  JDX[row];
						(*JFD).ncols_J_UL[row]     +=  1;

						if(column > row) // determining the numbers m1 and m2
							(*JFD).m2=maximum2((*JFD).m2, column-row);
						else
							(*JFD).m1=maximum2((*JFD).m1, row-column);
					}
					JDX[row] =  0.;
				}
			}

			JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, EqPar);
			for(J=0; J<nFields1D; J++){
				int row = Index_Slice_Fields1D(J, n);
				(*JFD).J_LL[J][column]=EqPar[J] - Fslice[row];
			}
			fill0_dvector(DX,  0, Ntotal);
		}

		for(column = 0; column < nFields1D; column ++){

			int column_1D = Index_Slice_Fields1D(column, n);
			DX[column_1D]=1.;
			get_Ws_from_M_and_k(par, Mj, DX, DWslice, DW1D);
		// 	fill0_derivs(DWslice, 0, N);
		

		for(jj=0; jj<=N2; jj++){
			double x2=par.grid_points_x2[jj];
			for(ii=0; ii<=N1; ii++){
				double x1=par.grid_points_x1[ii];

				get_W2D_at_Grid(par, ii, jj, Wslice, W_atGrid);
				get_W2D_at_Grid(par, ii, jj, DWslice, DW_atGrid);

				JDX_of_XY(par, x0, x1, x2, W_atGrid, W1D, DW_atGrid, DW1D, Eq);

				for(J=0; J<nFields; J++){
					int I = Index_Slice(par, J,ii,jj);
					(*JFD).J_UR[I][column] = Eq[J]- Fslice[I];
				}
			}
		}

		JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, EqPar);
		for(J=0; J<nFields1D; J++){
			int row = Index_Slice_Fields1D(J, n);
			(*JFD).J_LR[J][column] = EqPar[J]-Fslice[row];
		}
		fill0_dvector(DX,  0, Ntotal);
	}
	
	
	free_dvector(Fslice, 0, Ntotal);
	free_derivs(&Wslice,  0, N);
	free_derivs(&DWslice, 0, N);

	if(nFields1D>0) free_derivs(&W1D,  0, NFields1D);
	if(nFields1D>0)free_derivs(&DW1D, 0, NFields1D);

	free_derivs(&W_atGrid,  0, nFields-1);
	free_derivs(&DW_atGrid, 0, nFields-1);

	free_dvector(DX, 0, Ntotal);
	free_dvector(JDX, 0, Ntotal);
	free_dvector(Eq, 0, nFields-1);
	free_dvector(EqPar, 0, nFields1D);

	return;
}



