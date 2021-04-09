#include "2+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
int Newton_SDIRK_IG_bicgstab(parameters par, double x0, double *Mj, double *kj)
{
	int ntotal = par.nslice_total, Ntotal = ntotal-1, iter=0, j, output_bicgstab=0;
	double *F, *Dkj, norm;
	
	F     = dvector(0, Ntotal);
	Dkj   = dvector(0, Ntotal);
	
	fill0_dvector(kj, 0, Ntotal);	
	F_of_kj_IG(par, x0, Mj, kj, F);
	
	norm = norm2(F, 0, Ntotal);
	
	if (Newton_SDIRK_verb == 1){
		printf("\t Initial Residual: \t |F| = %e \n", norm);
		printf("\t -------------------------------------------------------------------------\n");
	}

	while(iter < Newton_SDIRK_itmin || (norm > Newton_SDIRK_tol  && iter < Newton_SDIRK_itmax)){
		double normres=-1.0;
		iter += 1;

		fill0_dvector(Dkj, 0, Ntotal);

		par.IG_bicgstab_findiff_tol = IG_bicgstab_findiff_decr*norm;
		output_bicgstab = bicgstab_findiff(par, x0, Mj, kj, Dkj, F, &normres);
		

		for (j = 0; j < ntotal; j++) 
			kj[j] -= Dkj[j];

		F_of_kj_IG(par, x0, Mj, kj, F);
		norm = norm2(F, 0, Ntotal);
		if (isinf(norm) || isnan(norm) || norm > 1.0e+10) {
				printf("\n\t No Convergence of the Newton_SDIRK Raphson Method. Now exiting to system.\n\n");
 				exit(1);
		}
		if (Newton_SDIRK_verb == 1){
			printf("\t Newton_SDIRK: iter = %3d \t Number of bicgstab-steps = %3d \t |F| = %e (%e) \n",
					iter, output_bicgstab, norm, Newton_SDIRK_tol);
		}	
	}
	if(norm > Newton_SDIRK_tol)
		printf(
		"\t Newton_SDIRK Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element.\n"
		);

	free_dvector(F,      0, Ntotal);
	free_dvector(Dkj,    0, Ntotal);

	return iter;
}
//----------------------------------------------------------
void get_BandMatrix(parameters par, double x0, double *Mj, double *kj, JFD_Components *JFD){

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
    
    Get_JFD_Matrix(par, x0, Mj, kj, JFD);

    
    
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
}
//-----------------------------------------------------------
void Get_JFD_Components(parameters par, JFD_Components *JFD){
// Calculates the remaining components of JFD
 // Note: For the band-matrix operations we use the 'Numerical Recipes in C'
 // which provide the corresponding routines for matrices and vectors with indices
 // starting from 1.
	int n=par.nslice, N = n-1,
	    m1=(*JFD).m1, m2=(*JFD).m2,
		j, k, row_K, col_K, col_J;
	double d=-1.;


	for(j=0; j<n; j++){
		row_K = j+1; // see Note
		for(k=0; k < (*JFD).ncols_J_UL[j]; k++){
			col_J = (*JFD).cols_J_UL[j][k];
			col_K = col_J+1; // see Note
			(*JFD).Ku_UL[row_K][m1+1+col_K-row_K] = (*JFD).K_UL[row_K][m1+1+col_K-row_K] = (*JFD).J_UL[j][k];	// writing K
		}
	}

	

	// Calculating the LU-decomposition of K; introducing Kl and iK,
	// see 'Numerical Recipes in C' pages 51-54
	// bandec((*JFD).K_UL, n, m1, m2, (*JFD).Kl_UL, (*JFD).iK_UL, &d);
	bandec((*JFD).Ku_UL, n, m1, m2, (*JFD).Kl_UL, (*JFD).iK_UL, &d);

	
	if(nFields1D>0){
	  
	  double *col_J_UR=dvector(0,N), *b_J_UR=dvector(0,N);
	  
	  int iF, jF;
	  for(iF=0; iF<nFields1D; iF++){
	    for(j=0; j<n; j++) b_J_UR[j] = col_J_UR[j] = (*JFD).J_UR[j][iF];
	    
	   // // banbks((*JFD).K_UL, n, (*JFD).m1, (*JFD).m2, (*JFD).Kl_UL, (*JFD).iK_UL, col_J_UR);
		banbks((*JFD).Ku_UL, n, (*JFD).m1, (*JFD).m2, (*JFD).Kl_UL, (*JFD).iK_UL, col_J_UR);
	//	mprove_band((*JFD).K_UL, (*JFD).Ku_UL, n, (*JFD).m1, (*JFD).m2, (*JFD).Kl_UL, (*JFD).iK_UL, b_J_UR, col_J_UR);
	    
	    for(j=0; j<n; j++) (*JFD).K_UR[j+1][iF+1] = col_J_UR[j];
	  }
	  
	  for(iF=0; iF<nFields1D; iF++){
	    for(jF=0; jF<nFields1D; jF++){
	    	(*JFD).Klu_LR[iF+1][jF+1] = (*JFD).K_LR[iF+1][jF+1] = (*JFD).J_LR[iF][jF];
	      
	      for(j=0; j<n; j++){
	      	(*JFD).K_LR[iF+1][jF+1]   -= (*JFD).J_LL[iF][j] * (*JFD).K_UR[j+1][jF+1];	
	      	(*JFD).Klu_LR[iF+1][jF+1] -= (*JFD).J_LL[iF][j] * (*JFD).K_UR[j+1][jF+1];	
	      }      
	    }
	  }
	  double d_LR;

	  // //ludcmp((*JFD).K_LR, nFields1D, (*JFD).indx_LR, &d_LR, 1);
	  ludcmp((*JFD).Klu_LR, nFields1D, (*JFD).indx_LR, &d_LR, 1);		  
	  
	  
	  free_dvector(col_J_UR,0,N);
	  free_dvector(b_J_UR,0,N);	  
	}

	return;
}
//----------------------------------------------------------
void free_bandMatrix(parameters par, JFD_Components *JFD){

  int n=par.nslice, N = n-1, NFields1D = nFields1D-1,       
        maxcol = nFields*STENCILSIZE;
//       
	  free_dmatrix((*JFD).J_UL, 0, N, 0, maxcol-1);
	  if(NFields1D>0){
	    free_dmatrix((*JFD).J_UR, 0, N, 0, NFields1D);
	    free_dmatrix((*JFD).J_LL, 0, NFields1D, 0, N);
	    free_dmatrix((*JFD).J_LR, 0, NFields1D, 0, NFields1D);
	    
	    free_dmatrix((*JFD).K_UR, 1, N+1,    1, nFields1D);
	    free_dmatrix((*JFD).K_LR, 1, nFields1D,    1, nFields1D);
	    free_dmatrix((*JFD).Klu_LR, 1, nFields1D,    1, nFields1D);
	    free_ivector((*JFD).indx_LR, 1, nFields1D);
	  }
	  
	  
	  free_imatrix((*JFD).cols_J_UL, 0, N, 0, maxcol-1);
	  free_ivector((*JFD).ncols_J_UL , 0, N);	  
	  free_dmatrix((*JFD).K_UL , 1, N+1,    1, (*JFD).m1+(*JFD).m2+1);
	  free_dmatrix((*JFD).Ku_UL , 1, N+1,    1, (*JFD).m1+(*JFD).m2+1);
	  free_dmatrix((*JFD).Kl_UL, 1, N+1,    1, (*JFD).m1);
	  free_ivector((*JFD).iK_UL, 1, N+1);  
}
//----------------------------------------------------------------------------------------------------------
void Get_JFD_Matrix(parameters par, double x0, double *Mj, double *kj, JFD_Components *JFD){ // Calculating JFD
	int N1=par.N1, n1=par.n1, 
	    N2=par.N2, n2=par.n2,
	    ntotal=par.nslice_total, Ntotal = ntotal-1, 
	    n=par.nslice, N = n-1, NFields1D = nFields1D-1,	    
	    I, i, j, row, column, mcol, ii, jj, J, i0, i1, j0, j1;
	
	
	double *JDX, *DX, *Null, *Eq=dvector(0, nFields-1), *EqPar=dvector(0,nFields1D);

	derivs Wslice, W1D, DWslice, DW1D, W_atGrid, DW_atGrid;
	allocate_derivs(&Wslice, 0, N);
	if(nFields1D>0) allocate_derivs(&W1D,    0, NFields1D);

	get_Ws_from_M_and_k(par, Mj, kj, Wslice, W1D);
	Derivatives_Slice(par, Wslice);
	
	allocate_derivs(&DWslice, 0, N);
	if(nFields1D>0) allocate_derivs(&DW1D,    0, NFields1D);
 	fill0_derivs(DWslice, 0, N);	
	fill0_derivs(DW1D,0, NFields1D);
	
	allocate_derivs(&W_atGrid,  0, nFields-1);
	allocate_derivs(&DW_atGrid, 0, nFields-1);

	DX=dvector(0, Ntotal);
	fill0_dvector(DX,  0, Ntotal);

	Null = dvector(0, Ntotal);
	fill0_dvector(Null, 0, Ntotal);

	JDX = dvector(0, Ntotal);
	fill0_dvector(JDX,  0, Ntotal);

	(*JFD).m1 = (*JFD).m2 = 0;

		for(column = 0; column < n; column ++){
			get_indices_from_Index_Slice(par, column, &I, &i, &j);
			double x1=par.grid_points_x1[i],
				   x2=par.grid_points_x2[j];

			DX[column]=1.;
			get_Ws_from_M_and_k(par, Null, DX, DWslice, DW1D);

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
					JDX[row]=Eq[J];

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
					JDX[row]=Eq[J];

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
				(*JFD).J_LL[J][column]=EqPar[J];
			}
			fill0_dvector(DX,  0, Ntotal);
		}

		for(column = 0; column < nFields1D; column ++){

			int column_1D = Index_Slice_Fields1D(column, n);
			DX[column_1D]=1.;
			get_Ws_from_M_and_k(par, Null, DX, DWslice, DW1D);
			fill0_derivs(DWslice, 0, N);

		

		for(jj=0; jj<=N2; jj++){
			double x2=par.grid_points_x2[jj];
			for(ii=0; ii<=N1; ii++){
				double x1=par.grid_points_x1[ii];

				get_W2D_at_Grid(par, ii, jj, Wslice, W_atGrid);
				get_W2D_at_Grid(par, ii, jj, DWslice, DW_atGrid);

				JDX_of_XY(par, x0, x1, x2, W_atGrid, W1D, DW_atGrid, DW1D, Eq);

				for(J=0; J<nFields; J++){
					int I = Index_Slice(par, J,ii,jj);
					(*JFD).J_UR[I][column] = Eq[J];
				}
			}
		}

		JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, EqPar);
		for(J=0; J<nFields1D; J++){
			(*JFD).J_LR[J][column] = EqPar[J];
		}
		fill0_dvector(DX,  0, Ntotal);
	}
	
	
	
	free_derivs(&Wslice,  0, N);
	free_derivs(&DWslice, 0, N);

	if(nFields1D>0) free_derivs(&W1D,  0, NFields1D);
	if(nFields1D>0)free_derivs(&DW1D, 0, NFields1D);

	free_derivs(&W_atGrid,  0, nFields-1);
	free_derivs(&DW_atGrid, 0, nFields-1);

	free_dvector(DX, 0, Ntotal);
	free_dvector(Null, 0, Ntotal);
	free_dvector(JDX, 0, Ntotal);
	free_dvector(Eq, 0, nFields-1);
	free_dvector(EqPar, 0, nFields1D);

	return;
}
// -------------------------------------------------------------------------------
void PreCond_FinDiff(parameters par, JFD_Components JFD, double *b, double *Y, double *res){
	int ntotal = par.nslice_total, Ntotal=ntotal-1, n=par.nslice, N = n-1, jF, j;
	double *c_U, *c_L, *b_L;
	
	c_U=dvector(0, N);
	copy_dvector(c_U, b, 0, N);
	
	// banbks(JFD.K_UL, n, JFD.m1, JFD.m2, JFD.Kl_UL, JFD.iK_UL, c_U);
	banbks(JFD.Ku_UL, n, JFD.m1, JFD.m2, JFD.Kl_UL, JFD.iK_UL, c_U);
	//mprove_band(JFD.K_UL, JFD.Ku_UL, n, JFD.m1, JFD.m2, JFD.Kl_UL, JFD.iK_UL, b, c_U);
	// mprove_band(JFD.K_UL, JFD.Ku_UL, n, JFD.m1, JFD.m2, JFD.Kl_UL, JFD.iK_UL, b, c_U);
	
	for(j=0; j<n; j++) Y[j] = c_U[j];
	
	
	
	if(nFields1D>0){
		c_L=dvector(1, nFields1D);
		b_L=dvector(1, nFields1D);

		for(jF=0; jF<nFields1D; jF++){
			c_L[jF+1] = b_L[jF+1] = b[Ntotal-jF];

			for(j=0; j<n; j++){
				c_L[jF+1]-= JFD.J_LL[jF][j]*c_U[j];
				b_L[jF+1]-= JFD.J_LL[jF][j]*c_U[j];
			}
		}

		// lubksb(JFD.K_LR, nFields1D, JFD.indx_LR, c_L, 1);
		lubksb(JFD.Klu_LR, nFields1D, JFD.indx_LR, c_L, 1);
		//mprove(JFD.K_LR, JFD.Klu_LR, nFields1D, JFD.indx_LR, b_L, c_L);

		for(jF=0; jF<nFields1D; jF++) Y[Ntotal-jF] = c_L[jF+1];

			for(j=0; j<n; j++)
				for(jF=0; jF<nFields1D; jF++)
					Y[j] -= JFD.K_UR[j+1][jF+1]*c_L[jF+1];  
	
	  	
	  	free_dvector(c_L, 1, nFields1D);
	  	free_dvector(b_L, 1, nFields1D);
	}
	free_dvector(c_U, 0, N);
	return;
}
// -------------------------------------------------------------------------------
int bicgstab_findiff(parameters par, double x0, double *Mj, double *kj, double *Dkj, double *F, double *normres){
	int ntotal=par.nslice_total, Ntotal = ntotal-1, iter, n;
	double alpha = 0., beta = 0., rho = 0., rho1 = 1., rhotol = 1e-50, omega = 1.0e-50, omegatol = 1e-50,
		*p, *ph, *q, *r, *rt, *s, *sh, *t;

	r  = dvector(0, Ntotal), rt = dvector(0, Ntotal), p  = dvector(0, Ntotal);
	// 	compute initial residual rt = r = p = F - J*DX
	J_times_Dkj_IG(par, x0, Mj, kj, Dkj, r);
	for (n = 0; n < ntotal; n++) rt[n] = r[n] = p[n] = F[n] - r[n];



	*normres = norm2(r, 0, Ntotal);
	if (*normres <= par.IG_bicgstab_findiff_tol){
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
	get_BandMatrix(par, x0, Mj, kj, &JFD);
	
	double resprecond;

	for (iter = 0; iter < IG_bicgstab_findiff_itmax; iter++) {
		rho = scalarproduct(rt, r, 0, Ntotal);
		if (fabs(rho) < rhotol)	  break;
	
		
		if (iter > 0){ // compute direction vector p
			beta = (rho/rho1)*(alpha/omega);
			for (n = 0; n < ntotal; n++) 
				p[n] = r[n] + beta * (p[n] - omega * q[n]);
		}
		// compute direction adjusting vector ph and scalar alpha :
		
		PreCond_FinDiff(par, JFD, p, ph, &resprecond); //J_FD*ph=p (solve for ph)
		J_times_Dkj_IG(par, x0, Mj, kj, ph, q); // q = J*ph
		alpha = rho/scalarproduct(rt, q, 0, Ntotal);
		
		for (n = 0; n < ntotal; n++) s[n] = r[n] - alpha * q[n];
		
		// early check of tolerance:
		*normres = norm2(s, 0, Ntotal);
	
		
		if (*normres <= par.IG_bicgstab_findiff_tol){
			for (n = 0; n < ntotal; n++)
				Dkj[n] += alpha * ph[n];
			break;
		}
		// compute stabilizer vector sh and scalar omega:

		PreCond_FinDiff(par, JFD, s, sh, &resprecond);//J_FD*sh=s (solve for sh)
		J_times_Dkj_IG(par, x0, Mj, kj, sh, t); // t = J*sh

		omega = scalarproduct(t, s, 0, Ntotal) / scalarproduct(t, t, 0, Ntotal);

		
		// compute new solution approximation:
		for (n = 0; n < ntotal; n++) {
			if(n<ntotal) Dkj[n] += alpha * ph[n] + omega * sh[n];
			r[n]   = s[n] - omega * t[n];
		}
		
		rho1 = rho;
		// are we done? 
		*normres = norm2(r, 0, Ntotal);
		if (IG_bicgstab_findiff_verb == 1) {printf("\t\t BiCGStab: iter = %2d  norm = %6.1e\r", iter+1, *normres);fflush(0); }
		if (*normres <= par.IG_bicgstab_findiff_tol)  break;
		
		if (fabs(omega) < omegatol)  break;
	}
	

	free_bandMatrix(par,&JFD);
	
	free_dvector(p,  0, ntotal-1); 		free_dvector(ph, 0, ntotal-1); 
	free_dvector(q,  0, ntotal-1); 		free_dvector(r,  0, ntotal-1);
	free_dvector(s,  0, ntotal-1); 		free_dvector(sh, 0, ntotal-1); 
	free_dvector(rt, 0, ntotal-1); 		free_dvector(t,  0, ntotal-1);
		
	/* iteration failed */
	if (iter > bicgstab_itmax) return  -1;
	
	/* breakdown */
	if (fabs(rho) < rhotol) return -10;
	if (fabs(omega) < omegatol) return -11;
	
	/* success! */
	if (IG_bicgstab_findiff_verb == 1) {printf("\t\t BiCGStab_FinDiff: iter = %2d  norm = %6.1e (%6.1e)\n", iter+1, *normres, par.IG_bicgstab_findiff_tol);fflush(0); }
	return ++iter;
}
