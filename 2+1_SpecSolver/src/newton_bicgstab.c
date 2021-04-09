#include "2+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
int newton(parameters par, double *X)
{	// Newton Raphson Method, see pages 1, 2
	int ntotal = par.ntotal, Ntotal=ntotal-1, iter=0, j, bicgstab_iter, FLAG=0/*-1*/;
	double *F, *DX, norm, normres;
	
	F     = dvector(0, Ntotal);
	DX    = dvector(0, Ntotal);
	

		
 	F_of_X(par, X, F);
  	// PrintVector(F, 0, Ntotal);

	norm = norm2(F, 0, Ntotal);
 	copy_dvector(DX, F, 0, Ntotal); 

 	double tol=Newton_tol;

 	if (Newton_verb == 1){
 		printf(" Newton Raphson Method: : SDIRK-BiCGStab\n");
 		printf(" Initial Residual: \t |F| = %e (%e)\n", norm, tol);
 		printf(" ------------------------------------------------------------------\n");
		// exit(1);
 	}
  	while(iter < Newton_itmin || (norm > tol  && iter < Newton_itmax)){
   	// while(iter < Newton_itmin || (FLAG < 0  && iter < Newton_itmax)){
 		if(norm < tol) FLAG++;
 		iter += 1;
	
 		fill0_dvector(DX, 0, Ntotal);

 		par.bicgstab_tol = bicgstab_decr/* *norm */;
	
 		
 		bicgstab_iter    = bicgstab(par, X, F, DX, &normres);	


 		for (j = 0; j < ntotal; j++) X[j] -= DX[j];

 		F_of_X(par, X, F);
 		norm = norm2(F, 0, Ntotal);
	
 		if (isinf(norm) || isnan(norm) || norm > 1.0e+10){
 			printf("\n No Convergence of the Newton Raphson Method. Now exiting to system.\n\n");
 			exit(1);
 		}
	
 		if (Newton_verb == 1){
 			printf(" Newton: iter = %1d \t Number of bicgstab-steps = %2d \t |F| = %6.1e (%6.1e) %d \n", 
 						 iter, bicgstab_iter, norm, tol, FLAG);			
 		}
 	}

 	if(norm > tol)
 		printf(" Newton Raphson Method failed to converge to prescribed tolerance.\n");


 	free_dvector(F,  0, Ntotal);
 	free_dvector(DX, 0, Ntotal);

 	return iter;
}
// -------------------------------------------------------------------------------
int bicgstab(parameters par, double *X, double *F, double *DX, double *normres)
{	// Biconjugate Gradient Stabilized Method
	int n, iter, ntotal = par.ntotal, Ntotal=ntotal-1, nslice_total=par.nslice_total, Nslice_total = nslice_total-1, N0=par.N0;
	double tol = par.bicgstab_tol, *r, *rt, *p, *ph, *q, *s, *sh, *t, alpha = 0., beta = 0., rho = 0., rho1 = 1., rhotol = 1e-50, omega = 0., omegatol = 1e-50;
	

	r  = dvector(0, Ntotal), rt = dvector(0, Ntotal), p  = dvector(0, Ntotal);
	// 	compute initial residual rt = r = p = F - J*DX        
	J_times_DX(par, X, DX, r);
	for (n = 0; n < ntotal; n++) rt[n] = r[n] = p[n] = F[n] - r[n];
	
	*normres = norm2(r, 0, Ntotal);
	if (*normres <= tol){
		free_dvector(r,  0, Ntotal);
		free_dvector(rt,  0, Ntotal);
		free_dvector(p,  0, Ntotal); 
		return 0;
	}

	
	derivs V;
	cheb_derivs chebV;

	allocate_derivs(&V, 0, Ntotal);
	allocate_cheb_derivs(&chebV, 0, Nslice_total, 0, N0);

	Get_V_From_X(par, X, V);
	get_chebV_x0(par, V, chebV);
	
	ph = dvector(0, Ntotal);
	q  = dvector(0, Ntotal);    
	s  = dvector(0, Ntotal);    
	sh = dvector(0, Ntotal);
	t  = dvector(0, Ntotal);

	for (iter = 0; iter < bicgstab_itmax; iter++) {
		rho = scalarproduct(rt, r, 0, Ntotal);
		if (fabs(rho) < rhotol) break;
		if (iter > 0){ // compute direction vector p
			beta = (rho/rho1)*(alpha/omega);
			for (n = 0; n < ntotal; n++) 
				p[n] = r[n] + beta * (p[n] - omega * q[n]);
		}
		// compute direction adjusting vector ph and scalar alpha :

		PreCond(par, chebV, p, ph); // solves J*ph=p approximately for ph

		J_times_DX(par, X, ph, q);  // q = J*ph
		
		alpha = rho/scalarproduct(rt, q, 0, Ntotal);
		for (n = 0; n < ntotal; n++) s[n] = r[n] - alpha * q[n];
		// early check of tolerance:
		*normres = norm2(s, 0, Ntotal);
		if (*normres <= tol) {
			for (n = 0; n < ntotal; n++)
				DX[n] += alpha * ph[n];
			break;
		}
		// compute stabilizer vector sh and scalar omega:
		
		PreCond(par, chebV, s, sh);  // solves J*sh=s approximately for sh
	
		J_times_DX(par, X, sh, t);    // t = J*sh
		omega = scalarproduct(t, s, 0, Ntotal) / scalarproduct(t, t, 0, Ntotal);
		// compute new solution approximation:
		for (n = 0; n < ntotal; n++) {
			DX[n] += alpha * ph[n] + omega * sh[n];
			r[n]   = s[n] - omega * t[n];
		}
		// are we done? 
		*normres = norm2(r, 0, Ntotal);
		if (bicgstab_verb == 1) {printf("\t SDIRK BiCGStab: iter = %2d  norm = %6.1e\r", iter+1, *normres);fflush(0); }
		if (*normres <= tol) break;
		rho1 = rho;
		if (fabs(omega) < omegatol) break;
	}

	free_dvector(p,      0, Ntotal);    free_dvector(ph, 0, Ntotal); 
	free_dvector(q,      0, Ntotal);    free_dvector(r,  0, Ntotal);
	free_dvector(s,      0, Ntotal);    free_dvector(sh, 0, Ntotal); 
	free_dvector(rt,     0, Ntotal);    free_dvector(t,  0, Ntotal);
	free_cheb_derivs(&chebV, 0, Nslice_total, 0, N0);
	free_derivs(&V,  0, Ntotal); 
	

	/* iteration failed */
	if (iter > bicgstab_itmax){ return -1;}
	
	/* breakdown */
	if (fabs(rho)   < rhotol)   return -10;
	if (fabs(omega) < omegatol) return -11;
	
	/* success! */
	if (bicgstab_verb == 1) {printf("\t SDIRK BiCGStab: iter = %2d  norm = %6.1e\n", iter+1, *normres);fflush(0); }
	return iter+1;
}
//-----------------------------------------------------------------------------

void get_chebV_x0(parameters par, derivs V, cheb_derivs chebV)
{
	derivs p, cp;
	int i0, N0=par.N0 , i1, N1=par.N1, i2, N2=par.N2, iField;
		
	allocate_derivs(&p, 0, N0);
	allocate_derivs(&cp, 0, N0);

	
	for(iField=0; iField<nFields; iField++){
		for(i1=0; i1<=N1; i1++){
			for(i2=0;i2<=N2;i2++){
				for(i0=0; i0<=N0; i0++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					p.X[i0] = V.X[Idx];
					p.X1[i0] = V.X1[Idx];
					p.X11[i0] = V.X11[Idx];
					p.X2[i0] = V.X2[Idx];
					p.X12[i0] = V.X12[Idx];
					p.X22[i0] = V.X22[Idx];

					p.Y[i0] = V.X[Idx];
					p.Y1[i0] = V.X1[Idx];
					p.Y11[i0] = V.X11[Idx];
					p.Y2[i0] = V.X2[Idx];
					p.Y12[i0] = V.X12[Idx];
					p.Y22[i0] = V.X22[Idx];
				}
				int Idx_slice=Index_Slice(par, iField, i1, i2);
				get_Chebyshev_Coefficients_derivs(p, cp, N0, par.grid_x0);


				for(i0=0; i0<=N0; i0++){
					chebV.X[Idx_slice][i0]   = cp.X[i0];
					chebV.X1[Idx_slice][i0]  = cp.X1[i0];
					chebV.X11[Idx_slice][i0] = cp.X11[i0];
					chebV.X2[Idx_slice][i0]  = cp.X2[i0];
					chebV.X12[Idx_slice][i0] = cp.X12[i0];
					chebV.X22[Idx_slice][i0] = cp.X22[i0];

					chebV.Y[Idx_slice][i0]	 = cp.Y[i0];
					chebV.Y1[Idx_slice][i0]	 = cp.Y1[i0];
					chebV.Y11[Idx_slice][i0] = cp.Y11[i0];	
					chebV.Y2[Idx_slice][i0]	 = cp.Y2[i0];
					chebV.Y12[Idx_slice][i0] = cp.Y12[i0];	
					chebV.Y22[Idx_slice][i0] = cp.Y22[i0];	
				}

			}
		}
	}

	for(iField=0; iField<nFields1D; iField++){
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_Fields1D(iField, par.ngrid, i0, N0);
			p.X[i0] = V.X[Idx];
			p.Y[i0] = V.X[Idx];
		}
		int Idx_slice=Index_Slice_Fields1D(iField, par.nslice);
		get_Chebyshev_Coefficients_derivs(p, cp, N0, par.grid_x0);
		for(i0=0; i0<=N0; i0++){
			chebV.X[Idx_slice][i0]   = cp.X[i0];
			chebV.Y[Idx_slice][i0]	 = cp.Y[i0];
		}
	}


	free_derivs(&p, 0, N0);
	free_derivs(&cp, 0, N0);
}
