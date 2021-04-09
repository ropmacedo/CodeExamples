#include "2+1_Free_Boundaries.h"
//-------------------------------------------------------------------------
void solve_equations(parameters par, double *X, derivs Sol){

	int newt_iter;
	
	fprintf(stdout, "\nSolving equations for t0 = %lf, t1 = %lf \n", par.t0, par.t1);

	if(par.SOLVER_METHOD == 0){
		newt_iter = newton_direct(par, X);
	}
	else if(par.SOLVER_METHOD == 1 || par.SOLVER_METHOD == 2)
		newt_iter = newton(par, X);
	else{
		fprintf(stderr, "Error in solve_equations: Method = %d not known\n", par.SOLVER_METHOD );
	}

	get_solution_from_X(par, X, Sol);

	return;
}
//-------------------------------------------------------------------------
void Get_V_From_X(parameters par, double *X, derivs V)
{
	int i, ntotal=par.ntotal;
	
	for(i=0; i<ntotal; i++) V.X[i]=X[i];
		
	Derivatives(par, V);
}
//----------------------------------------------------------------------------
void get_solution_from_X(parameters par, double *X, derivs W){
  
	int iField, i0, i1, i2, 
		N0=par.N0, N1=par.N1, N2=par.N2, 
		nslice=par.nslice, ngrid=par.ngrid, ntotal=par.ntotal;
	derivs V;
	allocate_derivs(&V, 0, ntotal-1);
  	Get_V_From_X(par, X, V);

	for(i0=0; i0<=N0; i0++){
  		double x0=par.grid_points_x0[i0], x0p1= 1.+x0;
  	
  		for(i1=0; i1<=N1; i1++){
			for(i2=0; i2<=N2; i2++){			
				for(iField=0; iField < nFields; iField++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2),
						Idx_slice = Index_Slice(par, iField, i1, i2);
						
						W.X[Idx]    = par.ID.X[Idx_slice]   + x0p1*V.X[Idx];
						W.X1[Idx]   = par.ID.X1[Idx_slice]  + x0p1*V.X1[Idx];
						W.X11[Idx]  = par.ID.X11[Idx_slice] + x0p1*V.X11[Idx];
						W.X2[Idx]   = par.ID.X2[Idx_slice]  + x0p1*V.X2[Idx];
						W.X12[Idx]  = par.ID.X12[Idx_slice] + x0p1*V.X12[Idx];
						W.X22[Idx]  = par.ID.X22[Idx_slice] + x0p1*V.X22[Idx];
						
						W.Y[Idx]    = V.X[Idx]          + x0p1*V.Y[Idx];
						W.Y1[Idx]   = V.X1[Idx]         + x0p1*V.Y1[Idx];
						W.Y11[Idx]  = V.X11[Idx]        + x0p1*V.Y11[Idx];
						W.Y2[Idx]   = V.X2[Idx]         + x0p1*V.Y2[Idx];
						W.Y12[Idx]  = V.X12[Idx]        + x0p1*V.Y12[Idx];
						W.Y22[Idx]  = V.X22[Idx]        + x0p1*V.Y22[Idx];
				
				
				
				}			
			}	
		}
		
		for(iField=0; iField < nFields1D; iField++){
			int Idx = Index_Fields1D(iField, ngrid, i0, N0),
				Idx_slice = Index_Slice_Fields1D(iField, nslice);
				
				W.X[Idx]    = par.ID.X[Idx_slice]   + x0p1*V.X[Idx];
				W.Y[Idx]    = V.X[Idx]          + x0p1*V.Y[Idx];
			
		}
	}
	
	free_derivs(&V,  0, ntotal-1);
}
