#include "2+1_Free_Boundaries.h"
//-------------------------------------------------
void get_W3D_at_Grid(parameters par, int i0, int i1, int i2, derivs ID, derivs V, derivs W){

	int iField, N0=par.N0, N1=par.N1, N2=par.N2;
	double x0=par.grid_points_x0[i0], x0p1=1.+x0; 

	for(iField=0; iField < nFields; iField++){
		int Idx = Index(iField, i0, N0, i1, N1, i2, N2),
			Idx_slice = Index_Slice(par, iField, i1, i2);
			
			
		W.X[iField]    = ID.X[Idx_slice]   + x0p1*V.X[Idx];
		W.X1[iField]   = ID.X1[Idx_slice]  + x0p1*V.X1[Idx];
		W.X11[iField]  = ID.X11[Idx_slice] + x0p1*V.X11[Idx];
		W.X2[iField]   = ID.X2[Idx_slice]  + x0p1*V.X2[Idx];
		W.X12[iField]  = ID.X12[Idx_slice] + x0p1*V.X12[Idx];
		W.X22[iField]  = ID.X22[Idx_slice] + x0p1*V.X22[Idx];
		
		W.Y[iField]    = V.X[Idx]          + x0p1*V.Y[Idx];
		W.Y1[iField]   = V.X1[Idx]         + x0p1*V.Y1[Idx];
		W.Y11[iField]  = V.X11[Idx]        + x0p1*V.Y11[Idx];
		W.Y2[iField]   = V.X2[Idx]         + x0p1*V.Y2[Idx];
		W.Y12[iField]  = V.X12[Idx]        + x0p1*V.Y12[Idx];
		W.Y22[iField]  = V.X22[Idx]        + x0p1*V.Y22[Idx];
	}

	return;
}
//-------------------------------------------------
void get_W3D_at_slice(parameters par, int i0, int i1, int i2, derivs ID, derivs V, derivs Wslice){

	int iField, N0=par.N0, N1=par.N1, N2=par.N2;
	double x0=par.grid_points_x0[i0], x0p1=1.+x0; 

	for(iField=0; iField < nFields; iField++){
		int Idx = Index(iField, i0, N0, i1, N1, i2, N2),
			Idx_slice = Index_Slice(par, iField, i1, i2);
						
		Wslice.X[Idx_slice]   = par.ID.X[Idx_slice]   + x0p1*V.X[Idx];
		Wslice.X1[Idx_slice]  = par.ID.X1[Idx_slice]  + x0p1*V.X1[Idx];
		Wslice.X11[Idx_slice] = par.ID.X11[Idx_slice] + x0p1*V.X11[Idx];
		Wslice.X2[Idx_slice]  = par.ID.X2[Idx_slice]  + x0p1*V.X2[Idx];
		Wslice.X12[Idx_slice] = par.ID.X12[Idx_slice] + x0p1*V.X12[Idx];
		Wslice.X22[Idx_slice] = par.ID.X22[Idx_slice] + x0p1*V.X22[Idx];
		
		Wslice.Y[Idx_slice]   = V.X[Idx]          + x0p1*V.Y[Idx];
		Wslice.Y1[Idx_slice]  = V.X1[Idx]         + x0p1*V.Y1[Idx];
		Wslice.Y11[Idx_slice] = V.X11[Idx]        + x0p1*V.Y11[Idx];
		Wslice.Y2[Idx_slice]  = V.X2[Idx]         + x0p1*V.Y2[Idx];
		Wslice.Y12[Idx_slice] = V.X12[Idx]        + x0p1*V.Y12[Idx];
		Wslice.Y22[Idx_slice] = V.X22[Idx]        + x0p1*V.Y22[Idx];
	}

	return;
}

//-------------------------------------------------
void get_W1D_at_Grid(parameters par, int i0, derivs ID, derivs V, derivs W1D){

	int iField, N0=par.N0, nslice=par.nslice, ngrid=par.ngrid;
	double x0=par.grid_points_x0[i0], x0p1=1.+x0; 

	for(iField=0; iField < nFields1D; iField++){
			int Idx = Index_Fields1D(iField, ngrid, i0, N0),
				Idx_slice = Index_Slice_Fields1D(iField, nslice);
						
		W1D.X[iField]    = ID.X[Idx_slice]   + x0p1*V.X[Idx];
		W1D.Y[iField]    = V.X[Idx]          + x0p1*V.Y[Idx];
	}
}

// -------------------------------------------------------------------------------
void F_of_X(parameters par, double *X, double *F)
{
	double *Eq;
	int iField, i0, N0=par.N0, i1, N1=par.N1, i2, N2=par.N2, 
		nslice=par.nslice, ntotal=par.ntotal, ngrid=par.ngrid, naux = maximum2(nFields, nFields1D);
		
	derivs V, W1D, W, Wslice; 
	
	allocate_derivs(&V, 0, ntotal-1);
  	Get_V_From_X(par, X, V);
  	
  	allocate_derivs(&W, 0, nFields-1);
  	if(nFields1D!=0) allocate_derivs(&W1D, 0, nFields1D-1);
  	allocate_derivs(&Wslice, 0, nslice-1);
  	
	Eq = dvector(0, naux-1);

 	
	for(i0=0; i0<=N0; i0++){
		double x0=par.grid_points_x0[i0];
		get_W1D_at_Grid(par, i0, par.ID, V, W1D);
		
		for(i1=0; i1<=N1; i1++){
			double x1 = par.grid_points_x1[i1];
			
			for(i2=0; i2<=N2; i2++){
				double x2 = par.grid_points_x2[i2];
				get_W3D_at_Grid( par, i0, i1, i2, par.ID, V, W);
				get_W3D_at_slice(par, i0, i1, i2, par.ID, V, Wslice);
				
				
				F_of_XY(par, x0, x1, x2, W, W1D, Eq);
				for(iField=0; iField < nFields; iField++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					F[Idx] = Eq[iField];
				}
			}
		}
		
		FBound_of_XY(par, x0, Wslice, W1D, Eq);
		for(iField=0; iField < nFields1D; iField++){
			int Idx = Index_Fields1D(iField, ngrid, i0, N0);
			F[Idx] = Eq[iField];
		}
		
		//**********************************************//
		//												//
		//	INCLUDE CONDITIONS FOR DOMAIN TRANSITION?	//
		//				(R.P.M 01/06/2020)				//
		//**********************************************//

	}

	free_derivs(&V,      0, ntotal-1);
	free_dvector(Eq,   0, naux-1);
	free_derivs(&W,      0, nFields-1);
	if(nFields1D!=0) free_derivs(&W1D,    0, nFields1D-1);
	free_derivs(&Wslice, 0, nslice-1);
	
	return;
}
//---------------------------------------------------------------------
void J_times_DX(parameters par, double *X, double *DX,  double *JDX){
	double *DEq;
	int iField, i0, N0=par.N0, i1, N1=par.N1, i2, N2=par.N2, 
		nslice=par.nslice, nslice_total=par.nslice_total, ntotal=par.ntotal, ngrid=par.ngrid, naux = maximum2(nFields, nFields1D);
		
	derivs V, W1D, W, Wslice, DV, DW1D, DW, DWslice, Null; 
	
	allocate_derivs(&V, 0, ntotal-1); 	
	allocate_derivs(&DV, 0, ntotal-1);
  	Get_V_From_X(par, X, V);  	
  	Get_V_From_X(par, DX, DV);

  	
  	allocate_derivs(&W, 0, nFields-1);    
  	allocate_derivs(&DW, 0, nFields-1);

  	if(nFields1D!=0) {
  		allocate_derivs(&W1D, 0, nFields1D-1);		
  		allocate_derivs(&DW1D, 0, nFields1D-1);
  	}
  	
  	allocate_derivs(&Null, 0, nslice_total-1);
  	fill0_derivs(Null, 0, nslice_total-1);

  	allocate_derivs(&Wslice, 0, nslice-1);	
  	allocate_derivs(&DWslice, 0, nslice-1);
  	
	DEq = dvector(0, naux-1);

 	
	for(i0=0; i0<=N0; i0++) {
		double x0=par.grid_points_x0[i0];
		get_W1D_at_Grid(par, i0, par.ID, V, W1D);
		get_W1D_at_Grid(par, i0, Null, DV, DW1D);

		
		for(i1=0; i1<=N1; i1++){
			double x1 = par.grid_points_x1[i1];
			
			for(i2=0; i2<=N2; i2++){
				double x2 = par.grid_points_x2[i2];

				get_W3D_at_Grid(par, i0, i1, i2, par.ID, V, W);
				get_W3D_at_Grid(par, i0, i1, i2, Null, DV, DW);

				get_W3D_at_slice(par, i0, i1, i2, par.ID, V, Wslice);
				get_W3D_at_slice(par, i0, i1, i2, Null, DV, DWslice);
				
				JDX_of_XY(par, x0, x1, x2, W, W1D, DW, DW1D, DEq); // Jaux contains normal derivatives at inner boundaries
				for(iField=0; iField < nFields; iField++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					JDX[Idx] = DEq[iField];
				}
			}
		}

		JDXBound_of_XY(par, x0, Wslice, W1D, DWslice, DW1D, DEq);
		for(iField=0; iField < nFields1D; iField++){
			int Idx = Index_Fields1D(iField, ngrid, i0, N0);
			JDX[Idx] = DEq[iField];
		}
		
		//**********************************************//
		//												//
		//	INCLUDE CONDITIONS FOR DOMAIN TRANSITION?	//
		//				(R.P.M 01/06/2020)				//
		//**********************************************//

	}

	free_derivs(&V,      0, ntotal-1); 
	free_derivs(&DV,      0, ntotal-1);
	
	free_dvector(DEq,   0, naux-1);
	free_derivs(&W,      0, nFields-1); 
	free_derivs(&DW,      0, nFields-1);
	if(nFields1D!=0){
		free_derivs(&W1D,    0, nFields1D-1);
		free_derivs(&DW1D,    0, nFields1D-1);
	}
	free_derivs(&Null, 0, nslice_total-1);
	free_derivs(&Wslice, 0, nslice-1); 
	free_derivs(&DWslice, 0, nslice-1);
	
	return;

}
