#include "2+1_Free_Boundaries.h"

void get_Initial_Data(parameters *par){


	if((*par).EXACT_SOL_CHECK == 1 || (*par).EXACT_SOL_SOLVE==1){
		get_ExactSolution_ID(par);
		
		return;
	}

	int nslice=(*par).nslice, N1=(*par).N1, N2=(*par).N2, i1, i2;	
	int l=(*par).ell;
	double g, g_r, g_rr, h, h_r, h_rr, g_rrr, h_rrr, alpha=(*par).alpha, kappa=(*par).kappa, kappa2=sqr(kappa);
	
	g=alpha*plgndr(l, 0, 0.);
	h=l==0? 0: l*alpha*plgndr(l-1, 0, 0.);

	//FREE DATA
	h_r  = 10;
	h_rr = -2.;
	
	//REGULARITY CONDITION
	switch(l){		
			case 0:
					//FOR f_1
					g_r= -h_r; 
					// EXTRA REGULARITY CONDITION FOR f_3 (IF WANTED)
					g_rr = -(2*(1+kappa2)*g_r - (9*kappa2*kappa2+10*kappa2+9)*g)/12.;	
					break;
					
			case 1:
					//FOR f_2 
					g_r  = -( h_r + 2*(1+kappa2)*h );
					// EXTRA REGULARITY CONDITION FOR f_3 (IF WANTED)
					g_rr = -( h_rr/2. + (1. + 13*kappa2/5. + kappa2*kappa2 )*h  );					
					break;

			
			case 2:	
					//REGULARITY CONDITION 
					g_r  = -9./2.*(1+kappa2)*g;
					g_rr = -(0.5*h_rr + (1+kappa2)*(2827./1890.*g_r + 17./6.*h_r) );
					break;

			
			case 3:						
					//REGULARITY CONDITION
					g_r  = -(h_r + 16.*(1+kappa2)*h/3.);
					g_rr = -(1+kappa2)*(4159./315.*g_r + 169./315.*h_r);
					break;
			
			case 4:
					//REGULARITY CONDITION
					g_r  = -137./18.*(1+kappa2)*g; 
					g_rr = -(0.5*h_rr + (1+kappa2)*(173706./52745.*g_r + 34./5.*h_r));
					break;
			default:
				g_r=0.;
				g_rr=1.;
				break;
		}
					
		
	for(i1=0; i1<=N1; i1++){
		double r=func_r_of_x1(*par, (*par).grid_points_x1[i1]);
		for(i2=0; i2<=N2; i2++){
			double x=func_x_of_x2(*par, (*par).grid_points_x2[i2]);
			
			int I0 = Index_Slice(*par, 0, i1, i2),
			    I1 = Index_Slice(*par , 1, i1, i2);
			
			g_rrr = cos(2*Pi*r/(*par).r1)*exp(r/(*par).r1);
			h_rrr = 0.;	
			
			(*par).ID.X[I0] = (0.5*g_rr + r*g_rrr/6.) * plgndr(l, 0, x);
			(*par).ID.X[I1] = (0.5*h_rr + r*h_rrr/6.) * plgndr(l, 0, x);	

		
		
		}		
	}
	
	
	// ID for 0+1 FIELDS
	int I1D_0 = Index_Slice_Fields1D(0,nslice),
		I1D_1 = Index_Slice_Fields1D(1,nslice);		
		
		(*par).ID.X[I1D_0] = g_r;
		(*par).ID.X[I1D_1] = h_r;
		
	
	Derivatives_Slice(*par, (*par).ID);
	return;
}
//--------------------------------------------------------------------------------------------------------
void get_Initial_Data_newTimeStep(parameters *par, derivs W){


	if(strcmp( (*par).grid_x0,"Radau_LHS")==0 || strcmp( (*par).grid_x0,"Gauss")==0){
		fprintf(stdout, "WARNING in get_Initial_Data_newTimeStep\nTime grid = %s\nFinal Grid point does not coincide with final time\n", (*par).grid_x0);
	}
      

	int i0, N0=(*par).N0, iField, i1, i2, N1=(*par).N1, N2=(*par).N2;
   
  
	i0=0;	//Final Time step
	for(iField=0; iField < nFields1D; iField++){	  
		int indx       = Index_Fields1D(iField, (*par).ngrid, i0, N0),
			indx_slice = Index_Slice_Fields1D(iField, (*par).nslice);
			
		(*par).ID.X[indx_slice] = W.X[indx];
	}
	
	for(iField=0; iField < nFields; iField++){
		for(i1=0; i1<=N1; i1++){
			for(i2=0; i2<=N2; i2++){			
				int indx       = Index(iField, i0, N0, i1, N1, i2, N2),
					indx_slice = Index_Slice(*par, iField, i1, i2);				
					
				
				(*par).ID.X[indx_slice] = W.X[indx];				
			}
		}
	}
	
	Derivatives_Slice(*par, (*par).ID);

	return;
}
