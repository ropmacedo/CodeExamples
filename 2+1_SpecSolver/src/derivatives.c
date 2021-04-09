#include "2+1_Free_Boundaries.h"


// -------------------------------------------------------------------------------
void Get_Differentiation_Matrices(parameters *par)
{
	int N0 = (*par).N0, N1 = (*par).N1, N2 = (*par).N2;
	
		(*par).D_x0 = dmatrix(0, N0, 0, N0);
		Chebyshev_Differentiation_Matrix(N0, (*par).D_x0, (*par).grid_x0);
		
		(*par).D_x1 = dmatrix(0, N1, 0, N1);
		Chebyshev_Differentiation_Matrix(N1, (*par).D_x1, (*par).grid_x1);
		
		(*par).D_x2 = dmatrix(0, N2, 0, N2);
		Chebyshev_Differentiation_Matrix(N2, (*par).D_x2, (*par).grid_x2);

}
// -------------------------------------------------------------------------------
void Free_Differentiation_Matrices(parameters *par)
{
	int N0 = (*par).N0, N1 = (*par).N1, N2 = (*par).N2;
	
	free_dmatrix((*par).D_x0, 0, N0, 0, N0);
	free_dmatrix((*par).D_x1, 0, N1, 0, N1);
	free_dmatrix((*par).D_x2, 0, N2, 0, N2);
}
// -------------------------------------------------------------------------------
void Derivatives_Slice(parameters par, derivs V)
{
	double *p, *d1p, *d2p, *q, *d1q, *d2q, *r, *d1r, *s, *d1s;
	int i1, i2, iField, N1=par.N1, N2=par.N2;
		
	  

	for(iField=0; iField < nFields; iField++){
		p   = dvector(0, N2); d1p = dvector(0, N2); d2p = dvector(0, N2);
		q   = dvector(0, N2); d1q = dvector(0, N2); d2q = dvector(0, N2);
		
		for(i1=0; i1<=N1; i1++){
			for(i2=0; i2<=N2; i2++){
				int Idx=Index_Slice(par, iField, i1, i2);
				p[i2] = V.X[Idx];
				q[i2] = V.Y[Idx];
			}
			Chebyshev_Collocations_Derivatives(   p, d1p, par.D_x2, N2);
			Chebyshev_Collocations_Derivatives( d1p, d2p, par.D_x2, N2);

			Chebyshev_Collocations_Derivatives(   q, d1q, par.D_x2, N2);
			Chebyshev_Collocations_Derivatives( d1q, d2q, par.D_x2, N2);
				
			for(i2=0; i2<=N2; i2++){
				int Idx=Index_Slice(par, iField, i1, i2);
				V.X2[Idx]  = d1p[i2];
				V.X22[Idx] = d2p[i2];

				V.Y2[Idx]  = d1q[i2];
				V.Y22[Idx] = d2q[i2];
			}
		}
		free_dvector(p,0, N2); free_dvector(d1p,0, N2); free_dvector(d2p,0, N2);
		free_dvector(q,0, N2); free_dvector(d1q,0, N2); free_dvector(d2q,0, N2);
	
		p = dvector(0, N1); d1p = dvector(0, N1); d2p = dvector(0, N1);
		q = dvector(0, N1); d1q = dvector(0, N1); d2q = dvector(0, N1);

		r = dvector(0, N1); d1r = dvector(0, N1);
		s = dvector(0, N1); d1s = dvector(0, N1);  
	
		for(i2=0; i2<=N2; i2++) {
			for(i1=0; i1<=N1; i1++){
				int Idx=Index_Slice(par, iField, i1, i2);
				p[i1] = V.X[Idx];
				r[i1] = V.X2[Idx];

				q[i1] = V.Y[Idx];
				s[i1] = V.X2[Idx]; 				
			}
			Chebyshev_Collocations_Derivatives(  p, d1p, par.D_x1, N1);
			Chebyshev_Collocations_Derivatives(d1p, d2p, par.D_x1, N1);

			Chebyshev_Collocations_Derivatives(  q, d1q, par.D_x1, N1);
			Chebyshev_Collocations_Derivatives(d1q, d2q, par.D_x1, N1);

			Chebyshev_Collocations_Derivatives(  r, d1r, par.D_x1, N1);
			Chebyshev_Collocations_Derivatives(  s, d1s, par.D_x1, N1);
				
			for(i1=0; i1<=N1; i1++){
				int Idx=Index_Slice(par, iField, i1, i2);
				V.X1[Idx]  = d1p[i1];
				V.X11[Idx] = d2p[i1];
				V.X12[Idx] = d1r[i1];

				V.Y1[Idx]  = d1q[i1];
				V.Y11[Idx] = d2q[i1];
				V.Y12[Idx] = d1s[i1];
			}
		}
		
	 free_dvector(p, 0, N1); free_dvector(d1p, 0, N1); free_dvector(d2p, 0, N1);
	 free_dvector(q, 0, N1); free_dvector(d1q, 0, N1); free_dvector(d2q, 0, N1);
	 free_dvector(  r, 0, N1); free_dvector(d1r, 0, N1);
	 free_dvector(  s, 0, N1); free_dvector(d1s, 0, N1);  		
	}
	
	  
}
// -------------------------------------------------------------------------------
void Derivatives(parameters par, derivs V)
{
	double *p, *d1p, *d2p, *q, *d1q, *d2q, *r, *d1r, *s, *d1s;
	int i0, i1, i2, iField, N0=par.N0, N1=par.N1, N2=par.N2, ngrid=par.ngrid;
		
	 

	for(iField=0; iField < nFields; iField++){
		for(i2=0; i2<=N2; i2++){
			p = dvector(0, N0); d1p = dvector(0, N0);	d2p = dvector(0, N0); 
			
			for(i1=0; i1<=N1; i1++){
				for(i0=0; i0<=N0; i0++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					p[i0]=V.X[Idx];				
				}
				Chebyshev_Collocations_Derivatives(p, d1p, par.D_x0, N0);
				for(i0=0; i0<=N0; i0++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					V.Y[Idx]=d1p[i0];				
				}				
			}
			free_dvector(p, 0, N0); free_dvector(d1p, 0, N0); free_dvector(d2p, 0, N0);
			
			p = dvector(0, N1); d1p = dvector(0, N1);	d2p = dvector(0, N1); 
			q = dvector(0, N1); d1q = dvector(0, N1);	d2q = dvector(0, N1);
			for(i0=0; i0<=N0; i0++){
				for(i1=0; i1<=N1; i1++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					p[i1]=V.X[Idx];
					q[i1]=V.Y[Idx];				
				}
				Chebyshev_Collocations_Derivatives(  p, d1p, par.D_x1, N1);
				Chebyshev_Collocations_Derivatives(d1p, d2p, par.D_x1, N1);
				
				Chebyshev_Collocations_Derivatives(  q, d1q, par.D_x1, N1);
				Chebyshev_Collocations_Derivatives(d1q, d2q, par.D_x1, N1);
				
				for(i1=0; i1<=N1; i1++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					V.X1[Idx]  = d1p[i1];
					V.X11[Idx] = d2p[i1];
					V.Y1[Idx]  = d1q[i1];				
					V.Y11[Idx] = d2q[i1];				
				}				
			}
			free_dvector(p, 0, N1); free_dvector(d1p, 0, N1); free_dvector(d2p, 0, N1);
			free_dvector(q, 0, N1); free_dvector(d1q, 0, N1); free_dvector(d2q, 0, N1);	
		}
		
		p = dvector(0, N2); d1p = dvector(0, N2);	d2p = dvector(0, N2); 
		q = dvector(0, N2); d1q = dvector(0, N2);	d2q = dvector(0, N2);
		r = dvector(0,N2); 	d1r = dvector(0,N2);
		s = dvector(0,N2); 	d1s = dvector(0,N2);
		
		for(i0=0; i0<=N0; i0++){
			for(i1=0; i1<=N1; i1++){
				for(i2=0; i2<=N2; i2++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					p[i2]=V.X[Idx];
					q[i2]=V.Y[Idx];
					r[i2]=V.X1[Idx];
					s[i2]=V.Y1[Idx];					
				}
				Chebyshev_Collocations_Derivatives(  p, d1p, par.D_x2, N2);
				Chebyshev_Collocations_Derivatives(d1p, d2p, par.D_x2, N2);
				
				Chebyshev_Collocations_Derivatives(  q, d1q, par.D_x2, N2);
				Chebyshev_Collocations_Derivatives(d1q, d2q, par.D_x2, N2);
				
				Chebyshev_Collocations_Derivatives(  r, d1r, par.D_x2, N2);
				Chebyshev_Collocations_Derivatives(  s, d1s, par.D_x2, N2);
				for(i2=0; i2<=N2; i2++){
					int Idx = Index(iField, i0, N0, i1, N1, i2, N2);
					V.X2[Idx]=d1p[i2];
					V.X22[Idx]=d2p[i2];					
					V.Y2[Idx]=d1q[i2];
					V.Y22[Idx]=d2q[i2];
					
					V.X12[Idx]=d1r[i2];
					V.Y12[Idx]=d1s[i2];		
				}
			}
		}
		free_dvector(p, 0, N2); free_dvector(d1p, 0, N2); free_dvector(d2p, 0, N2);
		free_dvector(q, 0, N2); free_dvector(d1q, 0, N2); free_dvector(d2q, 0, N2);	
		free_dvector(r, 0, N2); free_dvector(s, 0, N2);	
		free_dvector(d1r, 0, N2); free_dvector(d1s, 0, N2);	
	}
	
	
	p = dvector(0, N0); d1p = dvector(0, N0);
	for(iField=0; iField < nFields1D; iField++){
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_Fields1D(iField, ngrid, i0, N0);
			
			p[i0]=V.X[Idx];
		}
		Chebyshev_Collocations_Derivatives(  p, d1p, par.D_x0, N0);
		
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_Fields1D(iField, ngrid, i0, N0);
			
			V.Y[Idx]=d1p[i0];
		}
			
	}
	free_dvector(p, 0, N0); free_dvector(d1p, 0, N0);			
	return;
}
// -------------------------------------------------------------------------------
void Get_DerivativesFinDif(parameters par, derivs V)
{
	int i_field, i1, n1=par.n1, 
	    i2,  n2=par.n2;
	       
	for(i_field=0; i_field<nFields; i_field++)
	  for(i1=0; i1<n1; i1++)
	    for(i2=0; i2< n2; i2++)
	    	Get_DerivativesFinDif_grid( par, i_field, i1,  i2,   V);

	return;
}
// -------------------------------------------------------------------------------
void Get_DerivativesFinDif_grid(parameters par, int i_field, int i1, int i2, derivs V)
{
	int indx =  	   Index_Slice(par, i_field, i1,   i2),
	indx_x1_P1 = 	   Index_Slice(par, i_field, i1+1, i2),
	indx_x1_M1 = 	   Index_Slice(par, i_field, i1-1, i2),
	indx_x1_P2 = 	   Index_Slice(par, i_field, i1+2, i2),
	indx_x1_M2 = 	   Index_Slice(par, i_field, i1-2, i2),
	indx_x1_P3 = 	   Index_Slice(par, i_field, i1+3, i2),
	indx_x1_M3 = 	   Index_Slice(par, i_field, i1-3, i2),
	indx_x1_P4 = 	   Index_Slice(par, i_field, i1+4, i2),
	indx_x1_M4 = 	   Index_Slice(par, i_field, i1-4, i2),

	indx_x2_P1 = 	   Index_Slice(par, i_field, i1, i2+1),
	indx_x2_M1 = 	   Index_Slice(par, i_field, i1, i2-1),
	indx_x2_P2 = 	   Index_Slice(par, i_field, i1, i2+2),
	indx_x2_M2 = 	   Index_Slice(par, i_field, i1, i2-2),
	indx_x2_P3 = 	   Index_Slice(par, i_field, i1, i2+3),
	indx_x2_M3 = 	   Index_Slice(par, i_field, i1, i2-3),
	indx_x2_P4 = 	   Index_Slice(par, i_field, i1, i2+4),
	indx_x2_M4 = 	   Index_Slice(par, i_field, i1, i2-4),
	
	indx_x1_P3_x2_M3 = Index_Slice(par, i_field, i1+3, i2-3),
	indx_x1_P3_x2_M2 = Index_Slice(par, i_field, i1+3, i2-2),
	indx_x1_P3_x2_M1 = Index_Slice(par, i_field, i1+3, i2-1),
	indx_x1_P3_x2_P1 = Index_Slice(par, i_field, i1+3, i2+1),
	indx_x1_P3_x2_P2 = Index_Slice(par, i_field, i1+3, i2+2),
	indx_x1_P3_x2_P3 = Index_Slice(par, i_field, i1+3, i2+3),
	 
	indx_x1_P2_x2_M3 = Index_Slice(par, i_field, i1+2, i2-3),
	indx_x1_P2_x2_M2 = Index_Slice(par, i_field, i1+2, i2-2),
	indx_x1_P2_x2_M1 = Index_Slice(par, i_field, i1+2, i2-1),
	indx_x1_P2_x2_P1 = Index_Slice(par, i_field, i1+2, i2+1),
	indx_x1_P2_x2_P2 = Index_Slice(par, i_field, i1+2, i2+2),
	indx_x1_P2_x2_P3 = Index_Slice(par, i_field, i1+2, i2+3),
	
	indx_x1_P1_x2_M3 = Index_Slice(par, i_field, i1+1, i2-3),
	indx_x1_P1_x2_M2 = Index_Slice(par, i_field, i1+1, i2-2),
	indx_x1_P1_x2_M1 = Index_Slice(par, i_field, i1+1, i2-1),
	indx_x1_P1_x2_P1 = Index_Slice(par, i_field, i1+1, i2+1),
	indx_x1_P1_x2_P2 = Index_Slice(par, i_field, i1+1, i2+2),
	indx_x1_P1_x2_P3 = Index_Slice(par, i_field, i1+1, i2+3),

	indx_x1_M3_x2_M3 = Index_Slice(par, i_field, i1-3, i2-3),
	indx_x1_M3_x2_M2 = Index_Slice(par, i_field, i1-3, i2-2),
	indx_x1_M3_x2_M1 = Index_Slice(par, i_field, i1-3, i2-1),
	indx_x1_M3_x2_P1 = Index_Slice(par, i_field, i1-3, i2+1),
	indx_x1_M3_x2_P2 = Index_Slice(par, i_field, i1-3, i2+2),
	indx_x1_M3_x2_P3 = Index_Slice(par, i_field, i1-3, i2+3),
	
	indx_x1_M2_x2_M3 = Index_Slice(par, i_field, i1-2, i2-3),
	indx_x1_M2_x2_M2 = Index_Slice(par, i_field, i1-2, i2-2),
	indx_x1_M2_x2_M1 = Index_Slice(par, i_field, i1-2, i2-1),
	indx_x1_M2_x2_P1 = Index_Slice(par, i_field, i1-2, i2+1),
	indx_x1_M2_x2_P2 = Index_Slice(par, i_field, i1-2, i2+2),
	indx_x1_M2_x2_P3 = Index_Slice(par, i_field, i1-2, i2+3),
	
	indx_x1_M1_x2_M3 = Index_Slice(par, i_field, i1-1, i2-3),
	indx_x1_M1_x2_M2 = Index_Slice(par, i_field, i1-1, i2-2),
	indx_x1_M1_x2_M1 = Index_Slice(par, i_field, i1-1, i2-1),
	indx_x1_M1_x2_P1 = Index_Slice(par, i_field, i1-1, i2+1),
	indx_x1_M1_x2_P2 = Index_Slice(par, i_field, i1-1, i2+2),
	indx_x1_M1_x2_P3 = Index_Slice(par, i_field, i1-1, i2+3);		
	
	double X_da1, X_d2a1, X_d4a1, Y_da1, Y_d2a1, Y_d4a1, 
		   h1, h1_2, h1_4, cos1, cos1_2, cos1_3, sin1, sin1_2,sin1_3,sin1_4,
	       X_da2, X_d2a2, X_d4a2, Y_da2, Y_d2a2, Y_d4a2, 
	       h2, h2_2, h2_4, cos2, cos2_2, cos2_3, sin2, sin2_2,sin2_3,sin2_4,
	       X_da1a2, X_d2a1a2, X_d2a2a1, X_d2a1_2a2, Y_da1a2, Y_d2a1a2, Y_d2a2a1, Y_d2a1_2a2;
	
	// calculation of Derivatives w.r.t. x1-Direction --- x1 = cos(a1) ---(Finite Dif)
	h1   = par.angle_incr_x1;
	h1_2 = sqr(h1); 
	h1_4 = sqr(h1_2);
	
	cos1 = par.grid_points_x1[i1];
	cos1_2=sqr(cos1); 
	cos1_3=cos1*cos1_2; 
	sin1=sqrt(1.-cos1_2); 
	sin1_2=sqr(sin1);
	sin1_3=sin1*sin1_2;
	sin1_4=sqr(sin1_2);
	
	if( FD_ORDER ==2 ){
	    X_da1 =  (0.5*V.X[indx_x1_P1]- 0.5*V.X[indx_x1_M1])/(h1),
	    X_d2a1 = ( V.X[indx_x1_M1] - 2*V.X[indx] + V.X[indx_x1_P1] )/(h1_2);
	    X_d4a1 = (V.X[indx_x1_M2] -4*V.X[indx_x1_M1] + 6*V.X[indx] - 4*V.X[indx_x1_P1] + V.X[indx_x1_P2])/(h1_4);

	    Y_da1 =  (0.5*V.Y[indx_x1_P1]- 0.5*V.Y[indx_x1_M1])/(h1),
	    Y_d2a1 = ( V.Y[indx_x1_M1] - 2*V.Y[indx] + V.Y[indx_x1_P1] )/(h1_2);
	    Y_d4a1 = (V.Y[indx_x1_M2] -4*V.Y[indx_x1_M1] + 6*V.Y[indx] - 4*V.Y[indx_x1_P1] + V.Y[indx_x1_P2])/(h1_4);
	}		
	else if (FD_ORDER == 4){
		X_da1 =( (1./12)*V.X[indx_x1_M2] -(2./3)*V.X[indx_x1_M1] + (2./3)*V.X[indx_x1_P1] - (1./12)*V.X[indx_x1_P2])/(h1);
		X_d2a1=(-(1./12)*V.X[indx_x1_M2] +(4./3)*V.X[indx_x1_M1] -(5./2)*V.X[indx] + (4./3)*V.X[indx_x1_P1] - (1./12)*V.X[indx_x1_P2])/(h1_2);
		X_d4a1=(-(1./6)*V.X[indx_x1_M3] + 2*V.X[indx_x1_M2] - (13./2)*V.X[indx_x1_M1] +(28./3)*V.X[indx] - (13./2)*V.X[indx_x1_P1] + 2*V.X[indx_x1_P2] - (1./6)*V.X[indx_x1_P3])/(h1_4);

		Y_da1 =( (1./12)*V.Y[indx_x1_M2] -(2./3)*V.Y[indx_x1_M1] + (2./3)*V.Y[indx_x1_P1] - (1./12)*V.Y[indx_x1_P2])/(h1);
		Y_d2a1=(-(1./12)*V.Y[indx_x1_M2] +(4./3)*V.Y[indx_x1_M1] -(5./2)*V.Y[indx] + (4./3)*V.Y[indx_x1_P1] - (1./12)*V.Y[indx_x1_P2])/(h1_2);
		Y_d4a1=(-(1./6)*V.Y[indx_x1_M3] + 2*V.Y[indx_x1_M2] - (13./2)*V.Y[indx_x1_M1] +(28./3)*V.Y[indx] - (13./2)*V.Y[indx_x1_P1] + 2*V.Y[indx_x1_P2] - (1./6)*V.Y[indx_x1_P3])/(h1_4);
	}
	else if (FD_ORDER == 6){
		X_da1 =( -(1./60)*V.X[indx_x1_M3] + (3./20)*V.X[indx_x1_M2] -(3./4)*V.X[indx_x1_M1] + (3./4)*V.X[indx_x1_P1] - (3./20)*V.X[indx_x1_P2]+ (1./60)*V.X[indx_x1_P3])/(h1);
		X_d2a1=((1./90)*V.X[indx_x1_M3]-(3./20)*V.X[indx_x1_M2] + (3./2)*V.X[indx_x1_M1] -(49./18)*V.X[indx] + (3./2)*V.X[indx_x1_P1] - (3./20)*V.X[indx_x1_P2] + (1./90)*V.X[indx_x1_P3])/(h1_2);
		X_d4a1=((7./240)*V.X[indx_x1_M4]-(2./5)*V.X[indx_x1_M3] + (169./60)*V.X[indx_x1_M2] - (122./15)*V.X[indx_x1_M1] + (91./8)*V.X[indx] - (122./15)*V.X[indx_x1_P1] + (169./60)*V.X[indx_x1_P2] - (2./5)*V.X[indx_x1_P3]+(7./240)*V.X[indx_x1_P4])/(h1_4);

		Y_da1 =( -(1./60)*V.Y[indx_x1_M3] + (3./20)*V.Y[indx_x1_M2] -(3./4)*V.Y[indx_x1_M1] + (3./4)*V.Y[indx_x1_P1] - (3./20)*V.Y[indx_x1_P2]+ (1./60)*V.Y[indx_x1_P3])/(h1);
		Y_d2a1=((1./90)*V.Y[indx_x1_M3]-(3./20)*V.Y[indx_x1_M2] + (3./2)*V.Y[indx_x1_M1] -(49./18)*V.Y[indx] + (3./2)*V.Y[indx_x1_P1] - (3./20)*V.Y[indx_x1_P2] + (1./90)*V.Y[indx_x1_P3])/(h1_2);
		Y_d4a1=((7./240)*V.Y[indx_x1_M4]-(2./5)*V.Y[indx_x1_M3] + (169./60)*V.Y[indx_x1_M2] - (122./15)*V.Y[indx_x1_M1] + (91./8)*V.Y[indx] - (122./15)*V.Y[indx_x1_P1] + (169./60)*V.Y[indx_x1_P2] - (2./5)*V.Y[indx_x1_P3]+(7./240)*V.Y[indx_x1_P4])/(h1_4);
	}
	
	if(strcmp( par.grid_x1,"Fourier" ) ==0){
	    V.X1[indx]=X_da1;
	    V.X11[indx]=X_d2a1;

	    V.Y1[indx]=Y_da1;
	    V.Y11[indx]=Y_d2a1;
	}
	else{ 
	    if( fabs(fabs(cos1)-1.) < TINY){
		    V.X1[indx] = -X_d2a1/cos1;
		    V.X11[indx] = (X_d4a1 + X_d2a1)/3.;

		    V.Y1[indx] = -Y_d2a1/cos1;
		    V.Y11[indx] = (Y_d4a1 + Y_d2a1)/3.;
	    }
	    else
	    {
		    V.X1[indx] = -X_da1/sin1;
		    V.X11[indx] = (X_d2a1 + cos1*V.X1[indx])/sin1_2;

		    V.Y1[indx] = -Y_da1/sin1;
		    V.Y11[indx] = (Y_d2a1 + cos1*V.Y1[indx])/sin1_2;
	    }
	}
	
	// calculation of Derivatives w.r.t. x2-Direction --- x2 = cos(a2) ---(Finite Dif)
	h2   = par.angle_incr_x2;
	h2_2 = sqr(h2); 
	h2_4 = sqr(h2_2);
	
	cos2   = par.grid_points_x2[i2]; 
	cos2_2 = sqr(cos2); 
	cos2_3 = cos2*cos2_2; 
	sin2   = sqrt(1.-cos2_2); 
	sin2_2 = sqr(sin2); 
	sin2_3 = sin2*sin2_2; 
	sin2_4 = sqr(sin2_2);
	

	if( FD_ORDER ==2 ){
	    X_da2 =  (0.5*V.X[indx_x2_P1]- 0.5*V.X[indx_x2_M1])/(h2),
	    X_d2a2 = ( V.X[indx_x2_M1] - 2*V.X[indx] + V.X[indx_x2_P1] )/(h2_2);
	    X_d4a2 = (V.X[indx_x2_M2] -4*V.X[indx_x2_M1] + 6*V.X[indx] - 4*V.X[indx_x2_P1] + V.X[indx_x2_P2])/(h2_4);

	    Y_da2 =  (0.5*V.Y[indx_x2_P1]- 0.5*V.Y[indx_x2_M1])/(h2),
	    Y_d2a2 = ( V.Y[indx_x2_M1] - 2*V.Y[indx] + V.Y[indx_x2_P1] )/(h2_2);
	    Y_d4a2 = (V.Y[indx_x2_M2] -4*V.Y[indx_x2_M1] + 6*V.Y[indx] - 4*V.Y[indx_x2_P1] + V.Y[indx_x2_P2])/(h2_4);
	}		
	else if (FD_ORDER == 4){
		X_da2 =( (1./12)*V.X[indx_x2_M2] -(2./3)*V.X[indx_x2_M1] + (2./3)*V.X[indx_x2_P1] - (1./12)*V.X[indx_x2_P2])/(h2);
		X_d2a2=(-(1./12)*V.X[indx_x2_M2] +(4./3)*V.X[indx_x2_M1] -(5./2)*V.X[indx] + (4./3)*V.X[indx_x2_P1] - (1./12)*V.X[indx_x2_P2])/(h2_2);
		X_d4a2=(-(1./6)*V.X[indx_x2_M3] + 2*V.X[indx_x2_M2] - (13./2)*V.X[indx_x2_M1] +(28./3)*V.X[indx] - (13./2)*V.X[indx_x2_P1] + 2*V.X[indx_x2_P2] - (1./6)*V.X[indx_x2_P3])/(h2_4);

		Y_da2 =( (1./12)*V.Y[indx_x2_M2] -(2./3)*V.Y[indx_x2_M1] + (2./3)*V.Y[indx_x2_P1] - (1./12)*V.Y[indx_x2_P2])/(h2);
		Y_d2a2=(-(1./12)*V.Y[indx_x2_M2] +(4./3)*V.Y[indx_x2_M1] -(5./2)*V.Y[indx] + (4./3)*V.Y[indx_x2_P1] - (1./12)*V.Y[indx_x2_P2])/(h2_2);
		Y_d4a2=(-(1./6)*V.Y[indx_x2_M3] + 2*V.Y[indx_x2_M2] - (13./2)*V.Y[indx_x2_M1] +(28./3)*V.Y[indx] - (13./2)*V.Y[indx_x2_P1] + 2*V.Y[indx_x2_P2] - (1./6)*V.Y[indx_x2_P3])/(h2_4);
	}
	else if (FD_ORDER == 6){
		X_da2 =( -(1./60)*V.X[indx_x2_M3] + (3./20)*V.X[indx_x2_M2] -(3./4)*V.X[indx_x2_M1] + (3./4)*V.X[indx_x2_P1] - (3./20)*V.X[indx_x2_P2]+ (1./60)*V.X[indx_x2_P3])/(h2);
		X_d2a2=((1./90)*V.X[indx_x2_M3]-(3./20)*V.X[indx_x2_M2] + (3./2)*V.X[indx_x2_M1] -(49./18)*V.X[indx] + (3./2)*V.X[indx_x2_P1] - (3./20)*V.X[indx_x2_P2] + (1./90)*V.X[indx_x2_P3])/(h2_2);
		X_d4a2=((7./240)*V.X[indx_x2_M4]-(2./5)*V.X[indx_x2_M3] + (169./60)*V.X[indx_x2_M2] - (122./15)*V.X[indx_x2_M1] + (91./8)*V.X[indx] - (122./15)*V.X[indx_x2_P1] + (169./60)*V.X[indx_x2_P2] - (2./5)*V.X[indx_x2_P3]+(7./240)*V.X[indx_x2_P4])/(h2_4);

		Y_da2 =( -(1./60)*V.Y[indx_x2_M3] + (3./20)*V.Y[indx_x2_M2] -(3./4)*V.Y[indx_x2_M1] + (3./4)*V.Y[indx_x2_P1] - (3./20)*V.Y[indx_x2_P2]+ (1./60)*V.Y[indx_x2_P3])/(h2);
		Y_d2a2=((1./90)*V.Y[indx_x2_M3]-(3./20)*V.Y[indx_x2_M2] + (3./2)*V.Y[indx_x2_M1] -(49./18)*V.Y[indx] + (3./2)*V.Y[indx_x2_P1] - (3./20)*V.Y[indx_x2_P2] + (1./90)*V.Y[indx_x2_P3])/(h2_2);
		Y_d4a2=((7./240)*V.Y[indx_x2_M4]-(2./5)*V.Y[indx_x2_M3] + (169./60)*V.Y[indx_x2_M2] - (122./15)*V.Y[indx_x2_M1] + (91./8)*V.Y[indx] - (122./15)*V.Y[indx_x2_P1] + (169./60)*V.Y[indx_x2_P2] - (2./5)*V.Y[indx_x2_P3]+(7./240)*V.Y[indx_x2_P4])/(h2_4);
	}
	
	if(strcmp( par.grid_x2,"Fourier" ) ==0){
	    V.X2[indx]=X_da2;
	    V.X22[indx]=X_d2a2;

	    V.Y2[indx]=Y_da2;
	    V.Y22[indx]=Y_d2a2;
	}
	else{ 
	    if( fabs(fabs(cos2)-1.) < TINY){
		    V.X2[indx] = -X_d2a2/cos2;
		    V.X22[indx] = (X_d4a2 + X_d2a2)/3.;

		    V.Y2[indx] = -Y_d2a2/cos2;
		    V.Y22[indx] = (Y_d4a2 + Y_d2a2)/3.;
	    }
	    else
	    {
		    V.X2[indx] = -X_da2/sin2;
		    V.X22[indx] = (X_d2a2 + cos2*V.X2[indx])/sin2_2;

		    V.Y2[indx] = -Y_da2/sin2;
		    V.Y22[indx] = (Y_d2a2 + cos2*V.Y2[indx])/sin2_2;
	    }
	}
	
    // calculation of mixed Derivatives w.r.t. x1 and x2-Directions (Finite Dif)
	if( FD_ORDER ==2 ){
	  X_da1a2 = (
	          0.5*(0.5*V.X[indx_x1_P1_x2_P1]- 0.5*V.X[indx_x1_P1_x2_M1])/(h2) 
	        - 0.5*(0.5*V.X[indx_x1_M1_x2_P1]- 0.5*V.X[indx_x1_M1_x2_M1])/(h2) 
		)/(h1);
	  
	  X_d2a1a2 = ( 
		  (0.5*V.X[indx_x1_M1_x2_P1]- 0.5*V.X[indx_x1_M1_x2_M1])/(h2)
		- 2*(0.5*V.X[indx_x2_P1]- 0.5*V.X[indx_x2_M1])/(h2) 
		+ (0.5*V.X[indx_x1_P1_x2_P1]- 0.5*V.X[indx_x1_P1_x2_M1])/(h2)
		)/(h1_2);
		
	  X_d2a2a1 =  (
		    0.5*( V.X[indx_x1_P1_x2_M1] - 2*V.X[indx_x1_P1] + V.X[indx_x1_P1_x2_P1] )/(h2_2)
		  - 0.5*( V.X[indx_x1_M1_x2_M1] - 2*V.X[indx_x1_M1] + V.X[indx_x1_M1_x2_P1] )/(h2_2)
		  )/(h1);
	  
	  X_d2a1_2a2 = ( 
		     ( V.X[indx_x1_M1_x2_M1] - 2*V.X[indx_x1_M1] + V.X[indx_x1_M1_x2_P1] )/(h2_2) 
		    -2*( V.X[indx_x2_M1] - 2*V.X[indx] + V.X[indx_x2_P1] )/(h2_2)
		    +( V.X[indx_x1_P1_x2_M1] - 2*V.X[indx_x1_P1] + V.X[indx_x1_P1_x2_P1] )/(h2_2)
		  )/(h1_2);

	  Y_da1a2 = (
	          0.5*(0.5*V.Y[indx_x1_P1_x2_P1]- 0.5*V.Y[indx_x1_P1_x2_M1])/(h2) 
	        - 0.5*(0.5*V.Y[indx_x1_M1_x2_P1]- 0.5*V.Y[indx_x1_M1_x2_M1])/(h2) 
		)/(h1);
	  
	  Y_d2a1a2 = ( 
		  (0.5*V.Y[indx_x1_M1_x2_P1]- 0.5*V.Y[indx_x1_M1_x2_M1])/(h2)
		- 2*(0.5*V.Y[indx_x2_P1]- 0.5*V.Y[indx_x2_M1])/(h2) 
		+ (0.5*V.Y[indx_x1_P1_x2_P1]- 0.5*V.Y[indx_x1_P1_x2_M1])/(h2)
		)/(h1_2);
		
	  Y_d2a2a1 =  (
		    0.5*( V.Y[indx_x1_P1_x2_M1] - 2*V.Y[indx_x1_P1] + V.Y[indx_x1_P1_x2_P1] )/(h2_2)
		  - 0.5*( V.Y[indx_x1_M1_x2_M1] - 2*V.Y[indx_x1_M1] + V.Y[indx_x1_M1_x2_P1] )/(h2_2)
		  )/(h1);
	  
	  Y_d2a1_2a2 = ( 
		     ( V.Y[indx_x1_M1_x2_M1] - 2*V.Y[indx_x1_M1] + V.Y[indx_x1_M1_x2_P1] )/(h2_2) 
		    -2*( V.Y[indx_x2_M1] - 2*V.Y[indx] + V.Y[indx_x2_P1] )/(h2_2)
		    +( V.Y[indx_x1_P1_x2_M1] - 2*V.Y[indx_x1_P1] + V.Y[indx_x1_P1_x2_P1] )/(h2_2)
		  )/(h1_2);
	

	}		
	else if (FD_ORDER == 4){
		X_da1a2 =( 
		      (1./12)*( (1./12)*V.X[indx_x1_M2_x2_M2] -(2./3)*V.X[indx_x1_M2_x2_M1] + (2./3)*V.X[indx_x1_M2_x2_P1] - (1./12)*V.X[indx_x1_M2_x2_P2])/(h2)
		     -(2./3)*( (1./12)*V.X[indx_x1_M1_x2_M2] -(2./3)*V.X[indx_x1_M1_x2_M1] + (2./3)*V.X[indx_x1_M1_x2_P1] - (1./12)*V.X[indx_x1_M1_x2_P2])/(h2)
		    + (2./3)*( (1./12)*V.X[indx_x1_P1_x2_M2] -(2./3)*V.X[indx_x1_P1_x2_M1] + (2./3)*V.X[indx_x1_P1_x2_P1] - (1./12)*V.X[indx_x1_P1_x2_P2])/(h2)
		    - (1./12)*( (1./12)*V.X[indx_x1_P2_x2_M2] -(2./3)*V.X[indx_x1_P2_x2_M1] + (2./3)*V.X[indx_x1_P2_x2_P1] - (1./12)*V.X[indx_x1_P2_x2_P2])/(h2)
		    )/(h1);
		    
		X_d2a1a2=(
		      -(1./12)*( (1./12)*V.X[indx_x1_M2_x2_M2] -(2./3)*V.X[indx_x1_M2_x2_M1] + (2./3)*V.X[indx_x1_M2_x2_P1] - (1./12)*V.X[indx_x1_M2_x2_P2])/(h2)
		      +(4./3)*( (1./12)*V.X[indx_x1_M1_x2_M2] -(2./3)*V.X[indx_x1_M1_x2_M1] + (2./3)*V.X[indx_x1_M1_x2_P1] - (1./12)*V.X[indx_x1_M1_x2_P2])/(h2)
		      -(5./2)*( (1./12)*V.X[indx_x2_M2] -(2./3)*V.X[indx_x2_M1] + (2./3)*V.X[indx_x2_P1] - (1./12)*V.X[indx_x2_P2])/(h2) 
		      +(4./3)*( (1./12)*V.X[indx_x1_P1_x2_M2] -(2./3)*V.X[indx_x1_P1_x2_M1] + (2./3)*V.X[indx_x1_P1_x2_P1] - (1./12)*V.X[indx_x1_P1_x2_P2])/(h2)
		      - (1./12)*( (1./12)*V.X[indx_x1_P2_x2_M2] -(2./3)*V.X[indx_x1_P2_x2_M1] + (2./3)*V.X[indx_x1_P2_x2_P1] - (1./12)*V.X[indx_x1_P2_x2_P2])/(h2)
		    )/(h1_2);
		X_d2a2a1 =( 
			(1./12)*(-(1./12)*V.X[indx_x1_M2_x2_M2] +(4./3)*V.X[indx_x1_M2_x2_M1] -(5./2)*V.X[indx_x1_M2] + (4./3)*V.X[indx_x1_M2_x2_P1] - (1./12)*V.X[indx_x1_M2_x2_P2])/(h2_2)
			-(2./3)*(-(1./12)*V.X[indx_x1_M1_x2_M2] +(4./3)*V.X[indx_x1_M1_x2_M1] -(5./2)*V.X[indx_x1_M1] + (4./3)*V.X[indx_x1_M1_x2_P1] - (1./12)*V.X[indx_x1_M1_x2_P2])/(h2_2)
			+(2./3)*(-(1./12)*V.X[indx_x1_P1_x2_M2] +(4./3)*V.X[indx_x1_P1_x2_M1] -(5./2)*V.X[indx_x1_P1] + (4./3)*V.X[indx_x1_P1_x2_P1] - (1./12)*V.X[indx_x1_P1_x2_P2])/(h2_2)
			-(1./12)*(-(1./12)*V.X[indx_x1_P2_x2_M2] +(4./3)*V.X[indx_x1_P2_x2_M1] -(5./2)*V.X[indx_x1_P2] + (4./3)*V.X[indx_x1_P2_x2_P1] - (1./12)*V.X[indx_x1_P2_x2_P2])/(h2_2)
		      )/(h1);
		      
		X_d2a1_2a2 = ( 
			(-(1./12)*V.X[indx_x1_M1_x2_M2] +(4./3)*V.X[indx_x1_M1_x2_M1] -(5./2)*V.X[indx_x1_M1] + (4./3)*V.X[indx_x1_M1_x2_P1] - (1./12)*V.X[indx_x1_M1_x2_P2])/(h2_2) 
		     - 2*(-(1./12)*V.X[indx_x2_M2] +(4./3)*V.X[indx_x2_M1] -(5./2)*V.X[indx] + (4./3)*V.X[indx_x2_P1] - (1./12)*V.X[indx_x2_P2])/(h2_2) 
		     + (-(1./12)*V.X[indx_x1_P1_x2_M2] +(4./3)*V.X[indx_x1_P1_x2_M1] -(5./2)*V.X[indx_x1_P1] + (4./3)*V.X[indx_x1_P1_x2_P1] - (1./12)*V.X[indx_x1_P1_x2_P2])/(h2_2) 
		      )/(h1_2);

		Y_da1a2 =( 
		      (1./12)*( (1./12)*V.Y[indx_x1_M2_x2_M2] -(2./3)*V.Y[indx_x1_M2_x2_M1] + (2./3)*V.Y[indx_x1_M2_x2_P1] - (1./12)*V.Y[indx_x1_M2_x2_P2])/(h2)
		     -(2./3)*( (1./12)*V.Y[indx_x1_M1_x2_M2] -(2./3)*V.Y[indx_x1_M1_x2_M1] + (2./3)*V.Y[indx_x1_M1_x2_P1] - (1./12)*V.Y[indx_x1_M1_x2_P2])/(h2)
		    + (2./3)*( (1./12)*V.Y[indx_x1_P1_x2_M2] -(2./3)*V.Y[indx_x1_P1_x2_M1] + (2./3)*V.Y[indx_x1_P1_x2_P1] - (1./12)*V.Y[indx_x1_P1_x2_P2])/(h2)
		    - (1./12)*( (1./12)*V.Y[indx_x1_P2_x2_M2] -(2./3)*V.Y[indx_x1_P2_x2_M1] + (2./3)*V.Y[indx_x1_P2_x2_P1] - (1./12)*V.Y[indx_x1_P2_x2_P2])/(h2)
		    )/(h1);
		    
		Y_d2a1a2=(
		      -(1./12)*( (1./12)*V.Y[indx_x1_M2_x2_M2] -(2./3)*V.Y[indx_x1_M2_x2_M1] + (2./3)*V.Y[indx_x1_M2_x2_P1] - (1./12)*V.Y[indx_x1_M2_x2_P2])/(h2)
		      +(4./3)*( (1./12)*V.Y[indx_x1_M1_x2_M2] -(2./3)*V.Y[indx_x1_M1_x2_M1] + (2./3)*V.Y[indx_x1_M1_x2_P1] - (1./12)*V.Y[indx_x1_M1_x2_P2])/(h2)
		      -(5./2)*( (1./12)*V.Y[indx_x2_M2] -(2./3)*V.Y[indx_x2_M1] + (2./3)*V.Y[indx_x2_P1] - (1./12)*V.Y[indx_x2_P2])/(h2) 
		      +(4./3)*( (1./12)*V.Y[indx_x1_P1_x2_M2] -(2./3)*V.Y[indx_x1_P1_x2_M1] + (2./3)*V.Y[indx_x1_P1_x2_P1] - (1./12)*V.Y[indx_x1_P1_x2_P2])/(h2)
		      - (1./12)*( (1./12)*V.Y[indx_x1_P2_x2_M2] -(2./3)*V.Y[indx_x1_P2_x2_M1] + (2./3)*V.Y[indx_x1_P2_x2_P1] - (1./12)*V.Y[indx_x1_P2_x2_P2])/(h2)
		    )/(h1_2);
		Y_d2a2a1 =( 
			(1./12)*(-(1./12)*V.Y[indx_x1_M2_x2_M2] +(4./3)*V.Y[indx_x1_M2_x2_M1] -(5./2)*V.Y[indx_x1_M2] + (4./3)*V.Y[indx_x1_M2_x2_P1] - (1./12)*V.Y[indx_x1_M2_x2_P2])/(h2_2)
			-(2./3)*(-(1./12)*V.Y[indx_x1_M1_x2_M2] +(4./3)*V.Y[indx_x1_M1_x2_M1] -(5./2)*V.Y[indx_x1_M1] + (4./3)*V.Y[indx_x1_M1_x2_P1] - (1./12)*V.Y[indx_x1_M1_x2_P2])/(h2_2)
			+(2./3)*(-(1./12)*V.Y[indx_x1_P1_x2_M2] +(4./3)*V.Y[indx_x1_P1_x2_M1] -(5./2)*V.Y[indx_x1_P1] + (4./3)*V.Y[indx_x1_P1_x2_P1] - (1./12)*V.Y[indx_x1_P1_x2_P2])/(h2_2)
			-(1./12)*(-(1./12)*V.Y[indx_x1_P2_x2_M2] +(4./3)*V.Y[indx_x1_P2_x2_M1] -(5./2)*V.Y[indx_x1_P2] + (4./3)*V.Y[indx_x1_P2_x2_P1] - (1./12)*V.Y[indx_x1_P2_x2_P2])/(h2_2)
		      )/(h1);
		      
		Y_d2a1_2a2 = ( 
			(-(1./12)*V.Y[indx_x1_M1_x2_M2] +(4./3)*V.Y[indx_x1_M1_x2_M1] -(5./2)*V.Y[indx_x1_M1] + (4./3)*V.Y[indx_x1_M1_x2_P1] - (1./12)*V.Y[indx_x1_M1_x2_P2])/(h2_2) 
		     - 2*(-(1./12)*V.Y[indx_x2_M2] +(4./3)*V.Y[indx_x2_M1] -(5./2)*V.Y[indx] + (4./3)*V.Y[indx_x2_P1] - (1./12)*V.Y[indx_x2_P2])/(h2_2) 
		     + (-(1./12)*V.Y[indx_x1_P1_x2_M2] +(4./3)*V.Y[indx_x1_P1_x2_M1] -(5./2)*V.Y[indx_x1_P1] + (4./3)*V.Y[indx_x1_P1_x2_P1] - (1./12)*V.Y[indx_x1_P1_x2_P2])/(h2_2) 
		      )/(h1_2);
	}
	else if (FD_ORDER == 6){
		X_da1a2 =( 
			-(1./60)*( -(1./60)*V.X[indx_x1_M3_x2_M3] + (3./20)*V.X[indx_x1_M3_x2_M2] -(3./4)*V.X[indx_x1_M3_x2_M1] + (3./4)*V.X[indx_x1_M3_x2_P1] - (3./20)*V.X[indx_x1_M3_x2_P2]+ (1./60)*V.X[indx_x1_M3_x2_P3])/(h2) 
		        +(3./20)*( -(1./60)*V.X[indx_x1_M2_x2_M3] + (3./20)*V.X[indx_x1_M2_x2_M2] -(3./4)*V.X[indx_x1_M2_x2_M1] + (3./4)*V.X[indx_x1_M2_x2_P1] - (3./20)*V.X[indx_x1_M2_x2_P2]+ (1./60)*V.X[indx_x1_M2_x2_P3])/(h2)
		        -(3./4)*( -(1./60)*V.X[indx_x1_M1_x2_M3] + (3./20)*V.X[indx_x1_M1_x2_M2] -(3./4)*V.X[indx_x1_M1_x2_M1] + (3./4)*V.X[indx_x1_M1_x2_P1] - (3./20)*V.X[indx_x1_M1_x2_P2]+ (1./60)*V.X[indx_x1_M1_x2_P3])/(h2) 
		        +(3./4)*( -(1./60)*V.X[indx_x1_P1_x2_M3] + (3./20)*V.X[indx_x1_P1_x2_M2] -(3./4)*V.X[indx_x1_P1_x2_M1] + (3./4)*V.X[indx_x1_P1_x2_P1] - (3./20)*V.X[indx_x1_P1_x2_P2]+ (1./60)*V.X[indx_x1_P1_x2_P3])/(h2) 
		        -(3./20)*( -(1./60)*V.X[indx_x1_P2_x2_M3] + (3./20)*V.X[indx_x1_P2_x2_M2] -(3./4)*V.X[indx_x1_P2_x2_M1] + (3./4)*V.X[indx_x1_P2_x2_P1] - (3./20)*V.X[indx_x1_P2_x2_P2]+ (1./60)*V.X[indx_x1_P2_x2_P3])/(h2)
		        +(1./60)*( -(1./60)*V.X[indx_x1_P3_x2_M3] + (3./20)*V.X[indx_x1_P3_x2_M2] -(3./4)*V.X[indx_x1_P3_x2_M1] + (3./4)*V.X[indx_x1_P3_x2_P1] - (3./20)*V.X[indx_x1_P3_x2_P2]+ (1./60)*V.X[indx_x1_P3_x2_P3])/(h2)
		    )/(h1);	
		
		X_d2a1a2=(
		       (1./90)*( -(1./60)*V.X[indx_x1_M3_x2_M3] + (3./20)*V.X[indx_x1_M3_x2_M2] -(3./4)*V.X[indx_x1_M3_x2_M1] + (3./4)*V.X[indx_x1_M3_x2_P1] - (3./20)*V.X[indx_x1_M3_x2_P2]+ (1./60)*V.X[indx_x1_M3_x2_P3])/(h2)
		      -(3./20)*( -(1./60)*V.X[indx_x1_M2_x2_M3] + (3./20)*V.X[indx_x1_M2_x2_M2] -(3./4)*V.X[indx_x1_M2_x2_M1] + (3./4)*V.X[indx_x1_M2_x2_P1] - (3./20)*V.X[indx_x1_M2_x2_P2]+ (1./60)*V.X[indx_x1_M2_x2_P3])/(h2) 
		      +(3./2)*( -(1./60)*V.X[indx_x1_M1_x2_M3] + (3./20)*V.X[indx_x1_M1_x2_M2] -(3./4)*V.X[indx_x1_M1_x2_M1] + (3./4)*V.X[indx_x1_M1_x2_P1] - (3./20)*V.X[indx_x1_M1_x2_P2]+ (1./60)*V.X[indx_x1_M1_x2_P3])/(h2) 
		      -(49./18)*( -(1./60)*V.X[indx_x2_M3] + (3./20)*V.X[indx_x2_M2] -(3./4)*V.X[indx_x2_M1] + (3./4)*V.X[indx_x2_P1] - (3./20)*V.X[indx_x2_P2]+ (1./60)*V.X[indx_x2_P3])/(h2)
		      +(3./2)*( -(1./60)*V.X[indx_x1_P1_x2_M3] + (3./20)*V.X[indx_x1_P1_x2_M2] -(3./4)*V.X[indx_x1_P1_x2_M1] + (3./4)*V.X[indx_x1_P1_x2_P1] - (3./20)*V.X[indx_x1_P1_x2_P2]+ (1./60)*V.X[indx_x1_P1_x2_P3])/(h2) 
		      -(3./20)*( -(1./60)*V.X[indx_x1_P2_x2_M3] + (3./20)*V.X[indx_x1_P2_x2_M2] -(3./4)*V.X[indx_x1_P2_x2_M1] + (3./4)*V.X[indx_x1_P2_x2_P1] - (3./20)*V.X[indx_x1_P2_x2_P2]+ (1./60)*V.X[indx_x1_P2_x2_P3])/(h2) 
		      +(1./90)*( -(1./60)*V.X[indx_x1_P3_x2_M3] + (3./20)*V.X[indx_x1_P3_x2_M2] -(3./4)*V.X[indx_x1_P3_x2_M1] + (3./4)*V.X[indx_x1_P3_x2_P1] - (3./20)*V.X[indx_x1_P3_x2_P2]+ (1./60)*V.X[indx_x1_P3_x2_P3])/(h2)
		    )/(h1_2);
		    
	       X_d2a2a1 =( 
			-(1./60)*((1./90)*V.X[indx_x1_M3_x2_M3]-(3./20)*V.X[indx_x1_M3_x2_M2] + (3./2)*V.X[indx_x1_M3_x2_M1] -(49./18)*V.X[indx_x1_M3] + (3./2)*V.X[indx_x1_M3_x2_P1] - (3./20)*V.X[indx_x1_M3_x2_P2] + (1./90)*V.X[indx_x1_M3_x2_P3])/(h2_2) 
			+(3./20)*((1./90)*V.X[indx_x1_M2_x2_M3]-(3./20)*V.X[indx_x1_M2_x2_M2] + (3./2)*V.X[indx_x1_M2_x2_M1] -(49./18)*V.X[indx_x1_M2] + (3./2)*V.X[indx_x1_M2_x2_P1] - (3./20)*V.X[indx_x1_M2_x2_P2] + (1./90)*V.X[indx_x1_M2_x2_P3])/(h2_2) 
			-(3./4)*((1./90)*V.X[indx_x1_M1_x2_M3]-(3./20)*V.X[indx_x1_M1_x2_M2] + (3./2)*V.X[indx_x1_M1_x2_M1] -(49./18)*V.X[indx_x1_M1] + (3./2)*V.X[indx_x1_M1_x2_P1] - (3./20)*V.X[indx_x1_M1_x2_P2] + (1./90)*V.X[indx_x1_M1_x2_P3])/(h2_2)
			+(3./4)*((1./90)*V.X[indx_x1_P1_x2_M3]-(3./20)*V.X[indx_x1_P1_x2_M2] + (3./2)*V.X[indx_x1_P1_x2_M1] -(49./18)*V.X[indx_x1_P1] + (3./2)*V.X[indx_x1_P1_x2_P1] - (3./20)*V.X[indx_x1_P1_x2_P2] + (1./90)*V.X[indx_x1_P1_x2_P3])/(h2_2) 
			-(3./20)*((1./90)*V.X[indx_x1_P2_x2_M3]-(3./20)*V.X[indx_x1_P2_x2_M2] + (3./2)*V.X[indx_x1_P2_x2_M1] -(49./18)*V.X[indx_x1_P2] + (3./2)*V.X[indx_x1_P2_x2_P1] - (3./20)*V.X[indx_x1_P2_x2_P2] + (1./90)*V.X[indx_x1_P2_x2_P3])/(h2_2)
			+(1./60)*((1./90)*V.X[indx_x1_P3_x2_M3]-(3./20)*V.X[indx_x1_P3_x2_M2] + (3./2)*V.X[indx_x1_P3_x2_M1] -(49./18)*V.X[indx_x1_P3] + (3./2)*V.X[indx_x1_P3_x2_P1] - (3./20)*V.X[indx_x1_P3_x2_P2] + (1./90)*V.X[indx_x1_P3_x2_P3])/(h2_2)
		     )/(h1);
		
		X_d2a1_2a2=(
			(1./90)*((1./90)*V.X[indx_x1_M3_x2_M3]-(3./20)*V.X[indx_x1_M3_x2_M2] + (3./2)*V.X[indx_x1_M3_x2_M1] -(49./18)*V.X[indx_x1_M3] + (3./2)*V.X[indx_x1_M3_x2_P1] - (3./20)*V.X[indx_x1_M3_x2_P2] + (1./90)*V.X[indx_x1_M3_x2_P3])/(h2_2)
		       -(3./20)*((1./90)*V.X[indx_x1_M2_x2_M3]-(3./20)*V.X[indx_x1_M2_x2_M2] + (3./2)*V.X[indx_x1_M2_x2_M1] -(49./18)*V.X[indx_x1_M2] + (3./2)*V.X[indx_x1_M2_x2_P1] - (3./20)*V.X[indx_x1_M2_x2_P2] + (1./90)*V.X[indx_x1_M2_x2_P3])/(h2_2) 
		       +(3./2)*((1./90)*V.X[indx_x1_M1_x2_M3]-(3./20)*V.X[indx_x1_M1_x2_M2] + (3./2)*V.X[indx_x1_M1_x2_M1] -(49./18)*V.X[indx_x1_M1] + (3./2)*V.X[indx_x1_M1_x2_P1] - (3./20)*V.X[indx_x1_M1_x2_P2] + (1./90)*V.X[indx_x1_M1_x2_P3])/(h2_2) 
		       -(49./18)*((1./90)*V.X[indx_x2_M3]-(3./20)*V.X[indx_x2_M2] + (3./2)*V.X[indx_x2_M1] -(49./18)*V.X[indx] + (3./2)*V.X[indx_x2_P1] - (3./20)*V.X[indx_x2_P2] + (1./90)*V.X[indx_x2_P3])/(h2_2) 
		       +(3./2)*((1./90)*V.X[indx_x1_P1_x2_M3]-(3./20)*V.X[indx_x1_P1_x2_M2] + (3./2)*V.X[indx_x1_P1_x2_M1] -(49./18)*V.X[indx_x1_P1] + (3./2)*V.X[indx_x1_P1_x2_P1] - (3./20)*V.X[indx_x1_P1_x2_P2] + (1./90)*V.X[indx_x1_P1_x2_P3])/(h2_2) 
		       -(3./20)*((1./90)*V.X[indx_x1_P2_x2_M3]-(3./20)*V.X[indx_x1_P2_x2_M2] + (3./2)*V.X[indx_x1_P2_x2_M1] -(49./18)*V.X[indx_x1_P2] + (3./2)*V.X[indx_x1_P2_x2_P1] - (3./20)*V.X[indx_x1_P2_x2_P2] + (1./90)*V.X[indx_x1_P2_x2_P3])/(h2_2) 
		       +(1./90)*((1./90)*V.X[indx_x1_P3_x2_M3]-(3./20)*V.X[indx_x1_P3_x2_M2] + (3./2)*V.X[indx_x1_P3_x2_M1] -(49./18)*V.X[indx_x1_P3] + (3./2)*V.X[indx_x1_P3_x2_P1] - (3./20)*V.X[indx_x1_P3_x2_P2] + (1./90)*V.X[indx_x1_P3_x2_P3])/(h2_2)
		    )/(h1_2);

		Y_da1a2 =( 
			-(1./60)*( -(1./60)*V.Y[indx_x1_M3_x2_M3] + (3./20)*V.Y[indx_x1_M3_x2_M2] -(3./4)*V.Y[indx_x1_M3_x2_M1] + (3./4)*V.Y[indx_x1_M3_x2_P1] - (3./20)*V.Y[indx_x1_M3_x2_P2]+ (1./60)*V.Y[indx_x1_M3_x2_P3])/(h2) 
		        +(3./20)*( -(1./60)*V.Y[indx_x1_M2_x2_M3] + (3./20)*V.Y[indx_x1_M2_x2_M2] -(3./4)*V.Y[indx_x1_M2_x2_M1] + (3./4)*V.Y[indx_x1_M2_x2_P1] - (3./20)*V.Y[indx_x1_M2_x2_P2]+ (1./60)*V.Y[indx_x1_M2_x2_P3])/(h2)
		        -(3./4)*( -(1./60)*V.Y[indx_x1_M1_x2_M3] + (3./20)*V.Y[indx_x1_M1_x2_M2] -(3./4)*V.Y[indx_x1_M1_x2_M1] + (3./4)*V.Y[indx_x1_M1_x2_P1] - (3./20)*V.Y[indx_x1_M1_x2_P2]+ (1./60)*V.Y[indx_x1_M1_x2_P3])/(h2) 
		        +(3./4)*( -(1./60)*V.Y[indx_x1_P1_x2_M3] + (3./20)*V.Y[indx_x1_P1_x2_M2] -(3./4)*V.Y[indx_x1_P1_x2_M1] + (3./4)*V.Y[indx_x1_P1_x2_P1] - (3./20)*V.Y[indx_x1_P1_x2_P2]+ (1./60)*V.Y[indx_x1_P1_x2_P3])/(h2) 
		        -(3./20)*( -(1./60)*V.Y[indx_x1_P2_x2_M3] + (3./20)*V.Y[indx_x1_P2_x2_M2] -(3./4)*V.Y[indx_x1_P2_x2_M1] + (3./4)*V.Y[indx_x1_P2_x2_P1] - (3./20)*V.Y[indx_x1_P2_x2_P2]+ (1./60)*V.Y[indx_x1_P2_x2_P3])/(h2)
		        +(1./60)*( -(1./60)*V.Y[indx_x1_P3_x2_M3] + (3./20)*V.Y[indx_x1_P3_x2_M2] -(3./4)*V.Y[indx_x1_P3_x2_M1] + (3./4)*V.Y[indx_x1_P3_x2_P1] - (3./20)*V.Y[indx_x1_P3_x2_P2]+ (1./60)*V.Y[indx_x1_P3_x2_P3])/(h2)
		    )/(h1);	
		
		Y_d2a1a2=(
		       (1./90)*( -(1./60)*V.Y[indx_x1_M3_x2_M3] + (3./20)*V.Y[indx_x1_M3_x2_M2] -(3./4)*V.Y[indx_x1_M3_x2_M1] + (3./4)*V.Y[indx_x1_M3_x2_P1] - (3./20)*V.Y[indx_x1_M3_x2_P2]+ (1./60)*V.Y[indx_x1_M3_x2_P3])/(h2)
		      -(3./20)*( -(1./60)*V.Y[indx_x1_M2_x2_M3] + (3./20)*V.Y[indx_x1_M2_x2_M2] -(3./4)*V.Y[indx_x1_M2_x2_M1] + (3./4)*V.Y[indx_x1_M2_x2_P1] - (3./20)*V.Y[indx_x1_M2_x2_P2]+ (1./60)*V.Y[indx_x1_M2_x2_P3])/(h2) 
		      +(3./2)*( -(1./60)*V.Y[indx_x1_M1_x2_M3] + (3./20)*V.Y[indx_x1_M1_x2_M2] -(3./4)*V.Y[indx_x1_M1_x2_M1] + (3./4)*V.Y[indx_x1_M1_x2_P1] - (3./20)*V.Y[indx_x1_M1_x2_P2]+ (1./60)*V.Y[indx_x1_M1_x2_P3])/(h2) 
		      -(49./18)*( -(1./60)*V.Y[indx_x2_M3] + (3./20)*V.Y[indx_x2_M2] -(3./4)*V.Y[indx_x2_M1] + (3./4)*V.Y[indx_x2_P1] - (3./20)*V.Y[indx_x2_P2]+ (1./60)*V.Y[indx_x2_P3])/(h2)
		      +(3./2)*( -(1./60)*V.Y[indx_x1_P1_x2_M3] + (3./20)*V.Y[indx_x1_P1_x2_M2] -(3./4)*V.Y[indx_x1_P1_x2_M1] + (3./4)*V.Y[indx_x1_P1_x2_P1] - (3./20)*V.Y[indx_x1_P1_x2_P2]+ (1./60)*V.Y[indx_x1_P1_x2_P3])/(h2) 
		      -(3./20)*( -(1./60)*V.Y[indx_x1_P2_x2_M3] + (3./20)*V.Y[indx_x1_P2_x2_M2] -(3./4)*V.Y[indx_x1_P2_x2_M1] + (3./4)*V.Y[indx_x1_P2_x2_P1] - (3./20)*V.Y[indx_x1_P2_x2_P2]+ (1./60)*V.Y[indx_x1_P2_x2_P3])/(h2) 
		      +(1./90)*( -(1./60)*V.Y[indx_x1_P3_x2_M3] + (3./20)*V.Y[indx_x1_P3_x2_M2] -(3./4)*V.Y[indx_x1_P3_x2_M1] + (3./4)*V.Y[indx_x1_P3_x2_P1] - (3./20)*V.Y[indx_x1_P3_x2_P2]+ (1./60)*V.Y[indx_x1_P3_x2_P3])/(h2)
		    )/(h1_2);
		    
	    Y_d2a2a1 =( 
			-(1./60)*((1./90)*V.Y[indx_x1_M3_x2_M3]-(3./20)*V.Y[indx_x1_M3_x2_M2] + (3./2)*V.Y[indx_x1_M3_x2_M1] -(49./18)*V.Y[indx_x1_M3] + (3./2)*V.Y[indx_x1_M3_x2_P1] - (3./20)*V.Y[indx_x1_M3_x2_P2] + (1./90)*V.Y[indx_x1_M3_x2_P3])/(h2_2) 
			+(3./20)*((1./90)*V.Y[indx_x1_M2_x2_M3]-(3./20)*V.Y[indx_x1_M2_x2_M2] + (3./2)*V.Y[indx_x1_M2_x2_M1] -(49./18)*V.Y[indx_x1_M2] + (3./2)*V.Y[indx_x1_M2_x2_P1] - (3./20)*V.Y[indx_x1_M2_x2_P2] + (1./90)*V.Y[indx_x1_M2_x2_P3])/(h2_2) 
			-(3./4)*((1./90)*V.Y[indx_x1_M1_x2_M3]-(3./20)*V.Y[indx_x1_M1_x2_M2] + (3./2)*V.Y[indx_x1_M1_x2_M1] -(49./18)*V.Y[indx_x1_M1] + (3./2)*V.Y[indx_x1_M1_x2_P1] - (3./20)*V.Y[indx_x1_M1_x2_P2] + (1./90)*V.Y[indx_x1_M1_x2_P3])/(h2_2)
			+(3./4)*((1./90)*V.Y[indx_x1_P1_x2_M3]-(3./20)*V.Y[indx_x1_P1_x2_M2] + (3./2)*V.Y[indx_x1_P1_x2_M1] -(49./18)*V.Y[indx_x1_P1] + (3./2)*V.Y[indx_x1_P1_x2_P1] - (3./20)*V.Y[indx_x1_P1_x2_P2] + (1./90)*V.Y[indx_x1_P1_x2_P3])/(h2_2) 
			-(3./20)*((1./90)*V.Y[indx_x1_P2_x2_M3]-(3./20)*V.Y[indx_x1_P2_x2_M2] + (3./2)*V.Y[indx_x1_P2_x2_M1] -(49./18)*V.Y[indx_x1_P2] + (3./2)*V.Y[indx_x1_P2_x2_P1] - (3./20)*V.Y[indx_x1_P2_x2_P2] + (1./90)*V.Y[indx_x1_P2_x2_P3])/(h2_2)
			+(1./60)*((1./90)*V.Y[indx_x1_P3_x2_M3]-(3./20)*V.Y[indx_x1_P3_x2_M2] + (3./2)*V.Y[indx_x1_P3_x2_M1] -(49./18)*V.Y[indx_x1_P3] + (3./2)*V.Y[indx_x1_P3_x2_P1] - (3./20)*V.Y[indx_x1_P3_x2_P2] + (1./90)*V.Y[indx_x1_P3_x2_P3])/(h2_2)
		     )/(h1);
		
		Y_d2a1_2a2=(
			(1./90)*((1./90)*V.Y[indx_x1_M3_x2_M3]-(3./20)*V.Y[indx_x1_M3_x2_M2] + (3./2)*V.Y[indx_x1_M3_x2_M1] -(49./18)*V.Y[indx_x1_M3] + (3./2)*V.Y[indx_x1_M3_x2_P1] - (3./20)*V.Y[indx_x1_M3_x2_P2] + (1./90)*V.Y[indx_x1_M3_x2_P3])/(h2_2)
		       -(3./20)*((1./90)*V.Y[indx_x1_M2_x2_M3]-(3./20)*V.Y[indx_x1_M2_x2_M2] + (3./2)*V.Y[indx_x1_M2_x2_M1] -(49./18)*V.Y[indx_x1_M2] + (3./2)*V.Y[indx_x1_M2_x2_P1] - (3./20)*V.Y[indx_x1_M2_x2_P2] + (1./90)*V.Y[indx_x1_M2_x2_P3])/(h2_2) 
		       +(3./2)*((1./90)*V.Y[indx_x1_M1_x2_M3]-(3./20)*V.Y[indx_x1_M1_x2_M2] + (3./2)*V.Y[indx_x1_M1_x2_M1] -(49./18)*V.Y[indx_x1_M1] + (3./2)*V.Y[indx_x1_M1_x2_P1] - (3./20)*V.Y[indx_x1_M1_x2_P2] + (1./90)*V.Y[indx_x1_M1_x2_P3])/(h2_2) 
		       -(49./18)*((1./90)*V.Y[indx_x2_M3]-(3./20)*V.Y[indx_x2_M2] + (3./2)*V.Y[indx_x2_M1] -(49./18)*V.Y[indx] + (3./2)*V.Y[indx_x2_P1] - (3./20)*V.Y[indx_x2_P2] + (1./90)*V.Y[indx_x2_P3])/(h2_2) 
		       +(3./2)*((1./90)*V.Y[indx_x1_P1_x2_M3]-(3./20)*V.Y[indx_x1_P1_x2_M2] + (3./2)*V.Y[indx_x1_P1_x2_M1] -(49./18)*V.Y[indx_x1_P1] + (3./2)*V.Y[indx_x1_P1_x2_P1] - (3./20)*V.Y[indx_x1_P1_x2_P2] + (1./90)*V.Y[indx_x1_P1_x2_P3])/(h2_2) 
		       -(3./20)*((1./90)*V.Y[indx_x1_P2_x2_M3]-(3./20)*V.Y[indx_x1_P2_x2_M2] + (3./2)*V.Y[indx_x1_P2_x2_M1] -(49./18)*V.Y[indx_x1_P2] + (3./2)*V.Y[indx_x1_P2_x2_P1] - (3./20)*V.Y[indx_x1_P2_x2_P2] + (1./90)*V.Y[indx_x1_P2_x2_P3])/(h2_2) 
		       +(1./90)*((1./90)*V.Y[indx_x1_P3_x2_M3]-(3./20)*V.Y[indx_x1_P3_x2_M2] + (3./2)*V.Y[indx_x1_P3_x2_M1] -(49./18)*V.Y[indx_x1_P3] + (3./2)*V.Y[indx_x1_P3_x2_P1] - (3./20)*V.Y[indx_x1_P3_x2_P2] + (1./90)*V.Y[indx_x1_P3_x2_P3])/(h2_2)
		    )/(h1_2);

	}
	
	if(strcmp( par.grid_x1,"Fourier" ) ==0){
	  
	    if(strcmp( par.grid_x2,"Fourier" ) ==0){
	    	V.X12[indx]=X_da1a2;
	    	V.Y12[indx]=Y_da1a2;
	    }
	    else{
	      if( fabs(fabs(cos2)-1.) < TINY){
			V.X12[indx] = -X_d2a2a1/cos2;
			V.Y12[indx] = -Y_d2a2a1/cos2;
	      }
	      else{
			V.X12[indx] = -X_da1a2/sin2;
			V.Y12[indx] = -Y_da1a2/sin2;
	      }
	    }
	}
	else if(strcmp( par.grid_x2,"Fourier" ) ==0){
	  
	    if(strcmp( par.grid_x1,"Fourier" ) ==0){
	    	V.X12[indx] = X_da1a2;
	    	V.Y12[indx] = Y_da1a2;
	    }
	    else{
	      if( fabs(fabs(cos1)-1.) < TINY){
	      	V.X12[indx] = - X_d2a1a2/cos1;
	      	V.Y12[indx] = - Y_d2a1a2/cos1;
	      }
	      else{
	      	V.X12[indx] = -X_da1a2/sin1;
	      	V.Y12[indx] = -Y_da1a2/sin1;
	      }
	    }
	}
	else{
	    if( fabs(fabs(cos1)-1.) > TINY && fabs(fabs(cos2)-1.) > TINY ){
	      V.X12[indx] = X_da1a2/(sin1*sin2);
	      V.Y12[indx] = Y_da1a2/(sin1*sin2);
	    }
	    else if( fabs(fabs(cos1)-1.) < TINY){
	      V.X12[indx] = X_d2a1a2/(cos1*sin2);
	      V.Y12[indx] = Y_d2a1a2/(cos1*sin2);
	    }
	    else if( fabs(fabs(cos2)-1.) < TINY){
	      V.X12[indx] = X_d2a2a1/(cos2*sin1);
	      V.Y12[indx] = Y_d2a2a1/(cos2*sin1);
	    }
	    else
	    {
	      V.X12[indx] = X_d2a1_2a2/(cos1*cos2);	
	      V.Y12[indx] = Y_d2a2a1/(cos2*sin1);	    
	    }
	}

	return;
}
// // -------------------------------------------------------------------------------