#include "2+1_Free_Boundaries.h"

void Func_Delta(double kappa, double R, double *Delta, double *Delta_R){
	double k2=sqr(kappa);

	*Delta   = (1.-R)*(1-k2*R);
	*Delta_R = -(1-k2*R) -k2*(1.-R);

	return;
}
//----------------------------------------------------------------------------------------------------
void Func_Sigma(double kappa, double R, double x, double *Sigma, double *Sigma_R, double *Sigma_x){
	double k2=sqr(kappa), R2=sqr(R), x2=sqr(x);

	*Sigma   = 1. + k2*R2*x2;
	*Sigma_R = 2*k2*R*x2;
	*Sigma_x = 2*k2*R2*x;

	return;
}
//----------------------------------------------------------------------------------------------------
void Func_Sigma0(double kappa, double R, double *Sigma0, double *Sigma0_R){
	
	double k2=sqr(kappa), R2=sqr(R);
	*Sigma0   = 1. + k2*R2;
	*Sigma0_R = 2*k2*R;
	
	return;
}
//----------------------------------------------------------------------------------------------------
void Func_R(double t, double r, double *R, double *R_r, double *R_t){
	
	*R   =r*(1.-t);
	*R_r = 1.-t;
	*R_t = -r;

	return;
}

//----------------------------------------------------------------------------------------------------
void Func_C(double kappa, double r, double *C, double *C_r){

	double Delta, Delta_r, Sigma0, Sigma0_r;

	Func_Delta(kappa, r, &Delta, &Delta_r);
	Func_Sigma0(kappa, r, &Sigma0, &Sigma0_r);
	
	*C   = Sigma0/Delta;
	*C_r = Sigma0_r/Delta - Sigma0*Delta_r/sqr(Delta);
	
	return;
}
//----------------------------------------------------------------------------------------------------
void Load_ExactSolutions(parameters par, double ***c_f2[], double ***c_df2[]){

	int ik, i0, i1, i2;
	
	double  r1=par.r1;

	if(r1!=r1_read){
		fprintf(stderr, "Error in Load_ExactSolutions: function generation in Mathematica has r1=%lf, code has par.r1=%lf\n",r1_read,r1 );
		exit(1);
	}

	FILE *fr_cf2, *fr_cdf2;
	char fn_cf2[200], fn_cdf2[200];

	sprintf(fn_cf2,"ChebCoeff_ExactSolution/cf2_l%d.dat", par.ell); fr_cf2=fopen(fn_cf2,"r");
	sprintf(fn_cdf2,"ChebCoeff_ExactSolution/cdf2_l%d.dat", par.ell); fr_cdf2=fopen(fn_cdf2,"r");
	

	if(fr_cf2 == NULL || fr_cdf2 == NULL){
		fprintf(stderr, "Error in Load ExactSolutions: Cannot open file\n");
		exit(1);
	}

	for(ik=0; ik<=Nkappa_read; ik++){
		for(i0=0; i0<=N0_read; i0++){
			for(i1=0; i1<=N1_read; i1++){
				for(i2=0; i2<=N2_read; i2++){
					long double cf2_read, cdf2_read;
					fscanf(fr_cf2, "%Lf", &cf2_read);
					fscanf(fr_cdf2, "%Lf", &cdf2_read);	

					

					c_f2[ik][i0][i1][i2]=cf2_read;
					c_df2[ik][i0][i1][i2]=cdf2_read;					
				}
			}
		}
	}
	

	return;
}
//--------------------------------------------------------
void ExactSolution_interpolate(parameters par, double x0, double x1, double x2, double ***c_f2[], double ***c_df2[], double *f1, double *f2, double *df1, double *df2){

	
	int ik;
	double *f2_vec=dvector(0,Nkappa_read), *cf2_vec=dvector(0,Nkappa_read),
		   *df2_vec=dvector(0,Nkappa_read), *cdf2_vec=dvector(0,Nkappa_read),
		   *f2Norm_vec=dvector(0,Nkappa_read), *cf2Norm_vec=dvector(0,Nkappa_read);

	for(ik=0; ik<=Nkappa_read; ik++){
		f2Norm_vec[ik]   = Clenshaw_Chebyshev_3D( c_f2[ik],  N0_read, N1_read, N2_read, -1., 1., 1.);
		f2_vec[ik]   = Clenshaw_Chebyshev_3D( c_f2[ik],  N0_read, N1_read, N2_read, x0, x1, x2);
		df2_vec[ik]  = Clenshaw_Chebyshev_3D(c_df2[ik],  N0_read, N1_read, N2_read, x0, x1, x2);
	}
	Chebyshev_Coefficients( f2_vec,  cf2_vec, Nkappa_read, "Lobatto");
	Chebyshev_Coefficients(df2_vec, cdf2_vec, Nkappa_read, "Lobatto");
	Chebyshev_Coefficients( f2Norm_vec,  cf2Norm_vec, Nkappa_read, "Lobatto");

	//-----
	// FILE *fp=fopen("cheb_kappa.dat", "w");
	// for(ik=0; ik<=Nkappa_read; ik++) fprintf(fp, "%d\t%3.15e\t%3.15e\n", ik,f2_vec[ik], cf2_vec[ik] );
	// fclose(fp);
	// exit(1);
	//-------

	double kappa= par.kappa, xkappa = 2*kappa/0.9-1., f2Norm;
	
	f2Norm = Clenshaw_Chebyshev( cf2Norm_vec,  Nkappa_read, xkappa);
	*f2    = Clenshaw_Chebyshev( cf2_vec,  Nkappa_read, xkappa)/f2Norm;
	*df2   = Clenshaw_Chebyshev(cdf2_vec,  Nkappa_read, xkappa)/f2Norm;


	double t=func_t_of_x0(par, x0), one_minus_t=1.-t, one_minus_k4=1.-pow(kappa,4), one_minus_k2 = 1. - sqr(kappa), one_minus_k2_3=pow(one_minus_k2,3);
	switch(par.ell){		
			case 0:
				*f1  = 0.5*one_minus_t*one_minus_k4/f2Norm;
				*df1 = -0.5*one_minus_k4/f2Norm;
				break;

			case 1:
				*f1  =  one_minus_t*one_minus_k2_3/12./f2Norm;
				*df1 = -one_minus_k2_3/12./f2Norm;
				break;

			case 2:
				*f1  = 0.;
				*df1 = 0.;
				break;

			case 3:
				*f1  = 0.;
				*df1 = 0.;
				break;
		}


	free_dvector( f2_vec, 0, Nkappa_read);    free_dvector( cf2_vec, 0, Nkappa_read);
	free_dvector(df2_vec, 0, Nkappa_read);    free_dvector(cdf2_vec, 0, Nkappa_read);
	free_dvector(f2Norm_vec, 0, Nkappa_read); free_dvector(cf2Norm_vec, 0, Nkappa_read);
	return;
}
//--------------------------------------------------------------------------------
void get_ExactSolution_ID(parameters *par){
	int ik, i1, i2, N1=(*par).N1, N2=(*par).N2;
	double f1, df1, f2, df2;


	double ***c_f2[Nkappa_read+1], ***c_df2[Nkappa_read+1]; 
	for(ik=0; ik<=Nkappa_read; ik++){
		c_f2[ik]  = d3tensor(0, N0_read, 0, N1_read, 0, N2_read);
		c_df2[ik] = d3tensor(0, N0_read, 0, N1_read, 0, N2_read);
	}
	Load_ExactSolutions(*par, c_f2, c_df2);


	for(i1=0; i1<=N1; i1++){
			double x1=(*par).grid_points_x1[i1];
				   
			
			for(i2=0; i2<=N2; i2++){
				double x2=(*par).grid_points_x2[i2];
					   

				int I0_slice = Index_Slice(*par, 0, i1, i2),
					I1_slice = Index_Slice(*par , 1, i1, i2);
				int I1D_0_slice = Index_Slice_Fields1D(0,(*par).nslice),
					I1D_1_slice = Index_Slice_Fields1D(1,(*par).nslice);	

					ExactSolution_interpolate(*par, -1., x1 , x2, c_f2, c_df2, &f1, &f2, &df1, &df2);
					// ID for 2+1 Fields					
					
					(*par).ID.X[I0_slice] = f2;
					(*par).ID.X[I1_slice] = df2;

					if(i1==0){
						// ID for 0+1 FIELDS				
						(*par).ID.X[I1D_0_slice] = f1;
						(*par).ID.X[I1D_1_slice] = df1;
					}		
			}		
		}	
	Derivatives_Slice(*par, (*par).ID);
	for(ik=0; ik<=Nkappa_read; ik++){
		free_d3tensor(c_f2[ik],  0, N0_read, 0, N1_read, 0, N2_read);
		free_d3tensor(c_df2[ik],  0, N0_read, 0, N1_read, 0, N2_read);
	}
	return;
}
//------------------------------------------------------------------------------
void get_ExactSolution(parameters par, double *X){
	int ik, i0, i1, i2, N0=par.N0, N1=par.N1, N2=par.N2;
	double f1, df1, f2, df2;


	double ***c_f2[Nkappa_read+1], ***c_df2[Nkappa_read+1]; 
	for(ik=0; ik<=Nkappa_read; ik++){
		c_f2[ik]  = d3tensor(0, N0_read, 0, N1_read, 0, N2_read);
		c_df2[ik] = d3tensor(0, N0_read, 0, N1_read, 0, N2_read);
	}
	Load_ExactSolutions(par, c_f2, c_df2);


	for(i1=0; i1<=N1; i1++){
			double x1=par.grid_points_x1[i1];
				   
			
			for(i2=0; i2<=N2; i2++){
				double x2=par.grid_points_x2[i2];
					   

				int I0_slice = Index_Slice(par, 0, i1, i2),
					I1_slice = Index_Slice(par , 1, i1, i2);
				int I1D_0_slice = Index_Slice_Fields1D(0,par.nslice),
					I1D_1_slice = Index_Slice_Fields1D(1,par.nslice);	

				for(i0=0; i0<=N0; i0++){
					double  x0 = par.grid_points_x0[i0];
					
					int I0 = Index(0, i0, N0, i1, N1, i2, N2), 
						I1 = Index(1, i0, N0, i1, N1, i2, N2);
					int I1D_0 = Index_Fields1D(0, par.ngrid, i0, N0), 
						I1D_1 = Index_Fields1D(1, par.ngrid, i0, N0);

					ExactSolution_interpolate(par, x0, x1 , x2, c_f2, c_df2, &f1, &f2, &df1, &df2);
					
					X[I0]= (  f2 - par.ID.X[I0_slice] )/(1.+x0);
					X[I1]= ( df2 - par.ID.X[I1_slice] )/(1.+x0);

					if(i1==0){
						X[I1D_0]= (  f1 - par.ID.X[I1D_0_slice] )/(1.+x0);
						X[I1D_1]= ( df1 - par.ID.X[I1D_1_slice] )/(1.+x0);
					}
				}			
			}		
		}
	for(ik=0; ik<=Nkappa_read; ik++){
		free_d3tensor(c_f2[ik],  0, N0_read, 0, N1_read, 0, N2_read);
		free_d3tensor(c_df2[ik],  0, N0_read, 0, N1_read, 0, N2_read);
	}
	return;
}
//------------------------------------------------------------------------------

void get_LinearEquationCoefficients(double kappa, double t, double r, double x, 
	double *S, double *L1_0, double *L1_r, double *L1_x,
	double *L2_0, double *L2_r, double *L2_x, double *L2_rr, double *L2_rx, double *L2_xx){


	double r2=sqr(r), r4=sqr(r2), R, R2, R_r, R_t, Sigma0, Sigma0_R, Sigma, Sigma_R, Sigma_x, Delta, Delta_R, kappa2= sqr(kappa), one_minus_t=1.-t, one_minus_t_2 = sqr(one_minus_t), x2=sqr(x), sin2 = 1-x2, C, C_r, C2;

	Func_R(t, r, &R, &R_r, &R_t);
	R2=sqr(R);
	Func_Delta(kappa, R, &Delta, &Delta_R);
	Func_Sigma(kappa, R, x, &Sigma, &Sigma_R, &Sigma_x);
	Func_Sigma0(kappa, R, &Sigma0, &Sigma0_R);
	Func_C(kappa, r, &C, &C_r);
	C2 = sqr(C);

	*S = ( R2*kappa2*sin2/C2 + one_minus_t_2*Delta - 2*one_minus_t*Sigma0/C) / Sigma;
	

	*L1_r =  -2*r*( Sigma0 - R*r*kappa2*sin2/C  )/(C*Sigma);
	
	*L1_0 = ( -( 2*Delta + R*Delta_R  )*one_minus_t  +    (2* (1 - 2*Sigma0 + R*Sigma0_R )  + kappa2*(4 - r*C_r/C)*r*R*sin2/C  )/C  )/Sigma;
	
	*L1_x = 0;

	*L2_0 = (-R + kappa2*( -R*( 1 -2*R - 4*r/C) +( 6 - 2*r*C_r/C  )*r2*sin2/C2   ) )/Sigma;	
	
	*L2_r = r2*kappa2*(2*R   +(6 - r*C_r/C  )*r*sin2/C  )/(C*Sigma);  
	
	*L2_rr = r4*kappa2*sin2/(C2*Sigma);
	
	*L2_xx = (1. - x2) / Sigma;
	
	*L2_x = -2 * x / Sigma;
	
	*L2_rx = 0;

	return;
}
//--------------------------------------------------------
double get_Source_f1(parameters par, double t){

	int l=par.ell, l2=l*l;
	double Pl=plgndr(l, 0, t), Plp1=plgndr(l+1, 0, t), tp1=1.+t, tp1_2=sqr(tp1), t2=sqr(t), t3=t2*t, 
		   alpha=par.alpha, kappa=par.kappa, kappa_2=sqr(kappa), kappa_2_p1=1.+kappa_2, Source_f1;
	
	
	if(l==0 && t==1){
		Source_f1 = 1.;
		
	}
	else{
		Source_f1 = -Plp1*(l+1)*(1- 2*t - t2)/tp1_2 + Pl*( -(1+t2) - (1+2*t+3*t2)*l -(1+3*t+t2-t3)*l2     )/tp1_2;
	}
	
	return Source_f1*alpha*kappa_2_p1;

}
//-------------------------------------------------------------------------------
double get_Source_f2(parameters par, double t, double r, double x,  derivs W1D_t ){

	int l = par.ell;
	double alpha = par.alpha, kappa=par.kappa, 
		   f1 = W1D_t.X[0],
	  	   df1 = W1D_t.X[1], df1_t = W1D_t.Y[1], Source_f2;
		

	
	double S, L1_0, L2_0, L0_A, L0_B, Pl=plgndr(l, 0, t), Plp1=plgndr(l+1, 0, t), k2=sqr(kappa), k4=sqr(k2), k6=k4*k2,
		e=1.-t, e2=sqr(e), e3=e2*e, x2=sqr(x), sin2=1-x2,
		R, R_r, R_t, R2, Delta, Delta_R, Sigma, Sigma_R, Sigma_x, Sigma0, Sigma0_R, C, C2, C_r, delta, delta_r;



	Func_R(t, r, &R, &R_r, &R_t);
	R2=sqr(R);
	Func_Delta(kappa, R, &Delta, &Delta_R);
	Func_Delta(kappa, r, &delta, &delta_r);
	Func_Sigma(kappa, R, x, &Sigma, &Sigma_R, &Sigma_x);
	Func_Sigma0(kappa, R, &Sigma0, &Sigma0_R);
	Func_C(kappa, r, &C, &C_r);
	C2 = sqr(C);

	S =  e*( R*k2*sin2/(C2*Sigma) + 2*(1.+k2)/(C*delta) ) - 2*k2*sin2*R*e2/(C*Sigma) + e3*( k2*sin2*R - (1.+k2) )/Sigma; 
	
	L1_0 = ( 2 - r*C_r/C )*R*k2*sin2/(C2*Sigma) + 2*e*(Sigma0_R - k2*R)/(C*Sigma) + e2*( 2. - Delta_R + 2*k2*(1.-R*sin2)   )/Sigma;

	L2_0 = e*( 2*R*k2 -(1.+k2) )/Sigma + ( 2*k2*R + (2 - r*C_r/C)*r*k2*sin2/C   )/(C*Sigma) + l*(l+1)*R*e*k2*x2/Sigma;
	
	L0_A = -(x * x * k6 *  (t * t + 2 * t - 1) *   pow( (-1 + t),  2) *  (1 + k2) * pow(r, 0.7e1) - 0.2e1 *  t * k6 * ( (t * t) + x * x - 0.2e1) * pow(r, 0.6e1) + 0.3e1 * k4 *  (1 + k2) * (  pow( t,  4) * x * x - 0.4e1 *  (t * t) * x * x + (0.5e1 * x * x + 0.1e1 / 0.3e1) *  t - 0.4e1 / 0.3e1 * x * x - 0.1e1 / 0.3e1) * pow(r, 0.5e1) + (-0.6e1 *   pow( t,  3) * k4 + ((- (8 * k2 * k2) + (-k4 - 0.1e1) *  k2) * x * x +  k2 + 0.14e2 * k4 + k6) *  t + ( (2 * k2 * k2) + (k4 + 0.1e1) *  k2) * x * x -  k2 - 0.2e1 * k4 - k6) * pow(r, 0.4e1) + (0.3e1 * x * x *  k2 *  (1 + k2) *   pow( t,  4) - 0.12e2 * x * x *  k2 *  (1 + k2) *  (t * t) + 0.20e2 * x * x *  k2 *  (1 + k2) *  t +  (-3 * k2 * k2 - 3 * k2) * x * x -  (4 * k2) - 0.4e1 * k4) * pow(r, 0.3e1) + (- (6 *  pow( t,  3) * k2) + ((- (12 * k2) - 0.3e1 * k4 - 0.3e1) * x * x +  (18 * k2) + 0.3e1 * k4 + 0.3e1) *  t - 0.2e1 * ( k2 + k4 / 0.2e1 + 0.1e1 / 0.2e1) * (x - 0.1e1) * (x + 0.1e1)) * r * r +  (1 + k2) * (-0.3e1 +   pow( t,  4) * x * x - 0.4e1 *  (t * t) * x * x + (0.9e1 * x * x - 0.1e1) *  t) * r -  (2 *  pow( t,  3)) - 0.2e1 *  t * x * x +  (4 * t)) *  (l + 1) *  k2 * pow( k2 * r * r + 0.1e1, -0.3e1) / (0.1e1 + r * r *  k2 *   pow( (-1 + t),  2) * x * x) *   pow( (t + 1),  (-2));

	L0_B = -((l * l * pow(t, 0.3e1) + (-l * l - 0.3e1 * l - 0.1e1) * t * t + (-0.3e1 * l * l - 0.2e1 * l) * t - l * l - l - 0.1e1) * x * x *  (1 + k2) * pow(-0.1e1 + t, 0.2e1) * k6 * pow(r, 0.7e1) - ((l * l - l) * pow(t, 0.4e1) + (-l * l * x * x - 0.3e1 * l * x * x - l * l - 0.2e1 * x * x + 0.3e1 * l) * t * t + l * l * x * x + l * x * x + 0.2e1) * k6 * pow(r, 0.6e1) + 0.3e1 * k4 *  (1 + k2) * (l * l * pow(t, 0.5e1) * x * x - 0.3e1 * x * x * (l * l + l + 0.1e1 / 0.3e1) * pow(t, 0.4e1) + (0.4e1 * l * x * x + 0.2e1 / 0.3e1 * l * l + 0.2e1 * x * x + 0.2e1 / 0.3e1 * l) * pow(t, 0.3e1) + (0.10e2 / 0.3e1 * l * l * x * x + (-0.5e1 / 0.3e1 * x * x - 0.1e1 / 0.3e1) * l - 0.3e1 * x * x - 0.1e1 / 0.3e1) * t * t + (-l * l * x * x + l * x * x / 0.3e1 - 0.4e1 / 0.3e1 * l * l + 0.7e1 / 0.3e1 * x * x - l + 0.1e1 / 0.3e1) * t - l * l * x * x / 0.3e1 - l * x * x / 0.3e1 - 0.2e1 / 0.3e1 * l * l - x * x - 0.2e1 / 0.3e1 * l) * pow(r, 0.5e1) - 0.3e1 * (l *  k2 * (l - 0.1e1) * pow(t, 0.4e1) + (((-0.5e1 / 0.3e1 * x * x - 0.1e1 / 0.3e1) *  k2 - (x - 0.1e1) * (x + 0.1e1) * (k4 + 0.1e1) / 0.3e1) * l * l - 0.13e2 / 0.3e1 * ( k2 + 0.2e1 / 0.13e2 * k4 + 0.2e1 / 0.13e2) * (x - 0.1e1) * (x + 0.1e1) * l + (-0.8e1 / 0.3e1 * x * x + 0.2e1 / 0.3e1) *  k2 - (x - 0.1e1) * (x + 0.1e1) * (k4 + 0.1e1) / 0.3e1) * t * t + 0.2e1 / 0.3e1 * (x - 0.1e1) * ( k2 + k4 / 0.2e1 + 0.1e1 / 0.2e1) * (l + 0.1e1) * (x + 0.1e1) * t + ((0.5e1 / 0.3e1 * x * x - 0.2e1 / 0.3e1) *  k2 + (x - 0.1e1) * (x + 0.1e1) * (k4 + 0.1e1) / 0.3e1) * l * l + ((0.5e1 / 0.3e1 * x * x - 0.2e1 / 0.3e1) *  k2 + (x - 0.1e1) * (x + 0.1e1) * (k4 + 0.1e1) / 0.3e1) * l +  (2 * k2)) *  k2 * pow(r, 0.4e1) + 0.3e1 *  (1 + k2) * (l * l * pow(t, 0.5e1) * x * x - 0.3e1 * x * x * (l * l + l + 0.1e1 / 0.3e1) * pow(t, 0.4e1) + (0.4e1 * l * x * x + 0.4e1 / 0.3e1 * l * l + 0.2e1 * x * x + 0.4e1 / 0.3e1 * l) * pow(t, 0.3e1) + 0.8e1 / 0.3e1 * x * x * (l * l - 0.3e1 / 0.2e1 * l - 0.7e1 / 0.4e1) * t * t + ((-x * x - 0.8e1 / 0.3e1) * l * l - 0.4e1 / 0.3e1 * l + 0.2e1 * x * x + 0.4e1 / 0.3e1) * t + l * l * x * x / 0.3e1 + l * x * x / 0.3e1 - 0.4e1 / 0.3e1 * l * l - x * x - 0.4e1 / 0.3e1 * l) *  k2 * pow(r, 0.3e1) + (-0.3e1 * l *  k2 * (l - 0.1e1) * pow(t, 0.4e1) + (((0.5e1 * x * x + 0.1e1) *  k2 + (k4 + 0.1e1) * x * x - k4 - 0.1e1) * l * l + 0.17e2 * ( k2 + 0.4e1 / 0.17e2 * k4 + 0.4e1 / 0.17e2) * (x - 0.1e1) * (x + 0.1e1) * l + (0.12e2 * x * x - 0.6e1) *  k2 + (0.3e1 * k4 + 0.3e1) * x * x - 0.3e1 * k4 - 0.3e1) * t * t + 0.2e1 * (x - 0.1e1) * ( k2 + k4 / 0.2e1 + 0.1e1 / 0.2e1) * (l + 0.1e1) * (x + 0.1e1) * t + ((-0.5e1 * x * x + 0.2e1) *  k2 + (-k4 - 0.1e1) * x * x + k4 + 0.1e1) * l * l + ((-0.5e1 * x * x + 0.2e1) *  k2 + (-k4 - 0.1e1) * x * x + k4 + 0.1e1) * l -  (6 * k2)) * r * r +  (1 + k2) * (l * l * pow(t, 0.5e1) * x * x - 0.3e1 * x * x * (l * l + l + 0.1e1 / 0.3e1) * pow(t, 0.4e1) + (0.4e1 * l * x * x + 0.2e1 * l * l + 0.2e1 * x * x + 0.2e1 * l) * pow(t, 0.3e1) + (0.2e1 * l * l * x * x - 0.7e1 * l * x * x - 0.7e1 * x * x + l + 0.1e1) * t * t + (-l * l * x * x - l * x * x - 0.4e1 * l * l + x * x - l + 0.3e1) * t + l * l * x * x + l * x * x - 0.2e1 * l * l - x * x - 0.2e1 * l) * r + (-l * l + l) * pow(t, 0.4e1) + (l * l * x * x + 0.3e1 * l * x * x + l * l + 0.2e1 * x * x - 0.3e1 * l) * t * t - l * l * x * x - l * x * x - 0.2e1) *  k2 * pow( k2 * r * r + 0.1e1, -0.3e1) / (0.1e1 + r * r *  k2 * pow(-0.1e1 + t, 0.2e1) * x * x) * pow(t + 0.1e1, -0.2e1);

	Source_f2 = S*df1_t + L1_0*df1 + L2_0*f1 + (L0_A*Plp1 + L0_B*Pl)*alpha;

	
	
	return Source_f2;

}