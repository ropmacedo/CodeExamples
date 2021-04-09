#include "2+1_Free_Boundaries.h"


void output_solution_timeframes(parameters par, derivs Sol){


	int iFields, i0, i1, i2, N0=par.N0, N1=par.N1, N2=par.N2, 
		i0_max=par.i0_max, i1_max=par.i1_max, i2_max=par.i2_max;
	
	double ***f=d3tensor(0,N0, 0,N1, 0,N2), ***cf=d3tensor(0,N0, 0,N1, 0,N2);
	
	FILE *fp;
	char fn[250];
	
	for(iFields=0; iFields<nFields; iFields++){
		for(i0=0; i0<=N0; i0++){
			for(i1=0; i1<=N1; i1++){
				for(i2=0; i2<=N2; i2++){				
					int Idx = Index(iFields, i0, N0, i1, N1, i2, N2) ;
					f[i0][i1][i2] = Sol.X[Idx];
				}
			}
		}
		
		Chebyshev_Coefficients_3D(f, cf, N0, par.grid_x0, N1, par.grid_x1, N2, par.grid_x2);
		
		
		for(i0=0; i0<=par.i0_max; i0++){
			double x0=-1.+2.*i0/i0_max, t=func_t_of_x0(par, x0);
			
			sprintf(fn, "data/%s/Solution_field_%d_timeframe_%04d.dat", par.SimName, iFields, i0);


			fp=fopen(fn, "w");
			fprintf(fp, "#1:t\t 2:r\t 3:x\t 4:Solution\n");
			
			for(i1=0; i1<=par.i1_max; i1++){
				double x1=-1.+2.*i1/i1_max, r=func_r_of_x1(par, x1);
				
				for(i2=0; i2<=par.i2_max; i2++){
					double x2=-1.+2.*i2/i2_max, x=func_x_of_x2(par, x2);
					
					double output = Clenshaw_Chebyshev_3D(cf, N0, N1, N2, x0, x1, x2);
					
					fprintf(fp, "%3.20e\t%3.20e\t%3.20e\t%3.20e\n", t, r, x, output);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
		}	
		
	}


	free_d3tensor(f,0,N0, 0,N1, 0,N2);
	free_d3tensor(cf, 0,N0, 0,N1, 0,N2);
	
return;
}
//----------------------------------------------------------
void output_solution_fields1D(parameters par, derivs Sol){
	
	int iFields, i0, N0=par.N0, i0_max=par.i0_max;
	
	FILE *fp;
	char fn[250];
	
	double *f=dvector(0,N0), *cf=dvector(0,N0);
	
	for(iFields=0; iFields<nFields; iFields++){
		for(i0=0; i0<=N0; i0++){
			int Idx = Index_Fields1D(iFields, par.ngrid, i0, N0);
			
			f[i0] = Sol.X[Idx];
		}
		Chebyshev_Coefficients(f, cf, N0, par.grid_x0);
		
		sprintf(fn, "data/%s/ChebCoeff_field1D_%d.dat", par.SimName, iFields);
		fp=fopen(fn, "w");
		fprintf(fp, "#1:i0\t 2:c_i0\n");
		for(i0=0; i0<=N0; i0++) fprintf(fp, "%d\t %3.20e\n", i0, cf[i0]);
		fclose(fp);
		
		sprintf(fn, "data/%s/Solution_field1D_%d.dat", par.SimName, iFields);
		fp=fopen(fn, "w");
		fprintf(fp, "#1:t\t 2:Solution\n");
		
		for(i0=0; i0<=par.i0_max; i0++){
			double x0=-1.+2.*i0/i0_max, t=func_t_of_x0(par, x0);
			
			double output = Clenshaw_Chebyshev(cf, N0, x0);
					
			fprintf(fp, "%3.20e\t%3.20e\n", t, output);
		}
		fclose(fp);
		
		
		
	}
	
	
	free_dvector(f,0,N0);
	free_dvector(cf, 0,N0);
	
	
return;
}
//----------------------------------------------------------
void output_legendre_modes(parameters par, derivs Sol){

	int iField, i0, i1, i2, ii2, N0=par.N0, N1=par.N1, N2=par.N2, 
		i0_max=par.i0_max, i1_max=par.i1_max;
	
	
	double **H=dmatrix(0,N2, 0, N2);
	for(i2=0;i2<=N2;i2++){
		for(ii2=0;ii2<=N2;ii2++){
		H[i2][ii2]= Chebyshev_Definite_Integral_Collocations_Matrix(i2, ii2, N2);
		}		
	}
	
	int l_max = N2, l1;
	double **Pl = dmatrix(0,l_max, 0, N2);

	
	for(l1=0; l1<=l_max; l1++){
		for(i2=0;i2<=N2;i2++){
	 		double x2=par.grid_points_x2[i2];
	 		Pl[l1][i2] = plgndr(l1, 0, x2);
		}
	}
	
	double ***f=d3tensor(0, l_max, 0, N0 , 0 , N1), ***cf=d3tensor(0, l_max, 0, N0 , 0 , N1),
		   ***f0=d3tensor(0, l_max, 0, N1 , 0 , N0), ***cf1=d3tensor(0, l_max, 0, N0 , 0 , N1),
		   ***cf0=d3tensor(0, l_max, 0, N1 , 0 , N0),
		   ***output=d3tensor(0, l_max, 0, i0_max , 0 , i1_max);
	
	FILE *fp, *fp_cheb, *fp_scri;
	char fn[250], fn_cheb[250], fn_scri[250];

	double *f_scri=dvector(0,N1), *cf_scri=dvector(0,N1);
	
	for(iField=0; iField<nFields; iField++){

		sprintf(fn, "data/%s/LegendreSolution3D_field_%d.dat", par.SimName, iField);
		fp=fopen(fn, "w");
		fprintf(fp, "#1:t\t 2:r\t ");
		for(l1=0; l1<=l_max; l1++){
			fprintf(fp, "%d:Solution_ell=%d\t ", 3+l1, l1);

			

		
			for(i0=0; i0<=N0; i0++){
				for(i1=0; i1<=N1; i1++){
					f[l1][i0][i1]=0.;
					for(i2=0; i2<=N2; i2++){
						int Idx = Index(iField, i0, N0, i1, N1, i2, N2) ;
						
						for(ii2=0; ii2<=N2; ii2++){					
							f[l1][i0][i1]+=0.5*(2.*l1+1)*Pl[l1][ii2]*H[i2][ii2]*Sol.X[Idx];
							
						}											
					}
					f0[l1][i1][i0]=f[l1][i0][i1];
				}
				
				Chebyshev_Coefficients(f[l1][i0], cf1[l1][i0], N1, par.grid_x1);			
			}
			
			for(i1=0; i1<=N1; i1++){
				Chebyshev_Coefficients(f0[l1][i1], cf0[l1][i1], N0, par.grid_x0);			
			}		
			
			
			Chebyshev_Coefficients_2D(f[l1], cf[l1], N0, par.grid_x0, N1, par.grid_x1);
			
			for(i0=0; i0<=par.i0_max; i0++){
				double x0=-1.+2.*i0/i0_max;				
					for(i1=0; i1<=par.i1_max; i1++){
						double x1=-1.+2.*i1/i1_max;
						output[l1][i0][i1] = Clenshaw_Chebyshev_2D(cf[l1], N0, N1, x0, x1);
					}
				}

			for(i1=0;i1<=N1; i1++){
				double x0 = 1., x1 = par.grid_points_x1[i1];
				f_scri[i1]=Clenshaw_Chebyshev_2D(cf[l1], N0, N1, x0, x1);
			}
			Chebyshev_Coefficients(f_scri, cf_scri, N1, par.grid_x1);
			
			sprintf(fn_scri, "data/%s/ChebScri_field_%d_ell_%d.dat", par.SimName, iField, l1);
			fp_scri=fopen(fn_scri, "w");

			for(i1=0;i1<=N1; i1++)
				fprintf(fp_scri, "%d\t%3.15e\n", i1, cf_scri[i1]);

			fclose(fp_scri);

		}
		fprintf(fp, "\n");
		
		
		for(i0=0;i0<=i0_max; i0++){
			double x0=-1.+2.*i0/i0_max, t=func_t_of_x0(par, x0);
			for(i1=0;i1<=i1_max; i1++){
				double x1=-1.+2.*i1/i1_max, r=func_r_of_x1(par, x1);
				
				fprintf(fp, "%3.20e\t %3.20e\t ", t , r);
				for(l1=0; l1<=l_max; l1++){
					fprintf(fp, "%3.20e\t ", output[l1][i0][i1]);
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		
		
		for(l1=0; l1<=l_max; l1++){
			sprintf(fn, "data/%s/LegendreSolution2D_field_%d_ell_%d.dat", par.SimName, iField, l1);
			fp=fopen(fn, "w");
			
			sprintf(fn_cheb, "data/%s/Cheb_x0_LegendreSolution_field_%d_ell_%d.dat", par.SimName, iField, l1);
			fp_cheb=fopen(fn_cheb, "w");
			
			fprintf(fp_cheb, "#1:i0\t ");
			for(i1=N1; i1>=0; i1--){
				double r=func_r_of_x1(par, par.grid_points_x1[i1]);
				fprintf(fp_cheb, "%d:c_i0(r=%lf)\t ", 2+(N1-i1), r);
			}
			fprintf(fp_cheb, "\n");
			
			
			
			for(i0=0; i0<=N0; i0++){
				fprintf(fp_cheb, "%d\t ", i0);
				for(i1=N1; i1>=0; i1--){
					fprintf(fp_cheb, "%3.20e\t ",cf0[l1][i1][i0]);
				}
				fprintf(fp_cheb, "\n ");
			}
			fclose(fp_cheb);
			
			sprintf(fn_cheb, "data/%s/Cheb_x1_LegendreSolution_field_%d_ell_%d.dat", par.SimName, iField, l1);
			fp_cheb=fopen(fn_cheb, "w");
			fprintf(fp_cheb, "#1:i1\t ");
			for(i0=N0; i0>=0; i0--){
				double t=func_t_of_x0(par, par.grid_points_x0[i0]);
				fprintf(fp_cheb, "%d:c_i1(t=%lf)\t ", 2+(N0-i0), t);
			}
			fprintf(fp_cheb, "\n");
			
			for(i1=0; i1<=N1; i1++){
				fprintf(fp_cheb, "%d\t ", i1);
				for(i0=N0; i0>=0; i0--){
					fprintf(fp_cheb, "%3.20e\t ",cf1[l1][i0][i1]);
				}
				fprintf(fp_cheb, "\n ");
			}
			fclose(fp_cheb);
			
			fprintf(fp, "#1:t\t ");
			for(i1=0;i1<=i1_max; i1++){
				double x1=-1.+2.*i1/i1_max, r=func_r_of_x1(par, x1);
				fprintf(fp, "%d:Solution(r=%lf)\t ", 2+i1, r);
				
			}
			fprintf(fp, "\n");
			
			for(i0=0;i0<=i0_max; i0++){
				double x0=-1.+2.*i0/i0_max, t=func_t_of_x0(par, x0);
				fprintf(fp, "%3.20e\t ", t);
				
				for(i1=0;i1<=i1_max; i1++){
					fprintf(fp, "%3.20e\t ", output[l1][i0][i1]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);

			sprintf(fn, "data/%s/LegendreSolution2D_radialDep_field_%d_ell_%d.dat", par.SimName, iField, l1);
			fp=fopen(fn, "w");
			fprintf(fp, "#1:r\t ");
			for(i0=0;i0<=i0_max; i0++){
				double x0=-1.+2.*i0/i0_max, t=func_t_of_x0(par, x0);
				fprintf(fp, "%d:Solution(t=%lf)\t ", 2+i0, t);
				
			}
			fprintf(fp, "\n");
			
			for(i1=0;i1<=i1_max; i1++){
				double x1=-1.+2.*i1/i1_max, r=func_r_of_x1(par, x1);
				fprintf(fp, "%3.20e\t ", r);
				
				for(i0=0;i0<=i0_max; i0++){
					fprintf(fp, "%3.20e\t ", output[l1][i0][i1]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);

		}						
	}
	
	
	
	free_d3tensor( f, 0, l_max, 0, N0 , 0 , N1);
	free_d3tensor( f0, 0, l_max, 0, N1 , 0 , N0);
	free_d3tensor( cf0, 0, l_max, 0, N1 , 0 , N0);
	free_d3tensor( cf1, 0, l_max, 0, N0 , 0 , N1);
	free_d3tensor(cf, 0, l_max, 0, N0 , 0 , N1);
	free_d3tensor(output, 0, l_max, 0, i0_max , 0 , i1_max );
	
	free_dvector( f_scri, 0, N1);
	free_dvector(cf_scri, 0, N1);

	free_dmatrix(H, 0, N2, 0, N2);
	free_dmatrix(Pl, 0, l_max, 0, N2);

	return;	
}
