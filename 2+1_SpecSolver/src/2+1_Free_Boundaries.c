#include "2+1_Free_Boundaries.h"

//   -----------------------------------------------
double func(double t, double r, double x){

 return plgndr(1, 0, x)*exp(r*cos(2*Pi*r*t))*plgndr(3, 0, t);

}


int main() 
{
	parameters par;
	double *X;
	derivs Sol;

	printf("\n------------------------------------------------\n");fflush(0);	
	printf("   SCALAR FIELD ON KERR AROUND I^0-CYLINDER\n");fflush(0);	
	printf("------------------------------------------------\n\n");fflush(0);	
	
	
	int  nmax=1;

	clock_t start, end;
	double cpu_time_used_total;
	
	
	for(par.n_count=0; par.n_count<nmax; par.n_count++){
		
		set_parameters(nmax, &par);	

	
		X = dvector(0, par.ntotal-1);
		allocate_derivs(&Sol, 0, par.ntotal - 1);
		

		if(par.n_count==0) 
			get_Initial_Data(&par);
		else{
			get_Initial_Data_newTimeStep(&par, Sol);
			free_derivs(&Sol, 0, par.ntotal - 1);
		}


				cpu_time_used_total = 0;
				start = clock();
		get_initial_guess(par, X);
		
				end = clock();
				cpu_time_used_total = ((double) (end - start)) / CLOCKS_PER_SEC;
				printf("Total time = %lf\n", cpu_time_used_total);
		
				start = clock();
		solve_equations(par, X, Sol);
				end = clock();
				cpu_time_used_total = ((double) (end - start)) / CLOCKS_PER_SEC;
				printf("\nTotal time = %lf\n", cpu_time_used_total);


		
		// // //output_solution_timeframes(par, Sol);
		output_solution_fields1D(par, Sol);
		output_legendre_modes(par, Sol);
		
		free_dvector(X,  0, par.ntotal-1);
		free_parameters(&par);
	}
	free_derivs(&Sol, 0, par.ntotal - 1);


	printf("\a\n");
	return 1;
}
