#include "2+1_Free_Boundaries.h"

//-----------------------------------------------
void set_parameters(int nmax, parameters *par) 
{
	
	//----------------------------------------------------
	(*par).t_initial = 0.;
	(*par).t_final = 0.9;
	
	(*par).r0 = 0;
	(*par).r1 = 0.1;
	
	(*par).x0 = -1.;
	(*par).x1 = 1.;


	//SPIN PARAMETER
	(*par).kappa = 0.5;
	
	//ANGULAR MODE FOR ID
	(*par).ell = 2;

	//CONTRIBUTION FROM TAYLOR MODE 0 AT CYLINDER
	(*par).alpha = 1.;


	//SOLVE EQUATIONS WITH ID FROM EXACT SOLUTION IF 1
	(*par).EXACT_SOL_SOLVE = 0;

	//CHECK EQUATIONS AGAINST EXACT SOLUTION IF 1
	(*par).EXACT_SOL_CHECK  = 0;


	if((*par).EXACT_SOL_CHECK == 1 || (*par).EXACT_SOL_SOLVE == 1){
		//CONTRIBUTION FROM TAYLOR MODE 0 AT CYLINDER
		//NOTE: for exact solution with l=0, alpha = 1-kappa^2; 
		//      otherwise alpha = 0.
		//		Also, after normalisation f2Norm_l0 = (  2*log((1-r1*k^2)/(1-r1)) - r1*(1-k^2)( 2+r1*(1+k2)  )  )/(2*r1^3)
		double k=(*par).kappa, r1=(*par).r1, k2=sqr(k);
		(*par).alpha = (*par).ell == 0? (1.-k2)*2*pow(r1,3.)/(  2*log((1-r1*k2)/(1-r1)) -r1*(1-k2)*(2+r1*(1+k2)) ):0.;
	}
	
		
	
	// sprintf((*par).grid_x0,"Radau_RHS");
	sprintf((*par).grid_x0,"Gauss");
	sprintf((*par).grid_x1,"Lobatto");
	sprintf((*par).grid_x2,"Lobatto");
	
	(*par).SOLVER_METHOD = 2;
	//Spectral Resolution
	(*par).N0 = 20;
	(*par).N1 = 15;
	(*par).N2 = 4;
	
	//Output resolution
	(*par).i0_max = 100;
	(*par).i1_max = 100;
	(*par).i2_max = 20;
	
	
	sprintf((*par).SimName, "SDIRK_kappa%lf_l%d_N0_%d_N1_%d_N2_%d_tf%f",(*par).kappa, (*par).ell, (*par).N0,(*par).N1, (*par).N2,(*par).t_final );// Simulation Name
	//sprintf((*par).SimName, "test" );// Simulation Name

	
	//----------------------------------------------------

	double DT = ( (*par).t_final - (*par).t_initial) /nmax;
	(*par).t0 = (*par).t_initial + (*par).n_count*DT;
	(*par).t1 = (*par).t_initial + ((*par).n_count+1)*DT;
	
	int n0, n1, n2;
	(*par).n0 = n0 = total_grid_points((*par).grid_x0, (*par).N0);
	(*par).n1 = n1 = total_grid_points((*par).grid_x1, (*par).N1);
	(*par).n2 = n2 = total_grid_points((*par).grid_x2, (*par).N2);
	

	(*par).nslice =  nFields*n1*n2;
	(*par).nslice_total =  (*par).nslice + nFields1D;	
	
	(*par).ngrid = (*par).nslice*n0;
	(*par).ntotal =  (*par).nslice_total*n0;
	
	
	char mkdir[200];
	int func_out;	
	sprintf(mkdir, "mkdir -p data/%s", (*par).SimName);	
	func_out =  system(mkdir);

	sprintf((*par).SimName, "%s/%04d", (*par).SimName, (*par).n_count);
	sprintf(mkdir, "mkdir -p data/%s", (*par).SimName);	
	func_out =  system(mkdir);
	
	
	Get_grid_Arrays(par);
	Get_Differentiation_Matrices(par);
	
	allocate_derivs( &((*par).ID), 0, (*par).nslice_total - 1);
		
	(*par).SDIRK_s= 3;
	Butcher_Tableau(par);
		
		
	
	return;
}


//-----------------------------------------------
void free_parameters(parameters *par){

	Free_grid_Arrays(par);
	Free_Differentiation_Matrices(par);
	free_derivs(&((*par).ID), 0, (*par).nslice_total-1);
	free_Butcher_Tableau(par);

}
