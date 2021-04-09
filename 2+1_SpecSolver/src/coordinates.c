#include "2+1_Free_Boundaries.h"


double func_t_of_x0(parameters par, double x0){
 
 return par.t0 + 0.5*(par.t1-par.t0)*(1+x0);
}
//-----------------------------------------------
double func_r_of_x1(parameters par, double x1){
 
 return par.r0 + 0.5*(par.r1-par.r0)*(1+x1);
}
//-----------------------------------------------
double func_x_of_x2(parameters par, double x2){
 
 return par.x0 + 0.5*(par.x1-par.x0)*(1+x2);
}
//   -----------------------------------------------
void scale_derivs_x0x1x2_to_trx(parameters par, derivs W_x0x1x2, derivs W1D_x0, derivs W_trx, derivs W1D_t, int nsize, int nsize_1D){
  double t0=par.t0, t1=par.t1, fact = 2./(t1-t0),
	 	 r0=par.r0, r1=par.r1, facr = 2./(r1-r0),
	 	 x0=par.x0, x1=par.x1, facx = 2./(x1-x0);
	 
  
 int i;

 for(i=0; i<nsize; i++ ){
   		W_trx.X[i] = W_x0x1x2.X[i];
		W_trx.X1[i] = facr*W_x0x1x2.X1[i];
		W_trx.X11[i] = sqr(facr)*W_x0x1x2.X11[i];
		W_trx.X2[i] = facx*W_x0x1x2.X2[i];
		W_trx.X12[i] = facx*facr*W_x0x1x2.X12[i]; 
		W_trx.X22[i] = sqr(facx)*W_x0x1x2.X22[i]; 
		
		W_trx.Y[i] = fact*W_x0x1x2.Y[i];
		W_trx.Y1[i] = fact*facr*W_x0x1x2.Y1[i];
		W_trx.Y11[i] = fact*sqr(facr)*W_x0x1x2.Y11[i]; 
		W_trx.Y2[i] = fact*facx*W_x0x1x2.Y2[i];
		W_trx.Y12[i] = fact*facx*facr*W_x0x1x2.Y12[i]; 
		W_trx.Y22[i] = fact*sqr(facx)*W_x0x1x2.Y22[i];
 }
 
  for(i=0; i<nsize_1D; i++ ){
   		W1D_t.X[i] = W1D_x0.X[i];
		W1D_t.Y[i] = fact*W1D_x0.Y[i];
 }
 
  
}
//   -----------------------------------------------------------
int total_grid_points(char *grid, int N){
  int n;
  
  	  if(strcmp( grid,"Fourier" ) ==0)
	    n=2*N+1;
	  else
	    n=N+1;
	  
  
  return n;
}
//-----------------------------------------------
double get_grid_angle_incremente(char *grid, int N){
 	
	double arg=0;
	
	if(strcmp( grid,"Radau_RHS" ) ==0)
	  arg = 2*Pi/(2*N+1);	  
	else if(strcmp(grid,"Radau_LHS")==0)
	  arg = - 2*Pi/(2*N+1);
	else if(strcmp(grid,"Gauss")==0)
	  arg = Pi/(N+1);
	else if(strcmp( grid,"Lobatto")==0)
	  arg = Pi/N;
	else if(strcmp( grid,"Fourier")==0){
	  arg = 2.*Pi/(2*N+1);
	}
	else{ 
	  fprintf(stderr,
	  "Error in get_grid_angle_incremente: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto / Fourier\n input grid was: %s\n", grid);
	  exit(1);
	}
	
	return arg;
}
//-----------------------------------------------
double get_grid_point(char *grid, int i, int N){
 	
	double arg=0, x=0;
	
	if(strcmp( grid,"Radau_RHS" ) ==0){
	  arg = 2*Pi*i/(2*N+1);
	  x = cos(arg); 
	}
	else if(strcmp(grid,"Radau_LHS")==0){
	  arg = Pi - 2*Pi*i/(2*N+1);
	  x = cos(arg);
	}
	else if(strcmp(grid,"Gauss")==0){
	  arg = Pi*(i+0.5)/(N+1);
	  x = cos(arg);
	}
	else if(strcmp( grid,"Lobatto")==0){
	  arg = Pi*i/N;
	  x = cos(arg);
	}
	else if(strcmp( grid,"Fourier")==0){
	  x = 2.*Pi*i/(2*N+1);
	}
	else{ 
	  fprintf(stderr,
	  "Error in get_grid_point: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto / Fourier\n input grid was: %s\n", grid);
	  exit(1);
	}
	
	
	
	return x;
}
// -------------------------------------------------------------------------------
void Get_grid_Arrays(parameters *par)
{
	int i, N0=(*par).N0, N1=(*par).N1, N2=(*par).N2;
	
	(*par).grid_points_x0  = dvector(0, N0);
	(*par).grid_points_x1  = dvector(0, N1);
	(*par).grid_points_x2  = dvector(0, N2);
	
	(*par).angle_incr_x0 = get_grid_angle_incremente((*par).grid_x0, N0);
	for(i=0; i<=N0; i++){
			double x = get_grid_point((*par).grid_x0, i, N0);
			(*par).grid_points_x0[i]=x;		
	}
	(*par).angle_incr_x1 = get_grid_angle_incremente((*par).grid_x1, N1);
	for(i=0; i<=N1; i++){
			double x = get_grid_point((*par).grid_x1, i, N1);
			(*par).grid_points_x1[i]=x;		
	}
	(*par).angle_incr_x2 = get_grid_angle_incremente((*par).grid_x2, N2);
	for(i=0; i<=N2; i++){
			double x = get_grid_point((*par).grid_x2, i, N2);
			(*par).grid_points_x2[i]=x;		
	}

return;
}
// -------------------------------------------------------------------------------
void Free_grid_Arrays(parameters *par)
{
	int N0=(*par).N0, N1=(*par).N1, N2=(*par).N2;
	free_dvector((*par).grid_points_x0, 0, N0);
	free_dvector((*par).grid_points_x1, 0, N1);
	free_dvector((*par).grid_points_x2, 0, N2);

}