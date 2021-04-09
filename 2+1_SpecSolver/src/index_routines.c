#include "2+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
int Index_Boundary(char grid[], int i, int n)
{
	int N = n -1, i1 = i;
	if(i1 < 0){
		if(strcmp( grid,"Fourier" ) ==0){
		  i1 = i1 + n;		  
		}
		else if(strcmp( grid,"Lobatto")==0){
		  i1 = -i1;		  
		}
		else if(strcmp( grid,"Gauss")==0){
		  i1 = -(i1+1);		  
		}
		else if(strcmp( grid,"Radau_RHS" ) ==0){
		  i1 = -i1;		  
		}
		else if(strcmp( grid,"Radau_LHS" ) ==0){
		  i1 = -i1;		  
		}
		else{
		  printf("Error in Index_Boundary: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto/ Fourier\n grid was: %s\n", grid);
		  exit(1);
		}

	}
	else if(i1>N){
		if(strcmp( grid,"Fourier" ) ==0){
			i1 = i1 - n;
		}
		else if(strcmp( grid,"Lobatto")==0){
			i1 = 2*N-i1;
		}
		else if(strcmp( grid,"Gauss")==0){
		  	i1 = 2*N+1-i1;
		}
		else if(strcmp( grid,"Radau_RHS" ) ==0){
		 	i1 = 2*N+1-i1;
		}
		else if(strcmp( grid,"Radau_LHS" ) ==0){
		  	i1 = 2*N+1-i1;
		}
		else{
		  printf("Error in Index_Boundary: argument has to be: Radau_RHS / Radau_LHS / Gauss/ Lobatto/ Fourier\n grid was: %s\n", grid);
		  exit(1);
		}

	}
	else{
		i1=i;
	}

	return i1;
}
// -------------------------------------------------------------------------------
int Index_Slice(parameters par, int iField, int j1, int j2)
{		
	int i1, i2, indx, n1=par.n1, n2=par.n2;

	i1=Index_Boundary(par.grid_x1, j1, n1);
	i2=Index_Boundary(par.grid_x2, j2, n2);
	
	if (n1<=n2) indx = i1 + i2*n1;
	else  indx = i2 + i1*n2;
	
	return  iField + nFields*indx;
}

// -------------------------------------------------------------------------------
int Index_Slice_Fields1D(int iField, int nslice) 
{
	
	return nslice + iField;
}
// -------------------------------------------------------------------------------
int Index(int iField, int i0, int N0, int i1, int N1, int i2, int N2) 
{
	
	int n1=N1+1, n2=N2+1, n12=n1*n2, indx_t;
	
	if (n1<=n2) indx_t = i1 + i2*n1;
	else  indx_t = i2 + i1*n2;

	return iField + nFields*(indx_t + i0*n12);
}
// -------------------------------------------------------------------------------
int Index_Fields1D(int iField, int ngrid, int i0, int N0) 
{
	

	return ngrid + iField + nFields1D*i0;
}
//-------------------------------------------------------------------------------------------
void get_indices_from_Index_Slice(parameters par, int indx, int *iField, int *i1, int *i2)
{
  int n1=par.n1, n2=par.n2;
  
  	if(n1> n2){
  		int na,nb, n=n2;
  		*i1    = indx/(nFields*n);
		na    = indx - nFields*n*(*i1);
		*i2    = na/(nFields);
		nb = na - nFields*(*i2);
		*iField = nb;
	}
	else{
		int na,nb, n=n1;
		*i2    = indx/(nFields*n);
		na    = indx - nFields*n*(*i2);
		*i1    = na/(nFields);
		nb = na - nFields*(*i1);
		*iField = nb;
     }      	
}