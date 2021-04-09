#include "2+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
void Butcher_Tableau(parameters *par)
{    
	int s = (*par).SDIRK_s;
	double L = 0.4358665215084589994, L2 = L*L;

	(*par).SDIRK_A = dmatrix(1, s, 1, s);
	(*par).SDIRK_B = dvector(1, s);
	(*par).SDIRK_C = dvector(1, s);

	if (s==3){
		(*par).SDIRK_A[1][1] = L;            (*par).SDIRK_A[1][2]= 0;      (*par).SDIRK_A[1][3]= 0;
		(*par).SDIRK_A[2][1] = 0.5*(1-L);    (*par).SDIRK_A[2][2]= L;      (*par).SDIRK_A[2][3]= 0;
		(*par).SDIRK_A[3][1] = 0.25*(-6*L2 + 16*L-1.);
		(*par).SDIRK_A[3][2] = 0.25*(6*L2-20*L+5.);
		(*par).SDIRK_A[3][3] = L;
		

		
		(*par).SDIRK_B[1] = 0.25*(-6*L2+16*L-1.);
		(*par).SDIRK_B[2] = 0.25*(6*L2-20*L+5.);
		(*par).SDIRK_B[3] = L;
		
		(*par).SDIRK_C[1] = L;
		(*par).SDIRK_C[2] = 0.5 * ( 1.+L);
		(*par).SDIRK_C[3] = 1.;
		
		(*par).gamma = L;
	}
	else{
		(*par).SDIRK_A[1][1] = 1;
		(*par).SDIRK_B[1]    = 1;
		(*par).SDIRK_C[1]    = 1;
		(*par).gamma         = 1;
	}
}
// -------------------------------------------------------------------------------
void free_Butcher_Tableau(parameters *par)
{    
	int s = (*par).SDIRK_s;

	free_dmatrix((*par).SDIRK_A, 1, s, 1, s);
	free_dvector((*par).SDIRK_B, 1, s);
	free_dvector((*par).SDIRK_C, 1, s);
}
// -------------------------------------------------------------------------------
