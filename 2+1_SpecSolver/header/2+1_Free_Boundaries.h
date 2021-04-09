#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "utilities.h"
#include "spectral_utilities.h"

#define nDoms     1
#define nFields   2
#define nFields1D 2

#define Newton_SDIRK_itmin  1      
#define Newton_SDIRK_itmax  5      
#define Newton_SDIRK_tol    1.e-11 
#define Newton_SDIRK_verb   0

#define IG_bicgstab_findiff_decr     1.e-3   
#define IG_bicgstab_findiff_verb     0
#define IG_bicgstab_findiff_itmax    200


#define Newton_itmin      1    
#define Newton_itmax      10     
#define Newton_tol        5.e-10
#define Newton_verb       1

#define bicgstab_decr     1.e-10 
#define bicgstab_verb     1
#define bicgstab_itmax    100

#define Newton_LinSolve_SDIRK_itmin  1      
#define Newton_LinSolve_SDIRK_itmax  1 
#define Newton_LinSolve_SDIRK_tol 5.e-5
#define Newton_LinSolve_SDIRK_verb 0

#define SDIRK_bicgstab_findiff_decr     5.e-6
#define SDIRK_bicgstab_findiff_verb     0
#define SDIRK_bicgstab_findiff_itmax    200

#define debug_verb 0

#define FD_ORDER 4
#if FD_ORDER == 2
 #define STENCILSIZE 9 
#elif FD_ORDER == 4
 #define STENCILSIZE 25
#elif FD_ORDER == 6
 #define STENCILSIZE 49
#endif


//Parameters must coincide with function output in Mathematica
//Also grid type must coincide
#define N0_read 25
#define N1_read 25
#define N2_read 4
#define Nkappa_read 25
#define r1_read 0.1

typedef struct DERIVS{ 
	double *X, *X1, *X11, *X2, *X12, *X22,  
		   *Y, *Y1, *Y11, *Y2, *Y12, *Y22;
} derivs;

typedef struct CHEB_DERIVS{ 
	double **X, **X1, **X11, **X2, **X12, **X22, 
		   **Y, **Y1, **Y11, **Y2, **Y12, **Y22;
} cheb_derivs;

typedef struct PARAMETERS{ 	
	int N0, N1, N2, n0, n1, n2,
		nslice, nslice_total, ngrid , ntotal, n_count,
		i0_max, i1_max, i2_max, SDIRK_s, SOLVER_METHOD;
	
	double 	*grid_points_x0, *grid_points_x1, *grid_points_x2,
			angle_incr_x0, angle_incr_x1, angle_incr_x2, 
			**D_x0, **D_x1, **D_x2, 
			bicgstab_tol, IG_bicgstab_findiff_tol, SDIRK_bicgstab_findiff_tol,
			**SDIRK_A, *SDIRK_B, *SDIRK_C, gamma, gah;
		
	char	SimName[200], grid_x0[50], grid_x1[50],  grid_x2[50];
	
	derivs ID;
	
	
	double t_initial, t_final, t0, t1, r0, r1, x0, x1;
	
	//********* Problem specific parameters ********
	int ell, EXACT_SOL_CHECK, EXACT_SOL_SOLVE;
	double alpha;
	double kappa;
	//**********************************************
	
} parameters;

typedef struct JFD_COMPONENTS{ 
	double	**J_UL, **K_UL, **Ku_UL, **Kl_UL;
	int	*ncols_J_UL, **cols_J_UL, *iK_UL,
		m1, m2;
		
	double	**J_UR, **J_LL, **J_LR, **K_UR, **K_LR, **Klu_LR;
	int    *indx_LR;	
} JFD_Components;


// Subroutines in "2+1_Free_Boundaries.c"
int main();

//Subroutines in allocate_and_free_derivs.c
void allocate_derivs(derivs *V, int nl, int nh);
void fill0_derivs(derivs V, int nl, int nh);
void copy_derivs(derivs Vout, derivs Vin, int nl, int nh);
void free_derivs(derivs *V, int nl, int nh);
void get_Chebyshev_Coefficients_derivs(derivs psi, derivs c, int N, char *grid);
void allocate_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch);
void fill0_cheb_derivs(cheb_derivs V, int nrl, int nrh, int ncl, int nch);
void copy_cheb_derivs(cheb_derivs Vout, cheb_derivs Vin, int nrl, int nrh, int ncl, int nch);
void free_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch);

// Subroutines in "butcher_tableau.c"
void Butcher_Tableau(parameters *par);
void free_Butcher_Tableau(parameters *par);


// Subroutines in "coordinates.c"
double func_t_of_x0(parameters par, double x0);
double func_r_of_x1(parameters par, double x1);
double func_x_of_x2(parameters par, double x2);
void scale_derivs_x0x1x2_to_trx(parameters par, derivs W_x0x1x2, derivs W1D_x0, derivs W_trx, derivs W1D_t, int nsize, int nsize_1D);
int total_grid_points(char *grid, int N);
double get_grid_angle_incremente(char *grid, int N);
double get_grid_point(char *grid, int i, int N);
void Get_grid_Arrays(parameters *par);
void Free_grid_Arrays(parameters *par);

// Subroutines in "derivatives.c"
void Get_Differentiation_Matrices(parameters *par);
void Free_Differentiation_Matrices(parameters *par);
void Derivatives_Slice(parameters par, derivs V);
void Derivatives(parameters par, derivs V);
void Get_DerivativesFinDif(parameters par, derivs Vs);
void Get_DerivativesFinDif_grid(parameters par, int i_field, int i1, int i2, derivs V);

// Subroutines in "debug.c"
void pause();
void PrintVector(FILE *stream, double *J, int na, int nb, double threashold);
void PrintMatrix(FILE *stream, double **J, int la, int lb, int ca, int cb);

// Subroutines in "field_and_bound_eqs.c"
void F_of_XY(parameters par, double x0, double x1, double x2, derivs W, derivs W1D, double *F);
void FBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, double *F);
void JDX_of_XY(parameters par, double x0, double x1, double x2, derivs W, derivs W1D, derivs DW, derivs DW1D, double *JDX);
void JDXBound_of_XY(parameters par, double x0, derivs Wslice, derivs W1D, derivs DWslice, derivs DW1D, double *JDX);


// Subroutines in "func_and_jacobian.c"
void get_W3D_at_Grid(parameters par, int i0, int i1, int i2, derivs ID, derivs V, derivs W);
void get_W3D_at_slice(parameters par, int i0, int i1, int i2, derivs ID, derivs V, derivs Wslice);
void get_W1D_at_Grid(parameters par, int i0, derivs ID, derivs V, derivs W1D);
void F_of_X(parameters par, double *X, double *F);
void J_times_DX(parameters par, double *X, double *DX, double *JDX);

// Subroutines in "index_routines.c
int Index_Boundary(char grid[], int i, int n);
int Index_Slice(parameters par, int iField, int j1, int j2);
int Index_Slice_Fields1D(int iField, int nslice);
int Index(int iField, int i0, int N0, int i1, int N1, int i2, int N2);
int Index_Fields1D(int iField, int ngrid, int i0, int N0);
void get_indices_from_Index_Slice(parameters par, int indx, int *iField, int *i1, int *i2);

// Subroutines in "initial_data.c"
void get_Initial_Data(parameters *par);
void get_Initial_Data_newTimeStep(parameters *par, derivs W);


// Subroutines in "initial_guess.c"
void get_initial_guess(parameters par, double *X);
void SDIRK_Initial_Guess(parameters par, double *X);
void SDIRK_Step_IG(parameters par, double x0, double h, double *X_previous, double *X_next);
int Newton_SDIRK_IG(parameters par, double x0, double *Mj, double *kj);
void get_Ws_from_M_and_k(parameters par, double *Mj, double *kj, derivs Wslice, derivs W1D);
void get_W2D_at_Grid(parameters par, int i1, int i2, derivs Wslice, derivs W);
void F_of_kj_IG(parameters par, double x0, double *Mj, double *kj, /*derivs Wslice, derivs W1D, */double *F);
void J_times_Dkj_IG(parameters par, double x0,  double *Mj, double *kj, double *Dkj, /*derivs Wslice, derivs W1D,*/ double *JDqj);
void Jacobian_SDIRK_IG(parameters par, double x0, double *Mj, double *kj,/*derivs Wslice, derivs W1D,*/ double **J);
void DF_of_X_SDIRK_IG(parameters par, double x0, double *Mj, double *kj, /*derivs Wslice, derivs W1D,*/ double **J);

// Subroutines in "initial_guess_finn_diff_precondirioner.c"
int Newton_SDIRK_IG_bicgstab(parameters par, double x0, double *Mj, double *kj);
void get_BandMatrix(parameters par, double x0, double *Mj, double *kj, JFD_Components *JFD);
void free_bandMatrix(parameters par, JFD_Components *JFD);
void Get_JFD_Matrix(parameters par, double x0, double *Mj, double *kj, JFD_Components *JFD);
void Get_JFD_Components(parameters par, JFD_Components *JFD);
void PreCond_FinDiff(parameters par, JFD_Components JFD, double *b, double *Y, double *res);
int bicgstab_findiff(parameters par, double x0, double *Mj, double *kj, double *Dkj, double *F, double *normres);


// Subroutines in "newton_direct.c"
void Jacobian(parameters par, double *X, double **J);
void DF_of_X(parameters par, double *X, double **J);
int newton_direct(parameters par, double *X);

// Subroutines in "newton_bicgstab.c"
int newton(parameters par, double *X);
void get_chebV_x0(parameters par, derivs V, cheb_derivs chebV);
int bicgstab(parameters par, double *X, double *F, double *DX, double *normres);

// Subroutines in "output.c"
void output_solution_timeframes(parameters par, derivs Sol);
void output_solution_fields1D(parameters par, derivs Sol);
void output_legendre_modes(parameters par, derivs Sol);

// Subroutines in "parameteres.c"
void set_parameters(int nmax, parameters *par);
void free_parameters(parameters *par);


// Subroutines in "problem_specific_functions.c"
void Func_Delta(double kappa, double R, double *Delta, double *Delta_R);
void Func_Sigma(double kappa, double R, double x, double *Sigma, double *Sigma_R, double *Sigma_x);
void Func_Sigma0(double kappa, double R, double *Sigma0, double *Sigma0_R);
void Func_R(double t, double r, double *R, double *R_r, double *R_t);
void Func_C(double kappa, double r, double *C, double *C_r);
void Load_ExactSolutions(parameters par, double ***c_f2[], double ***c_df2[]);
void ExactSolution_interpolate(parameters par, double x0, double x1, double x2, double ***c_f2[], double ***c_df2[], double *f1, double *f2, double *df1, double *df2);
void get_ExactSolution_ID(parameters *par);
void get_ExactSolution(parameters par, double *X);
double get_Source_f1(parameters par, double t);
double get_Source_f2(parameters par, double t, double r, double x,  derivs W1D_t );
void get_LinearEquationCoefficients(double kappa, double t, double r, double x, 
	double *S, double *L1_0, double *L1_r, double *L1_x,
	double *L2_0, double *L2_r, double *L2_x, double *L2_rr, double *L2_rx, double *L2_xx);

// Subroutines in "SDIRK_preconditioner.c"
int idx_order_time_flow(int j, int N0);
void PreCond(parameters par, cheb_derivs chebV, double *F, double *X);
void get_chebF(parameters par, double *F, double **chebF);
void SDIRK_Step(parameters par, cheb_derivs chebV, double **chebF, double x0, double h, double *Xn, double *Xnp1);
void Jacobian_LinSolve_SDIRK(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double  **J);
int Newton_LinSolve_SDIRK(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj);
void LinSolve_SDIRK(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj);
void Get_Wslice_W1D_Fslice(parameters par, double x0, double **chebF, cheb_derivs chebV, derivs Wslice, derivs W1D, double *Fslice);
void F_of_kj(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F);

// Subroutines in "SDIRK_preconditioner_fin_diff_preconditioner.c"
int LinSolve_SDIRK_bicgstab(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj);
int SDIRK_bicgstab_findiff(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F, double *normres);
void get_SDIRK_BandMatrix(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F, JFD_Components *JFD);
void Get_SDIRK_JFD_Matrix(parameters par, cheb_derivs chebV, double **chebF, double x0, double *Mj, double *kj, double *F, JFD_Components *JFD);


// Subroutines in solve_equations.c
void Get_V_From_X(parameters par, double *X, derivs V);
void get_solution_from_X(parameters par, double *X, derivs W);
void solve_equations(parameters par, double *X, derivs Sol);







