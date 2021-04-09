#include "2+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
void allocate_derivs(derivs *V, int nl, int nh)
{	
	(*V).X   = dvector(nl, nh);
	(*V).X1  = dvector(nl, nh);
	(*V).X11 = dvector(nl, nh);
	(*V).X2   = dvector(nl, nh);
	(*V).X12  = dvector(nl, nh);
	(*V).X22 = dvector(nl, nh);
	
	(*V).Y   = dvector(nl, nh);
	(*V).Y1  = dvector(nl, nh);
	(*V).Y11 = dvector(nl, nh);
	(*V).Y2   = dvector(nl, nh);
	(*V).Y12  = dvector(nl, nh);
	(*V).Y22 = dvector(nl, nh);
}
// -------------------------------------------------------------------------------
void fill0_derivs(derivs V, int nl, int nh)
{	
	fill0_dvector(V.X,   nl, nh);
	fill0_dvector(V.X1,  nl, nh);
	fill0_dvector(V.X11, nl, nh);
	fill0_dvector(V.X2,   nl, nh);
	fill0_dvector(V.X12,  nl, nh);
	fill0_dvector(V.X22, nl, nh);
	
	fill0_dvector(V.Y,   nl, nh);
	fill0_dvector(V.Y1,  nl, nh);
	fill0_dvector(V.Y11, nl, nh);
	fill0_dvector(V.Y2,   nl, nh);
	fill0_dvector(V.Y12,  nl, nh);
	fill0_dvector(V.Y22, nl, nh);
}
// -------------------------------------------------------------------------------
void copy_derivs(derivs Vout, derivs Vin, int nl, int nh)
{	
	copy_dvector(Vout.X,   Vin.X,   nl, nh);
	copy_dvector(Vout.X1,  Vin.X1,  nl, nh);
	copy_dvector(Vout.X11, Vin.X11, nl, nh);
	copy_dvector(Vout.X2,   Vin.X2,   nl, nh);
	copy_dvector(Vout.X12,  Vin.X12,  nl, nh);
	copy_dvector(Vout.X22, Vin.X22, nl, nh);
	
	copy_dvector(Vout.Y,   Vin.Y,   nl, nh);
	copy_dvector(Vout.Y1,  Vin.Y1,  nl, nh);
	copy_dvector(Vout.Y11, Vin.Y11, nl, nh);
	copy_dvector(Vout.Y2,   Vin.Y2,   nl, nh);
	copy_dvector(Vout.Y12,  Vin.Y12,  nl, nh);
	copy_dvector(Vout.Y22, Vin.Y22, nl, nh);
}
// -------------------------------------------------------------------------------
void free_derivs(derivs *V, int nl, int nh)
{	
	free_dvector((*V).X,   nl, nh);
	free_dvector((*V).X1,  nl, nh);
	free_dvector((*V).X11, nl, nh);
	free_dvector((*V).X2,   nl, nh);
	free_dvector((*V).X12,  nl, nh);
	free_dvector((*V).X22, nl, nh);	

	free_dvector((*V).Y,   nl, nh);
	free_dvector((*V).Y1,  nl, nh);
	free_dvector((*V).Y11, nl, nh);
	free_dvector((*V).Y2,   nl, nh);
	free_dvector((*V).Y12,  nl, nh);
	free_dvector((*V).Y22, nl, nh);
}
// -------------------------------------------------------------------------------
void get_Chebyshev_Coefficients_derivs(derivs psi, derivs c, int N, char *grid){
	Chebyshev_Coefficients(psi.X,   c.X,   N, grid);
	Chebyshev_Coefficients(psi.X1,  c.X1,  N, grid);
	Chebyshev_Coefficients(psi.X11, c.X11, N, grid);
	Chebyshev_Coefficients(psi.X2,  c.X2,  N, grid);
	Chebyshev_Coefficients(psi.X12, c.X12, N, grid);
	Chebyshev_Coefficients(psi.X22, c.X22, N, grid);


	Chebyshev_Coefficients(psi.Y,   c.Y,   N, grid);
	Chebyshev_Coefficients(psi.Y1,  c.Y1,  N, grid);
	Chebyshev_Coefficients(psi.Y11, c.Y11, N, grid);
	Chebyshev_Coefficients(psi.Y2,  c.Y2,  N, grid);
	Chebyshev_Coefficients(psi.Y12, c.Y12, N, grid);
	Chebyshev_Coefficients(psi.Y22, c.Y22, N, grid);
}
// -------------------------------------------------------------------------------
void allocate_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch)
{	
	(*V).X   = dmatrix(nrl, nrh, ncl, nch);
	(*V).X1  = dmatrix(nrl, nrh, ncl, nch);
	(*V).X11 = dmatrix(nrl, nrh, ncl, nch);
	(*V).X2  = dmatrix(nrl, nrh, ncl, nch);
	(*V).X12 = dmatrix(nrl, nrh, ncl, nch);
	(*V).X22 = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y   = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y1  = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y11 = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y2  = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y12 = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y22 = dmatrix(nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
void fill0_cheb_derivs(cheb_derivs V, int nrl, int nrh, int ncl, int nch)
{	
	fill0_dmatrix(V.X,   nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X1,  nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X11, nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X2,   nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X12,  nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X22, nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y,   nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y1,  nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y11, nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y2,   nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y12,  nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y22, nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
void copy_cheb_derivs(cheb_derivs Vout, cheb_derivs Vin, int nrl, int nrh, int ncl, int nch)
{	
	copy_dmatrix(Vout.X,   Vin.X,   nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X1,  Vin.X1,  nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X11, Vin.X11, nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X2,   Vin.X2,   nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X12,  Vin.X12,  nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X22, Vin.X22, nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y,   Vin.Y,   nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y1,  Vin.Y1,  nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y11, Vin.Y11, nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y2,   Vin.Y2,   nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y12,  Vin.Y12,  nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y22, Vin.Y22, nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
void free_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch)
{	
	free_dmatrix((*V).X,   nrl, nrh, ncl, nch);
	free_dmatrix((*V).X1,  nrl, nrh, ncl, nch);
	free_dmatrix((*V).X11, nrl, nrh, ncl, nch);
	free_dmatrix((*V).X2,   nrl, nrh, ncl, nch);
	free_dmatrix((*V).X12,  nrl, nrh, ncl, nch);
	free_dmatrix((*V).X22, nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y,   nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y1,  nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y11, nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y2,   nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y12,  nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y22, nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// void allocate_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch)
// {	
// 	(*V).X   = dmatrix(nrl, nrh, ncl, nch);
// 	(*V).X1  = dmatrix(nrl, nrh, ncl, nch);
// 	(*V).X11 = dmatrix(nrl, nrh, ncl, nch);
// 	(*V).Y   = dmatrix(nrl, nrh, ncl, nch);
// 	(*V).Y1  = dmatrix(nrl, nrh, ncl, nch);
// 	(*V).Y11 = dmatrix(nrl, nrh, ncl, nch);
// }
// // -------------------------------------------------------------------------------
// void fill0_cheb_derivs(cheb_derivs V, int nrl, int nrh, int ncl, int nch)
// {	
// 	fill0_dmatrix(V.X,   nrl, nrh, ncl, nch);
// 	fill0_dmatrix(V.X1,  nrl, nrh, ncl, nch);
// 	fill0_dmatrix(V.X11, nrl, nrh, ncl, nch);
// 	fill0_dmatrix(V.Y,   nrl, nrh, ncl, nch);
// 	fill0_dmatrix(V.Y1,  nrl, nrh, ncl, nch);
// 	fill0_dmatrix(V.Y11, nrl, nrh, ncl, nch);
// }
// // -------------------------------------------------------------------------------
// void copy_cheb_derivs(cheb_derivs Vout, cheb_derivs Vin, int nrl, int nrh, int ncl, int nch)
// {	
// 	copy_dmatrix(Vout.X,   Vin.X,   nrl, nrh, ncl, nch);
// 	copy_dmatrix(Vout.X1,  Vin.X1,  nrl, nrh, ncl, nch);
// 	copy_dmatrix(Vout.X11, Vin.X11, nrl, nrh, ncl, nch);
// 	copy_dmatrix(Vout.Y,   Vin.Y,   nrl, nrh, ncl, nch);
// 	copy_dmatrix(Vout.Y1,  Vin.Y1,  nrl, nrh, ncl, nch);
// 	copy_dmatrix(Vout.Y11, Vin.Y11, nrl, nrh, ncl, nch);
// }
// // -------------------------------------------------------------------------------
// void free_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch)
// {	
// 	free_dmatrix((*V).X,   nrl, nrh, ncl, nch);
// 	free_dmatrix((*V).X1,  nrl, nrh, ncl, nch);
// 	free_dmatrix((*V).X11, nrl, nrh, ncl, nch);
// 	free_dmatrix((*V).Y,   nrl, nrh, ncl, nch);
// 	free_dmatrix((*V).Y1,  nrl, nrh, ncl, nch);
// 	free_dmatrix((*V).Y11, nrl, nrh, ncl, nch);
// }
// // -------------------------------------------------------------------------------
