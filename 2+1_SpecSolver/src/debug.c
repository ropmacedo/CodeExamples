#include "2+1_Free_Boundaries.h"

void pause(){
	char pause;
	printf("\nPress any key to continue\n(or 'q' to quit)\n");
	// __fpurge(stdin);
   	fpurge(stdin);
	int ret;
	ret = scanf("%c", &pause);
	if (ret == EOF)
	{
	  printf("...Error in 'Pause'\n");
	  exit(1);
	}
	if(pause=='q')
		exit(1);
	// __fpurge(stdin);
	fpurge(stdin);
	printf("\n");
}
//------------------------------------------------------------------------------
void PrintVector(FILE *stream, double *J, int na, int nb, double threashold){
 int i;
 
 for(i=na; i<=nb; i++){
 if(fabs(J[i])>threashold)
   fprintf(stream, "vec[%d]=%3.15e\n",i, J[i]);
//    pause();
 }
}
//------------------------------------------------------------------------------
void PrintMatrix(FILE *stream, double **J, int la, int lb, int ca, int cb){
 int i, j;
 
 for(i=la; i<=lb; i++){
   for(j=ca; j<=cb; j++){
   	double out = /*fabs(J[i][j])<5.e-9? 0.:*/J[i][j];
     fprintf(stream, "%3.2e ", out);
   }
   fprintf(stream, "\n");
 }
}
//------------------------------------------------------------------------------
