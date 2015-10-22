#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP SeqHD_REAL2(SEXP Rdata, SEXP dupe)
{
   int col, i, k, N, M, pos;
   SEXP Rdimdata;
   double *data, buffer=0;
   Rdimdata = getAttrib(Rdata, R_DimSymbol);
   N = INTEGER(Rdimdata)[0];
   M = INTEGER(Rdimdata)[1];
   SEXP Rdata_copy;
   
   if(INTEGER(dupe)[0]==1)
   {
     Rdata_copy = PROTECT(duplicate(Rdata));
     data = REAL(Rdata_copy);
     //printf(" %d \n", INTEGER(dupe)[0]);
   }
   else
   {
     //printf(" %d \n", INTEGER(dupe)[0]);
     data = REAL(Rdata);
   }
      
   for(k=0;k<M;k++)
   {
      col=N*k;
      for(i=0;i<N;i++)
      {
         pos = i+col;
         if(ISNA(data[pos]))
         {data[pos]= buffer;}
         else{buffer=data[pos];}
      }
   }
   /*
    for(k=0;k<M;k++)
   {
      col=N*k;
      for(i=0;i<N;i++)
      {
         pos = i+col;
         printf("%f \t", data[pos]);
      }
      printf("\n");
   }
   */
   if(INTEGER(dupe)[0]==1)
   {
     //printf("return %d \n", INTEGER(dupe)[0]);
     UNPROTECT(1);
     return Rdata_copy;
   }
   else
   {
     //printf("return %d \n", INTEGER(dupe)[0]);
     return Rdata;
   }
}

