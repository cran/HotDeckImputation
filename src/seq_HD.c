#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
/*
SEXP SeqHD_REAL(SEXP Rdata)
{
   int col, i, k, N, M, pos;
   SEXP Rdimdata;
   double *data, buffer=0;
   Rdimdata = getAttrib(Rdata, R_DimSymbol);
   N = INTEGER(Rdimdata)[0];
   M = INTEGER(Rdimdata)[1];
   
   data = REAL(Rdata);
   
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
   return R_NilValue;
}
*/

SEXP SeqHD_REAL(SEXP Rdata, SEXP RInitialvalue, SEXP RNAValues)
{
  int col, i, k, N, M, pos;
  SEXP Rdimdata;
  double *data, *NAValues, *initval;
  Rdimdata = getAttrib(Rdata, R_DimSymbol);
  N = INTEGER(Rdimdata)[0];
  M = INTEGER(Rdimdata)[1];
  initval = REAL(RInitialvalue);
  NAValues = REAL(RNAValues);
  
  int *buffer_vals = malloc(M * sizeof(int));
  
  for(k=0;k<M;k++)
  {
    buffer_vals[k]=initval[k];
  }
  
  data = REAL(Rdata);
  
  for(k=0;k<M;k++)
  {
    col=N*k;
    if(ISNA(NAValues[k]))
    {
      for(i=0;i<N;i++)
      {
        pos = i+col;
        //printf("data %f, NA %f \n", data[pos],NAValues[k]);
        if(ISNA(data[pos]))
        {
          //printf("equal \n");
          data[pos]= buffer_vals[k];
        }
        else{
          //printf("not equal \n");
          buffer_vals[k]=data[pos];}
      }
    }
    else
    {
      for(i=0;i<N;i++)
      {
        pos = i+col;
        //printf("data %f, NA %f \n", data[pos],NAValues[k]);
        if(data[pos]==NAValues[k])
        {
          //printf("equal \n");
          data[pos]= buffer_vals[k];
        }
        else{
          //printf("not equal \n");
          buffer_vals[k]=data[pos];}
      }
    }
  }
  
  free(buffer_vals);
  return R_NilValue;
}

SEXP SeqHD_INTEGER(SEXP Rdata, SEXP RInitialvalue, SEXP RNAValues)
{
   int col, i, k, N, M, pos;
   SEXP Rdimdata;
   int *data, *NAValues, *initval;
   Rdimdata = getAttrib(Rdata, R_DimSymbol);
   N = INTEGER(Rdimdata)[0];
   M = INTEGER(Rdimdata)[1];
   initval = INTEGER(RInitialvalue);
   NAValues = INTEGER(RNAValues);
   
   int *buffer_vals = malloc(M * sizeof(int));
   
   for(k=0;k<M;k++)
   {
      buffer_vals[k]=initval[k];
   }
   
   data = INTEGER(Rdata);
   
   for(k=0;k<M;k++)
   {
      col=N*k;
      for(i=0;i<N;i++)
      {
         pos = i+col;
         if(data[pos]==NAValues[k])
         {data[pos]= buffer_vals[k];}
         else{buffer_vals[k]=data[pos];}
      }
   }
   
   free(buffer_vals);
   return R_NilValue;
}

SEXP CPS_SeqHD_INTEGER(SEXP Rdata, SEXP Rcovariates, SEXP RImputeables, SEXP RInitialvalue, SEXP RNAValues)
{
   int n , ncov, nimp, ncells, ncell_nequal;
   ncov = LENGTH(Rcovariates);
   nimp = LENGTH(RImputeables);
   n = INTEGER(getAttrib(Rdata, R_DimSymbol))[0];
   
   int i, covcounter, impcounter, pos_cov, pos_imp, cellcounter, pos_cell;
   int *data, *initval, *cov, *imp, *NAValues;
   data = INTEGER(Rdata);
   cov = INTEGER(Rcovariates);
   imp = INTEGER(RImputeables);
   initval = INTEGER(RInitialvalue);
   NAValues = INTEGER(RNAValues);
   
   //use first object encountered to set buffer and class values
   i=0;
   ncells = 1;
   int *adjustment_cells =   malloc(ncells* ncov * sizeof(int));
   int *buffer_vals = malloc(ncells * nimp * sizeof(int));
   
   for(covcounter = 0; covcounter< ncov; covcounter++)
   {
      pos_cov = i + cov[covcounter]*n;
      adjustment_cells[covcounter + ncov*(ncells-1)] = data[pos_cov];
   }
 
   for(impcounter = 0; impcounter< nimp; impcounter++)
   {
      pos_imp = i + imp[impcounter]*n;
      if(data[pos_imp] == NA_INTEGER)
      {
         data[pos_imp] = initval[impcounter];
      }
      buffer_vals[impcounter + nimp*(ncells-1)] = data[pos_imp];
   }

   
   // 1st class has been found, buffervals have been set for this class
   // if the first object had missing values, these have now been imputed
   for(i=1;i<n;i++)
   {
      ncell_nequal = 0;
      //check, if this object already falls into any of the created adjustment cells
      for(cellcounter =0; cellcounter < ncells; cellcounter++)
      {
         for(covcounter = 0; covcounter< ncov; covcounter++)
         {
            pos_cov = i + cov[covcounter]*n; // this is for data
            pos_cell = covcounter + ncov*(cellcounter); //this is for the adjustment cell
            if(adjustment_cells[pos_cell] != data[pos_cov])
            {
               ncell_nequal++;
               break;
            }
         }
         if(covcounter == ncov) // no break encountered; all covariates equal, this object belongs into this cell
         {
            for(impcounter = 0; impcounter< nimp; impcounter++)
            {
               pos_imp = i + imp[impcounter]*n;
               if(data[pos_imp] == NAValues[impcounter])
               {
                  data[pos_imp] = buffer_vals[impcounter + nimp*(cellcounter)];
               }
               else
               {
                  buffer_vals[impcounter + nimp*(cellcounter)] = data[pos_imp];
               }
            }
            break;
         }
      }
      if(ncell_nequal == ncells) // no cell found for this object, new cell must be created and initialized
      {
         ncells++;
         adjustment_cells = realloc (adjustment_cells, ncells* ncov * sizeof(int));
         buffer_vals = realloc (buffer_vals, ncells * nimp * sizeof(int));
         
         //set the adjustment cell values
         for(covcounter = 0; covcounter< ncov; covcounter++)
         {
            pos_cov = i + cov[covcounter]*n;
            adjustment_cells[covcounter + ncov*(ncells-1)] = data[pos_cov];
         }
         
         //set the buffer Values 
         for(impcounter = 0; impcounter< nimp; impcounter++)
         {
            pos_imp = i + imp[impcounter]*n;
            if(data[pos_imp] == NAValues[impcounter])
            {
               data[pos_imp] = initval[impcounter];
            }
            buffer_vals[impcounter + nimp*(ncells-1)] = data[pos_imp];
         }
      }
   }
   
   free(adjustment_cells);
   free(buffer_vals);
   return R_NilValue;
}


SEXP deepcopy_c(SEXP Rdata)
{
  SEXP Rdata_copy = duplicate(Rdata);
  return Rdata_copy;
}
