/*
   Do not delete!
   File name		VAM.c
	Part of:		   HotDeckImputation (GNU R contributed package)
	Author:			Dieter William Joenssen
	Copyright:		Dieter William Joenssen
	Email:			Dieter.Joenssen@TU-Ilmenau.de
	Created:		   07 July 2013
	Last Update: 	10 August 2013
	Description:	C functions to perform Vogel's approximation method of matching
*/

//Begin preamble
#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void c_modifVAM(int *distribution, double *distances, int *supply, int *demand, int *nrow_s, int *ncol_d);

void c_set_value(double *vector, double value, int *vector_length);

void c_set_value_int(int *vector, int value, int *vector_length);

void c_VAM(int *distribution, double *distances, int *supply, int *demand, int *nrow_s, int *ncol_d);

//End preamble

//Begin code

void c_modifVAM(int *distribution, double *distances, int *supply, int *demand, int *nrow_s, int *ncol_d)
{
   int j=-1,zerocounter_demand=0,sum_demand = 0,current_pos=0;
   int col_offset=0;
   int i=-1,zerocounter_supply=0,sum_supply =0;
   double maxdelta_cols=-1;
   int maxdelta_cols_index=-1;
   int remaining_col=-1;
   int remaining_row=-1;
   double col_min,row_min;
   int index_col_min,index_row_min;
   
   //initialize column minimums and their position to zero
   double *col_mins1 = malloc(ncol_d[0] * sizeof(double));
   int *col_mins1_index = malloc(ncol_d[0] * sizeof(int));
   double *col_mins2 = malloc(ncol_d[0] * sizeof(double));
   int *col_mins2_index = malloc(ncol_d[0] * sizeof(int));
   double *col_delta = malloc(ncol_d[0] * sizeof(double));
   
   //test if there are any zero supplys or demand
   //proper initalization of the 
   for(i= 0; i< nrow_s[0]; i++)
   {
      if(supply[i]==0)
      {zerocounter_supply++;continue;}
      sum_supply= sum_supply+supply[i];
   }
   for(j=0; j< ncol_d[0];j++)
   {
      if(demand[j]==0)
      {zerocounter_demand++;continue;}
      sum_demand=sum_demand+demand[j];
   }
   //set column minimums to DBL_MAX
   c_set_value(col_mins1, DBL_MAX, ncol_d);
   c_set_value(col_mins2, DBL_MAX, ncol_d);
   //set delta to negative DBL_MAX
   c_set_value(col_delta,(-1*DBL_MAX), ncol_d);
   maxdelta_cols= (-1*DBL_MAX);
   
   while(((sum_demand!=0)&(sum_supply!=0)&(zerocounter_supply< (nrow_s[0]-1))&(zerocounter_demand< (ncol_d[0]-1))))
   {
      //get the column minimums of 1st and 2nd order
      for(j=0; j< ncol_d[0];j++)
      {
         if(demand[j]==0)
         {continue;}
         
         col_offset=j * nrow_s[0];
         
         for(i= 0; i< nrow_s[0]; i++)
         {
            if(supply[i]==0)
            {continue;}
            
            current_pos = i + col_offset;
            
            //if we find a new minimum, then the set the vectors for minimum positions and values
            //to this new minimum, the old minimum and minimum positions are then the minimum of order 2
            if(distances[current_pos]<=col_mins1[j])
            {
               col_mins2[j]=col_mins1[j];
               col_mins2_index[j]=col_mins1_index[j];
               col_mins1[j]=distances[current_pos];
               col_mins1_index[j]=i;
               continue;
            }
            //if not smaller than minimum, it might be smaller than the minimum of order 2
            //check and replace
            if(distances[current_pos]<=col_mins2[j])
            {
               col_mins2[j]=distances[current_pos];
               col_mins2_index[j]=i;
               continue;
            }
         }
         //with the found minimums, calculate their difference
         col_delta[j]=col_mins2[j]-col_mins1[j];
         //find the position of the maximum delta
         if(col_delta[j]>=maxdelta_cols)
         {
            maxdelta_cols=col_delta[j];
            maxdelta_cols_index= j;
         }
      }
      
      //we looked through all columns and found our maximum delta
      //we saved the column minumum for that colum so we can directly allocate
      //the maximum number of supply or demand
      col_offset=(maxdelta_cols_index * nrow_s[0]);
      current_pos = col_mins1_index[maxdelta_cols_index] + col_offset;
      
      //if supply == demand for this row col combination
      if(supply[col_mins1_index[maxdelta_cols_index]]==demand[maxdelta_cols_index])
      {
         distribution[current_pos]=demand[maxdelta_cols_index];
         
         sum_demand= sum_demand-demand[maxdelta_cols_index];
         zerocounter_demand++;
                  
         sum_supply= sum_supply- demand[maxdelta_cols_index];
         zerocounter_supply++;
         supply[col_mins1_index[maxdelta_cols_index]]=0;
         demand[maxdelta_cols_index]=0;
      }
      else
      {
         //if supply > demand allocate all demand, that demand is now zero with supply remaining
         if(supply[col_mins1_index[maxdelta_cols_index]]>demand[maxdelta_cols_index])
         {
            distribution[current_pos]=demand[maxdelta_cols_index];
            
            sum_demand= sum_demand-demand[maxdelta_cols_index];
            zerocounter_demand++;
                        
            sum_supply= sum_supply- demand[maxdelta_cols_index];
            supply[col_mins1_index[maxdelta_cols_index]]=supply[col_mins1_index[maxdelta_cols_index]]-demand[maxdelta_cols_index];
            demand[maxdelta_cols_index]=0;
         }
         else
         {
            //this shouldn't happen in this application, but just in case
            //if supply < demand allocate all supply, that supply is now zero with demand remaining
            distribution[current_pos]=supply[col_mins1_index[maxdelta_cols_index]];
         
            sum_demand= sum_demand-supply[col_mins1_index[maxdelta_cols_index]];
            demand[maxdelta_cols_index]=demand[maxdelta_cols_index]-supply[col_mins1_index[maxdelta_cols_index]];
            
            sum_supply= sum_supply- demand[maxdelta_cols_index];
            zerocounter_supply++;
            supply[col_mins1_index[maxdelta_cols_index]]=0;
         }
      }
   
      //set column minimums to DBL_MAX
      c_set_value(col_mins1, DBL_MAX, ncol_d);
      c_set_value(col_mins2, DBL_MAX, ncol_d);
      //set delta to negative DBL_MAX
      c_set_value(col_delta,(-1*DBL_MAX), ncol_d);
      maxdelta_cols= (-1*DBL_MAX);
      
   }
   
   //apparently either the supply or demand has beenexhausted, or we have encountered that:
   //only one column is left, one row is left or both
   //if supply is all done or demand is, then do nothing
   //printf("sum_supply = %d\n",sum_supply);
   //printf("sum_demand = %d\n",sum_demand);
   //printf("zerocounter_demand = %d\n",zerocounter_demand);
   //printf("zerocounter_supply = %d\n",zerocounter_supply);
   
   if((sum_demand==0)|(sum_supply==0))
   {
      //do nothing
   }
   else
   {
      //if only one col is left
      if(zerocounter_demand == (ncol_d[0]-1))
      {
         
         //find the col where only value is left
         for(j=0; j< ncol_d[0];j++)
         {
            if(demand[j]!=0)
            {
               remaining_col = j;
               col_offset = j* nrow_s[0];
               break;
            }
         }
         //printf("remaining_col = %d\n",remaining_col);
         
         //check if only one row is left
         if(zerocounter_supply== (nrow_s[0]-1))
         {
            //find that row
            for(i=0; i< nrow_s[0];i++)
            {
               if(supply[i]!=0)
               {remaining_row = i;}
            }
            //printf("remaining_row = %d\n",remaining_row);
            
            current_pos=remaining_row+col_offset;
            //now allocate everything to that position and you are done
            if(demand[remaining_col]==supply[remaining_row])
            {
               distribution[current_pos]=demand[remaining_col];
         
               sum_demand= sum_demand-demand[remaining_col];
               zerocounter_demand++;
                        
               sum_supply= sum_supply- demand[remaining_col];
               zerocounter_supply++;
               supply[remaining_row]=0;
               demand[remaining_col]=0;
            }
            else
            {
               if(demand[remaining_col]<supply[remaining_row])
               {
                  distribution[current_pos]=demand[remaining_col];
                  sum_demand= sum_demand-demand[remaining_col];
                  zerocounter_demand++;
                           
                  sum_supply= sum_supply- demand[remaining_col];
                  supply[remaining_row]=supply[remaining_row]-demand[remaining_col];
                  demand[remaining_col]=0;
               }
               else
               {
                  //demand > supply
                  distribution[current_pos]=supply[remaining_row];
         
                  sum_demand= sum_demand-supply[remaining_row];
                  demand[remaining_col]=demand[remaining_col]-supply[remaining_row];
                  
                  sum_supply= sum_supply- supply[remaining_row];
                  zerocounter_supply++;
                  supply[remaining_row]=0;
               }
            }
         }
         else
         {
            //itteratively find the minimum value and allocate till done
            //while(zerocounter_demand!=ncol_d[0])
            while((sum_demand!=0)&(sum_supply!=0)&(zerocounter_supply< (nrow_s[0]-1))&(zerocounter_demand< (ncol_d[0]-1)))
            {
               col_min=DBL_MAX;
               index_col_min=-1;
               //find row with minimum
               for(i=0;i<nrow_s[0];i++)
               {
                  if(supply[i]==0)
                  {continue;}
                  
                  current_pos=i+col_offset;
                  //test if value is minimum
                  if(distances[current_pos]<=col_min)
                  {
                     col_min=distances[current_pos];
                     index_col_min = i;
                  }
               }
               
               current_pos=index_col_min+col_offset;
               //allocate supply and demand
               if(demand[remaining_col]==supply[index_col_min])
               {
                  distribution[current_pos]=demand[remaining_col];
            
                  sum_demand= sum_demand-demand[remaining_col];
                  zerocounter_demand++;
                           
                  sum_supply= sum_supply- demand[remaining_col];
                  zerocounter_supply++;
                  supply[index_col_min]=0;
                  demand[remaining_col]=0;
               }
               else
               {
                  if(demand[remaining_col]<supply[index_col_min])
                  {
                     distribution[current_pos]=demand[remaining_col];
                     sum_demand= sum_demand-demand[remaining_col];
                     zerocounter_demand++;
                              
                     sum_supply= sum_supply- demand[remaining_col];
                     supply[index_col_min]=supply[index_col_min]-demand[remaining_col];
                     demand[remaining_col]=0;
                  }
                  else
                  {
                     //demand > supply
                     distribution[current_pos]=supply[index_col_min];
            
                     sum_demand= sum_demand-supply[index_col_min];
                     demand[remaining_col]=demand[remaining_col]-supply[index_col_min];
                     
                     sum_supply= sum_supply- supply[index_col_min];
                     zerocounter_supply++;
                     supply[index_col_min]=0;
                  }
               }
            }
            
         }
         
      }
      else
      {
         //only one row is left
         //souldn't happen in this application,but just to be shure
         //find the row where only value is left
         for(i=0; i< nrow_s[0];i++)
         {
            if(supply[i]!=0)
            {remaining_row = i;break;}
         }
         //allocate the remaining supply to the demand itteratevly
         while(((sum_demand!=0)&(sum_supply!=0)&(zerocounter_supply< (nrow_s[0]))&(zerocounter_demand< (ncol_d[0]))))
         {
            row_min=DBL_MAX;
            index_row_min=-1;
            //find col with minimum
            for(j=0;j<ncol_d[0];j++)
            {
               if(demand[j]==0)
               {continue;}
               col_offset=j * nrow_s[0];
               current_pos=remaining_row+col_offset;
               //test if value is minimum
               if(distances[current_pos]<=row_min)
               {
                  row_min=distances[current_pos];
                  index_row_min = j;
               }
            }
            
            col_offset=index_row_min * nrow_s[0];
            current_pos=remaining_row+col_offset;
            //allocate supply and demand
            if(demand[index_row_min]==supply[remaining_row])
            {
               distribution[current_pos]=demand[index_row_min];
         
               sum_demand= sum_demand-demand[index_row_min];
               zerocounter_demand++;
                        
               sum_supply= sum_supply- demand[index_row_min];
               zerocounter_supply++;
               supply[remaining_row]=0;
               demand[index_row_min]=0;
            }
            else
            {
               if(demand[index_row_min]<supply[remaining_row])
               {
                  distribution[current_pos]=demand[index_row_min];
                  sum_demand= sum_demand-demand[index_row_min];
                  zerocounter_demand++;
                           
                  sum_supply= sum_supply- demand[index_row_min];
                  supply[remaining_row]=supply[remaining_row]-demand[index_row_min];
                  demand[index_row_min]=0;
               }
               else
               {
                  //demand > supply
                  distribution[current_pos]=supply[remaining_row];
         
                  sum_demand= sum_demand-supply[remaining_row];
                  demand[index_row_min]=demand[index_row_min]-supply[remaining_row];
                  
                  sum_supply= sum_supply- supply[remaining_row];
                  zerocounter_supply++;
                  supply[remaining_row]=0;
               }
            }
         }
         
      }
   }

   free(col_mins1);
   free(col_mins1_index);
   free(col_mins2);
   free(col_mins2_index);
   free(col_delta);
}

void c_set_value(double *vector, double value, int *vector_length)
{
   int i;
   for( i = 0; i < vector_length[0]; i++)
	{
      vector[i]=value;
	}
}

void c_set_value_int(int *vector, int value, int *vector_length)
{
   int i;
   for( i = 0; i < vector_length[0]; i++)
   {
      vector[i]=value;
	}
}

void c_VAM(int *distribution, double *distances, int *supply, int *demand, int *nrow_s, int *ncol_d)
{
   
   //initialize column minimums
   double *col_mins1 = malloc(ncol_d[0] * sizeof(double));
   int *col_mins1_index = malloc(ncol_d[0] * sizeof(int));
   double *col_mins2 = malloc(ncol_d[0] * sizeof(double));
   //int *col_mins2_index = malloc(ncol_d[0] * sizeof(int));
   //double *col_delta = malloc(ncol_d[0] * sizeof(double));
   
   //initialize row minimums 
   double *row_mins1 = malloc(nrow_s[0] * sizeof(double));
   int *row_mins1_index = malloc(nrow_s[0] * sizeof(int));
   double *row_mins2 = malloc(nrow_s[0] * sizeof(double));
   //int *row_mins2_index = malloc(nrow_s[0] * sizeof(int));
   //double *row_delta = malloc(nrow_s[0] * sizeof(double));
   
   //test if there are any zero supplys or demand
   int zerocounter_demand=0,sum_demand=0;
   int zerocounter_supply=0,sum_supply=0;
   int i=0,j=0, current_pos=0, col_offset=0;
   double max_delta,current_delta;
   int min_position_max_delta_row=-1,min_position_max_delta_col=-1;
   
   for(j=0; j< ncol_d[0];j++)
   {
      if(demand[j]==0)
      {zerocounter_demand++;continue;}
      sum_demand=sum_demand+demand[j];
   }
      
   for(i= 0; i< nrow_s[0]; i++)
   {
      if(supply[i]==0)
      {zerocounter_supply++;continue;}
      sum_supply= sum_supply+supply[i];
   }
   
   if(sum_supply == sum_demand)
   {
   //while demand and is not zero and there is more than one row and column left (if only one row or colum left, extra rules needbe applyed)
   while(((sum_demand!=0)&(sum_supply!=0)&(zerocounter_supply< (nrow_s[0]-1))&(zerocounter_demand< (ncol_d[0]-1))))
   {
      //initialize decision variables
         //set values to dbl_max -1 or -1*dbl_max
            c_set_value(col_mins1,DBL_MAX, ncol_d);
            c_set_value_int(col_mins1_index,-1, ncol_d);
            c_set_value(col_mins2,DBL_MAX, ncol_d);
            //c_set_value_int(col_mins2_index,-1, ncol_d);
            //c_set_value(col_delta,(-1*DBL_MAX), ncol_d);
            
            c_set_value(row_mins1,DBL_MAX, ncol_d);
            c_set_value_int(row_mins1_index,-1, ncol_d);
            c_set_value(row_mins2,DBL_MAX, ncol_d);
            //c_set_value_int(row_mins2_index,-1, ncol_d);
            //c_set_value(row_delta,(-1*DBL_MAX), ncol_d);
      
         max_delta = 0;
         min_position_max_delta_row = -1;
         min_position_max_delta_col = -1;
      ////////
      //get column and row minimums of 1st and 2nd order
      //going col one then col two etc
      for(j=0; j < ncol_d[0]; j++)
      {
         //if demand==0 skip this column
         if(demand[j]==0)
         {continue;}
         
         //calculate which is the current offset for the column#
         col_offset=j * nrow_s[0];
         
         //go through the rows in the current column (so we are going down and then over)
         for(i=0;i<nrow_s[0]; i++)
         {
            //if supply ==0, skip this row
            if(supply[i]==0)
            {continue;}
            
            //current position in the "matrix"
            current_pos = i + col_offset;
            
            //evaluate for row mins (i) first
            //if we find a new minimum, then the set the vectors for minimum positions and values
            //to this new minimum, the old minimum and minimum positions are then the minimum of order 2
            if(distances[current_pos]<=row_mins1[i])
            {
               row_mins2[i]=row_mins1[i];
               //row_mins2_index[i]=row_mins1_index[i];
               row_mins1[i]=distances[current_pos];
               row_mins1_index[i]=j;
            }
            else
            {
               //if not smaller than minimum, it might be smaller than the minimum of order 2
               //check and replace
               if(distances[current_pos]<=row_mins2[i])
               {
                  row_mins2[i]=distances[current_pos];
                  //row_mins2_index[i]=current_pos;
               }
            }
            
            //now evaluate col mins (j)
            //if we find a new minimum, then the set the vectors for minimum positions and values
            //to this new minimum, the old minimum and minimum positions are then the minimum of order 2
            if(distances[current_pos]<=col_mins1[j])
            {
               col_mins2[j]=col_mins1[j];
               //col_mins2_index[j]=col_mins1_index[j];
               col_mins1[j]=distances[current_pos];
               col_mins1_index[j]=i;
            }
            else
            {
               //if not smaller than minimum, it might be smaller than the minimum of order 2
               //check and replace
               if(distances[current_pos]<=col_mins2[j])
               {
                  col_mins2[j]=distances[current_pos];
                  //col_mins2_index[j]=current_pos;
               }
            }
         }
         
      }
      
      //col and row minimums of 1st and 2nd order are found
      //now determine which offers the largest maximum delta
      for(j=0; j < ncol_d[0]; j++)
      {
         //if demand==0 skip this column
         if(demand[j]==0)
         {continue;}
         
         current_delta = col_mins2[j] -col_mins1[j];
         if(current_delta >= max_delta)
         {
            min_position_max_delta_row= col_mins1_index[j];
            min_position_max_delta_col= j;
            max_delta = current_delta;
         }
      }
      for(i=0;i<nrow_s[0]; i++)
      {
         //if supply ==0, skip this row
         if(supply[i]==0)
         {continue;}
         
         current_delta = row_mins2[i] - row_mins1[i];
         if(current_delta >= max_delta)
         {
            min_position_max_delta_row= i;
            min_position_max_delta_col= row_mins1_index[i];
            max_delta = current_delta;
         }
      }
      
      //we have now determined the largest delta, and where the minimum is in that row or column
      // perform the allocation as constrained by either supply or demand
      
      //calculate which is the current offset for the column#
         col_offset= min_position_max_delta_col * nrow_s[0];
      //current position in the "matrix"
         current_pos = min_position_max_delta_row + col_offset;
         
      //if supply == demand, easiest case, allocate all
      if(supply[min_position_max_delta_row] == demand[min_position_max_delta_col])
      {
         distribution[current_pos]=demand[min_position_max_delta_col];
         
         sum_demand= sum_demand-demand[min_position_max_delta_col];
         zerocounter_demand++;
         
         sum_supply= sum_supply- demand[min_position_max_delta_col];
         zerocounter_supply++;
         
         supply[min_position_max_delta_row]=0;
         demand[min_position_max_delta_col]=0;
      }
      else
      {
         //obviously supply =/= demand
         //if supply > demand, all demand may be satisfied, but supply will remain
         if(supply[min_position_max_delta_row] > demand[min_position_max_delta_col])
         {
            distribution[current_pos] = demand[min_position_max_delta_col];
            
            sum_demand= sum_demand-demand[min_position_max_delta_col];
            zerocounter_demand++;
            
            sum_supply= sum_supply- demand[min_position_max_delta_col];
            
            supply[min_position_max_delta_row]=supply[min_position_max_delta_row]-demand[min_position_max_delta_col];
            demand[min_position_max_delta_col]=0;
         }
         else
         {
            //if supply < demand, all supply may be allocated, but demand will remain
            distribution[current_pos] = supply[min_position_max_delta_row];
            
            sum_demand= sum_demand-supply[min_position_max_delta_row];
                        
            sum_supply= sum_supply- supply[min_position_max_delta_row];
            zerocounter_supply++;
            
            demand[min_position_max_delta_col]=demand[min_position_max_delta_col]-supply[min_position_max_delta_row];
            supply[min_position_max_delta_row]=0;
            
         }
      }
      
   }
   
   //if we only have one row or one column left, only the minimum needs to be found
   //so which is it rows or columns
   if((zerocounter_supply==(nrow_s[0]-1))&(zerocounter_demand== (ncol_d[0]-1)))
   {
      //if only one row and one column remains, then we need to only allocate all supply and demand,
      //transportation problem is defined as balanced by definition
      //find the nonzero row/col
      for(j=0; j< ncol_d[0];j++)
      {
         if(demand[j]==0)
         {continue;}
         min_position_max_delta_col=j;
         break;
      }
      for(i=0;i<nrow_s[0]; i++)
      {
         if(supply[i]==0)
         {continue;}
         min_position_max_delta_row=i;
         break;
      }
      
      col_offset= min_position_max_delta_col * nrow_s[0];
      current_pos = min_position_max_delta_row + col_offset;
      
      //allocate everything
      distribution[current_pos] = demand[min_position_max_delta_col];
      sum_demand= sum_demand-demand[min_position_max_delta_col];
      zerocounter_demand++;
      
      sum_supply= sum_supply- demand[min_position_max_delta_col];
      zerocounter_supply++;
      
      supply[min_position_max_delta_row]=0;
      demand[min_position_max_delta_col]=0;
   }
   else
   {
      //if we only have one row left allocate the rest of the supply
      if(zerocounter_supply==(nrow_s[0]-1))
      {
         //find which row doesn't have zero
         for(i=0;i<nrow_s[0]; i++)
         {
            if(supply[i]==0)
            {continue;}
            min_position_max_delta_row=i;
            break;
         }
         
         while(zerocounter_demand!=(ncol_d[0]-1))
         {
            
            row_mins1[min_position_max_delta_row]= DBL_MAX;
            min_position_max_delta_col=-1;
            
            //find minimum in the remaining row
            for(j=0; j< ncol_d[0];j++)
            {
               if(demand[j]==0)
               {continue;}
               
               //calculate which is the current offset for the column#
                  col_offset= j * nrow_s[0];
               //current position in the "matrix"
                  current_pos = min_position_max_delta_row + col_offset;
               
               if(distances[current_pos]<=row_mins1[min_position_max_delta_row])
               {
                  row_mins1[min_position_max_delta_row]=distances[current_pos];
                  min_position_max_delta_col=j;
               }
               
            }
            
            col_offset= min_position_max_delta_col * nrow_s[0];
            current_pos = min_position_max_delta_row + col_offset;
            //allocate to that minimum
            if(supply[min_position_max_delta_row] == demand[min_position_max_delta_col])
            {
               distribution[current_pos]=demand[min_position_max_delta_col];
               
               sum_demand= sum_demand-demand[min_position_max_delta_col];
               zerocounter_demand++;
               
               sum_supply= sum_supply- demand[min_position_max_delta_col];
               zerocounter_supply++;
               
               supply[min_position_max_delta_row]=0;
               demand[min_position_max_delta_col]=0;
            }
            else
            {
               //obviously supply =/= demand
               //if supply > demand, all demand may be satisfied, but supply will remain
               if(supply[min_position_max_delta_row] > demand[min_position_max_delta_col])
               {
                  distribution[current_pos] = demand[min_position_max_delta_col];
                  
                  sum_demand= sum_demand-demand[min_position_max_delta_col];
                  zerocounter_demand++;
                  
                  sum_supply= sum_supply- demand[min_position_max_delta_col];
                  
                  supply[min_position_max_delta_row]=supply[min_position_max_delta_row]-demand[min_position_max_delta_col];
                  demand[min_position_max_delta_col]=0;
               }
               else
               {
                  //if supply < demand, all supply may be allocated, but demand will remain
                  distribution[current_pos] = supply[min_position_max_delta_row];
                  
                  sum_demand= sum_demand-supply[min_position_max_delta_row];
                              
                  sum_supply= sum_supply- supply[min_position_max_delta_row];
                  zerocounter_supply++;
                  
                  demand[min_position_max_delta_col]=demand[min_position_max_delta_col]-supply[min_position_max_delta_row];
                  supply[min_position_max_delta_row]=0;
                  
               }
            }
         }
      }
      else
      {
         //otherwise only one column is left, allocate the remaining demand
         //find which col doesn't have zero
         for(j=0;i<ncol_d[0]; j++)
         {
            if(demand[j]==0)
            {continue;}
            min_position_max_delta_col=j;
            break;
         }
         
         while(zerocounter_supply!=(nrow_s[0]-1))
         {
            col_mins1[min_position_max_delta_col]= DBL_MAX;
            min_position_max_delta_row=-1;
            
            //calculate which is the current offset for the column#
                  col_offset= min_position_max_delta_col * nrow_s[0];
            //find the minimum in the remaining column
            for(i=0;i<nrow_s[0]; i++)
            {
               if(supply[i]==0)
               {continue;}
               
               //current position in the "matrix"
                  current_pos = i + col_offset;
               if(distances[current_pos]<=col_mins1[min_position_max_delta_col])
               {
                  col_mins1[min_position_max_delta_col]=distances[current_pos];
                  min_position_max_delta_row=i;
               }
               
            }
            col_offset= min_position_max_delta_col * nrow_s[0];
            current_pos = min_position_max_delta_row + col_offset;
            //allocate to that minimum
            if(supply[min_position_max_delta_row] == demand[min_position_max_delta_col])
            {
               distribution[current_pos]=demand[min_position_max_delta_col];
               
               sum_demand= sum_demand-demand[min_position_max_delta_col];
               zerocounter_demand++;
               
               sum_supply= sum_supply- demand[min_position_max_delta_col];
               zerocounter_supply++;
               
               supply[min_position_max_delta_row]=0;
               demand[min_position_max_delta_col]=0;
            }
            else
            {
               //obviously supply =/= demand
               //if supply > demand, all demand may be satisfied, but supply will remain
               if(supply[min_position_max_delta_row] > demand[min_position_max_delta_col])
               {
                  distribution[current_pos] = demand[min_position_max_delta_col];
                  
                  sum_demand= sum_demand-demand[min_position_max_delta_col];
                  zerocounter_demand++;
                  
                  sum_supply= sum_supply- demand[min_position_max_delta_col];
                  
                  supply[min_position_max_delta_row]=supply[min_position_max_delta_row]-demand[min_position_max_delta_col];
                  demand[min_position_max_delta_col]=0;
               }
               else
               {
                  //if supply < demand, all supply may be allocated, but demand will remain
                  distribution[current_pos] = supply[min_position_max_delta_row];
                  
                  sum_demand= sum_demand-supply[min_position_max_delta_row];
                              
                  sum_supply= sum_supply- supply[min_position_max_delta_row];
                  zerocounter_supply++;
                  
                  demand[min_position_max_delta_col]=demand[min_position_max_delta_col]-supply[min_position_max_delta_row];
                  supply[min_position_max_delta_row]=0;
                  
               }
            }
            
         }
      }
      
   }
   
   }
   else
   {
      //printf("\n sums of supply =/= demand, transport problem must be balanced\n");
      //printf("\n sum_supply = %d =/= sum_demand = %d\n",sum_supply, sum_demand);
   }
   
   //down the memory-hole
   free(col_mins1);
   free(col_mins1_index);
   free(col_mins2);

   free(row_mins1);
   free(row_mins1_index);
   free(row_mins2);
}


//End code
