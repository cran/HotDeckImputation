/*
   Do not delete!
   File name		minheuristics.c
	Part of:		   HotDeckImputation (GNU R contributed package)
	Author:			Dieter William Joenssen
	Copyright:		Dieter William Joenssen
	Email:			Dieter.Joenssen@TU-Ilmenau.de
	Created:		   10 August 2013
	Last Update: 	10 August 2013
	Description:	C functions to perform various minimum selection heuristics
*/

//Begin preamble
#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


// simple recipient donor matching based on min distance in each recipient instance
// is dependent on order of imputation
void c_match_d_r( int *recipient_donor_list, double *distance, int *recipient_order,
					double *donor_limit, int * n_donors, int *n_recipients);//

void c_colmin(int *distribution, double *distance, int *supply, int *demand, int *nrow_s, int *ncol_d);//

void c_matmin(int *distribution, double *distance, int *supply, int *demand, int *nrow_s, int *ncol_d);//
//End preamble

//Begin code

// simple recipient donor matching based on min distance in each recipient instance
// is dependent on order of imputation
void c_match_d_r(int *recipient_donor_list, double *distance, int *recipient_order, double *donor_limit, int * n_donors, int *n_recipients)
{

   // initialize variables
	int counter_donor, counter_recipient, current_recipient,current_pos_m;
	int current_min_dist_donor;
	double current_min_dist, current_dist;

	// for each recipient (i.e columb of the distance matrix)
	for(counter_recipient = 0; counter_recipient <n_recipients[0]; counter_recipient++)
	{
	
		// this is relevant in case of random recipient order
		// counter_recipient != actual recipient considered
		current_recipient = recipient_order[counter_recipient];

		// reset min dist for each recipient
		current_min_dist = DBL_MAX;
		// reset which donor has min dist
		current_min_dist_donor = -1;
		
      //precalc multiplication
      current_pos_m=(n_donors[0] * current_recipient);
      
		// repeat for each donor (i.e. row of distance matrix)
		for(counter_donor = 0; counter_donor< n_donors[0]; counter_donor++)
		{
			// if allocation may be made (i.e. donor limit not exhausted)
			if(donor_limit[counter_donor] > 0)
			{
				// the distance of donor recipient pair currently considered
				current_dist = distance[counter_donor + current_pos_m];

				// if the distance between the recipient and the current donor canidate is
				// smaller than previous smallest distance
				if(current_dist < current_min_dist)
				{
					// change min dist to current dist for future comparisons
					current_min_dist = current_dist;
					// change best donor index to this donor
					current_min_dist_donor = counter_donor;
				}
			}
		}

		// if a best donor is found
		if(current_min_dist_donor != -1)
		{			
			//reduce this donors limit by one
			if(donor_limit[current_min_dist_donor]!=DBL_MAX)
			{donor_limit[current_min_dist_donor] = donor_limit[current_min_dist_donor] - 1;}
			
			//save the pairing
			recipient_donor_list[current_recipient] = current_min_dist_donor+1;
		}
	}
	//still missing: resolve if no best donor is found, i.e. only inf available
   //wehn would this even happen?
}

void c_colmin(int *distribution, double *distance, int *supply, int *demand, int *nrow_s, int *ncol_d)
{
   int i=-1,j=-1,current_pos=0,zerocounter_supply=0,zerocounter_demand=0;
   int minval_i=-1;
   double minval = DBL_MAX;
   
   for(i=0;i< nrow_s[0];i++)
   {if(supply[i]==0){zerocounter_supply++;}}
   for(j=0;j<ncol_d[0];j++)
   {if(demand[j]==0){zerocounter_demand++;}}
   
   j=0;
   while((j<ncol_d[0])&(zerocounter_demand<ncol_d[0])&(zerocounter_supply<nrow_s[0]))
   {
      if(demand[j]==0)
      {
         j++;
         continue;
      }
      //find col minimum
      for(i=0;i< nrow_s[0];i++)
      {
         
         if(supply[i]==0)
         {continue;}
         current_pos=i+(nrow_s[0]*j);
         if(distance[current_pos]<=minval)
         {
            minval=distance[current_pos];
            minval_i=i;
         }
      }
      
      //minimum found, allocate all possible demand
      current_pos=minval_i+(nrow_s[0]*j);
      if(supply[minval_i]==demand[j])
      {
         distribution[current_pos]= supply[minval_i];
         supply[minval_i]=0;
         demand[j]=0;
         zerocounter_supply++;
         zerocounter_demand++;
         j++;
      }
      else
      {
         if(supply[minval_i]>demand[j])
         {
            distribution[current_pos]= demand[j];
            supply[minval_i]=supply[minval_i]-demand[j];
            demand[j]=0;
            zerocounter_demand++;
            j++;
         }
         else
         {
            distribution[current_pos]= supply[minval_i];
            demand[j]=demand[j]-supply[minval_i];
            zerocounter_supply++;
            supply[minval_i]=0;
         }
      }
      minval = DBL_MAX;
      minval_i=-1; 
   }

}

void c_matmin(int *distribution, double *distance, int *supply, int *demand, int *nrow_s, int *ncol_d)
{
   int i=-1,j=-1,zerocounter_supply=0,zerocounter_demand=0,sum_supply = 0,sum_demand = 0,current_pos=0;
   int minval_i=-1,minval_j=-1,pos_optimum=-1;
   double minval = DBL_MAX;
   
   for(i=0;i< nrow_s[0];i++)
   {
      if(supply[i]==0)
      {zerocounter_supply++;continue;}
      sum_supply = sum_supply + supply[i];
   }
   for(j=0;j<ncol_d[0];j++)
   {
      if(demand[j]==0)
      {zerocounter_demand++;continue;}
      sum_demand = sum_demand + demand[j];
   }
   
   while((sum_demand!=0)&(sum_supply!=0)&(zerocounter_demand<ncol_d[0])&(zerocounter_supply<nrow_s[0]))
   {
      minval = DBL_MAX;
      minval_i=minval_j=-1;
      
      //find minimum of remaining values
      for(i=0;i< nrow_s[0];i++)
      {
         if(supply[i]==0)
         {continue;}
         for(j=0;j<ncol_d[0];j++)
         {
            if(demand[j]==0)
            {continue;}
            current_pos=i+(nrow_s[0]*j);
            if(distance[current_pos]<=minval)
            {
               minval=distance[current_pos];
               minval_i=i;
               minval_j=j;
            }
         }
      }
      
      //continue here with supply "==" ">" "<" demand logic
      pos_optimum = minval_i+(nrow_s[0]*minval_j);
      //if supply == demand
      if(supply[minval_i]==demand[minval_j])
      {
         distribution[pos_optimum] = demand[minval_j];
         //update sums and zerocounters
         zerocounter_demand++;
         zerocounter_supply++;
         sum_demand = sum_demand - demand[minval_j];
         sum_supply = sum_supply - demand[minval_j];
         //update supply and demand values
         supply[minval_i]=0;
         demand[minval_j] = 0;
      }
      else
      {
         //if supply is > demand
         if(supply[minval_i]>demand[minval_j])
         {
            distribution[pos_optimum] = demand[minval_j];
            //update sums and zerocounters
            zerocounter_demand++;
            sum_demand = sum_demand - demand[minval_j];
            sum_supply = sum_supply - demand[minval_j];
            //update supply and demand values
            supply[minval_i]=supply[minval_i]-demand[minval_j];
            demand[minval_j] = 0;
         }
         else
         {
            distribution[pos_optimum] = supply[minval_i];
            //update sums and zerocounters
            zerocounter_supply++;
            sum_demand = sum_demand - supply[minval_i];
            sum_supply = sum_supply - supply[minval_i];
            //update supply and demand values
            demand[minval_j] = demand[minval_j] - supply[minval_i];
            supply[minval_i]=0;
         }
      }
   }

}

//End code
