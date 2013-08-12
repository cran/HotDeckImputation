/*
   Do not delete!
	File name		Distmatrix.c
	Part of:		   HotDeckImputation (GNU R contributed package)
	Author:			Dieter William Joenssen
	Copyright:		Dieter William Joenssen
	Email:			Dieter.Joenssen@TU-Ilmenau.de
	Created:		   14 May 2013
	Last Update: 	10 August 2013
	Description:	C functions to usefull for Hot deck Imputation
*/

//Begin preamble
#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void c_impute_mean( double *data, int *n, int *m);//

void c_test_NA( double *data, int *min, int *nin);//

void reweight_data( double *data, int *n, int *m, double *weights);//

void c_dist_hot_deck( double *distances, double *data, int *n, int *m, int *donors,
						int *recipients, int *n_donors, int *n_recipients,double *p);//

void c_impute_NN_HD( double *data, int *n, int *m, int *recipients, int *n_recipients,//
						int *donors);
void c_colMeans( double *data, int *n, int *m, double *means);//

//End preamble

//Begin code

//function that imputes the mean of the complete cases in any attribute for the missing values in that attribute
void c_impute_mean(double *data, int *n, int *m)
{
	int current_n, current_m,current_pos_m;
	double mean = 0.0;
	int available_case_counter = 0;
	
	for( current_m = 0; current_m < m[0]; current_m++)
	{
		//find mean
      current_pos_m =(n[0]*current_m);
		for(current_n = 0; current_n < n[0]; current_n++)
		{
			if( !ISNA( data[current_n+current_pos_m] ))
			{
				mean = mean + data[current_n+current_pos_m];
				available_case_counter++;
			}
		}
		mean = mean / (available_case_counter);
		available_case_counter = 0;
		//impute mean
		for(current_n = 0; current_n < n[0]; current_n++)
		{
			if(ISNA(data[current_n+current_pos_m]))
			{
				data[current_n+current_pos_m]= mean;
			}
		}
		mean=0.0;
	}
	
}

//function that returns the missing data indicator matrix
// if is.na(X) then 1 else 0
void c_test_NA(double *data, int *min, int *nin)
{
	int n=nin[0], m=min[0];
	int current_n,current_m,current_pos_m;

	//loop for 1. variable from which distance is calculated
   
	for(current_m=0;current_m < m;current_m++)
	{
      current_pos_m=(n*current_m);
      for(current_n=0; current_n<n;current_n++)
		{
			if(ISNA(data[current_n+current_pos_m]))
			{data[current_n+current_pos_m]=1;}
			else
			{data[current_n+current_pos_m]=0;}
		
		}

	}
}

//function that reweights the data before distance calculation
//this reduces the required computations for calculating a weighted distance
//at the expense of cloneing the data matrix
void reweight_data(double *data, int *n, int *m,double *weights )
{
	int current_n,current_m,current_pos_m;

   for(current_m=0;current_m < m[0];current_m++)	
	{
      
      current_pos_m=(n[0]*current_m);
      
		for(current_n=0; current_n<n[0];current_n++)
		{
			data[current_n+current_pos_m]=data[current_n+current_pos_m] * weights[current_m];
	
		}

	}
}

// calculates the distances required for hot deck imputation
// minkovski_factor needs to be implemented, currently only manhattan distance used
/*
void c_dist_hot_deck(
				double *distances,
				double *data,
				int *n,
				int *m,
				int *donors,
				int *recipients,
				int *n_donors,
				int *n_recipients)
{
	int current_donor,current_recipient,current_m,current_pos_m;
	int i = 0;
	int m_denom = m[0];
	for(current_recipient = 0;current_recipient<n_recipients[0];current_recipient++)
	{
		for(current_donor=0; current_donor<n_donors[0];current_donor++)
		{
			for(current_m=0;current_m < m[0];current_m++)
			{
            current_pos_m=(n[0] * current_m);
				if(ISNA(data[recipients[current_recipient] + current_pos_m]) |
				ISNA(data[donors[current_donor] + current_pos_m]))
				{
					m_denom = m_denom -1;
				}
				else
				{
					distances[i]= distances[i] + 
					fabs(data[recipients[current_recipient] + current_pos_m] - 
					data[donors[current_donor] + current_pos_m]);
				}
			}
			if(m_denom == 0)
			{distances[i] = DBL_MAX;}
			else
			{distances[i] = distances[i] * m[0] / m_denom;}
			i++;
			m_denom = m[0];
		}
	}
}
*/
//calculates required distance matrix
void c_dist_hot_deck(
      		double *distances,
				double *data,
				int *n,
				int *m,
				int *donors,
				int *recipients,
				int *n_donors,
				int *n_recipients,
            double *p)
{
	int current_donor,current_recipient,current_m,current_pos_m;
	int i = 0;
	int m_denom = m[0];
   
   
   
   if(p[0]==1)//manhattan distance
   {
      for(current_recipient = 0;current_recipient<n_recipients[0];current_recipient++)
      {
   		for(current_donor=0; current_donor<n_donors[0];current_donor++)
   		{
   			for(current_m=0;current_m < m[0];current_m++)
   			{
               current_pos_m=(n[0] * current_m);
   				if(ISNA(data[recipients[current_recipient] + current_pos_m]) |
   				ISNA(data[donors[current_donor] + current_pos_m]))
   				{
   					m_denom = m_denom -1;
   				}
   				else
   				{
   					distances[i]= distances[i] + 
   					fabs(data[recipients[current_recipient] + current_pos_m] - 
   					data[donors[current_donor] + current_pos_m]);
   				}
   			}
   			if(m_denom == 0)
   			{distances[i] = DBL_MAX;}
   			else
   			{distances[i] = distances[i] * m[0] / m_denom;}
   			i++;
   			m_denom = m[0];
   		}
   	}
   }
   if(p[0]==2)//euklidean
   {
      for(current_recipient = 0;current_recipient<n_recipients[0];current_recipient++)
      {
   		for(current_donor=0; current_donor<n_donors[0];current_donor++)
   		{
   			for(current_m=0;current_m < m[0];current_m++)
   			{
               current_pos_m=(n[0] * current_m);
   				if(ISNA(data[recipients[current_recipient] + current_pos_m]) |
   				ISNA(data[donors[current_donor] + current_pos_m]))
   				{
   					m_denom = m_denom -1;
   				}
   				else
   				{
   					distances[i]= distances[i] + 
   					(data[recipients[current_recipient] + current_pos_m] - 
   					data[donors[current_donor] + current_pos_m]) *
                  (data[recipients[current_recipient] + current_pos_m] - 
      				data[donors[current_donor] + current_pos_m]);
   				}
   			}
   			if(m_denom == 0)
   			{distances[i] = DBL_MAX;}
   			else
   			{distances[i] = distances[i] * m[0] / m_denom;}
   			i++;
   			m_denom = m[0];
   		}
   	}
   }
   if((p[0]!=2)&(p[0]!=1))
   {
      for(current_recipient = 0;current_recipient<n_recipients[0];current_recipient++)
      {
      	for(current_donor=0; current_donor<n_donors[0];current_donor++)
   		{
   			for(current_m=0;current_m < m[0];current_m++)
   			{
               current_pos_m=(n[0] * current_m);
   				if(ISNA(data[recipients[current_recipient] + current_pos_m]) |
   				ISNA(data[donors[current_donor] + current_pos_m]))
   				{
   					m_denom = m_denom -1;
   				}
   				else
   				{
   					distances[i]= distances[i] + 
   					pow(fabs((data[recipients[current_recipient] + current_pos_m] - 
   					data[donors[current_donor] + current_pos_m])),p[0]);
   				}
   			}
   			if(m_denom == 0)
   			{distances[i] = DBL_MAX;}
   			else
   			{distances[i] = distances[i] * m[0] / m_denom;}
   			i++;
   			m_denom = m[0];
   		}
   	}
   }
}

//performs actual imputation of values
void c_impute_NN_HD(double *data, int *n, int *m, int *recipients, int *n_recipients, int *donors)
{
   int counter, current_m,current_pos_m,recip_pos;
   for(current_m = 0; current_m < m[0];current_m++)
   {
      //precalc the offset due to variable position
      current_pos_m=(n[0] * current_m);
      for(counter = 0; counter<n_recipients[0];counter++)
      {
         //precalc total recipient position
         recip_pos=recipients[counter] + current_pos_m;
         if(ISNA(data[recip_pos]))
         {
            data[recip_pos] = data[donors[counter] + current_pos_m];
         }
      }
   }
   
/*   legacy code
   for(counter = 0; counter<n_recipients[0];counter++)
   {
      for(current_m = 0; current_m < m[0];current_m++)
      {
         if(ISNA(data[recipients[counter] + (n[0] * current_m)]))
         {
            data[recipients[counter] + (n[0] * current_m)] = data[donors[counter] + (n[0] * current_m)];
         }
      }
   }
   */
}

//calculates the columb means of complete cases
void c_colMeans(double *data, int *n, int *m, double *means)
{
   int current_n, current_m,current_pos_m;
   int available_case_counter = 0;
	
	for( current_m = 0; current_m < m[0]; current_m++)
	{
      current_pos_m=(n[0]*current_m);
		//find mean
		for(current_n = 0; current_n < n[0]; current_n++)
		{
			if( !ISNA( data[current_n+current_pos_m] ))
			{
				means[current_m] = means[current_m] + data[current_n+current_pos_m];
				available_case_counter++;
			}
		}
		means[current_m] = means[current_m] / (available_case_counter);
      available_case_counter = 0;
	}
}

//End code
