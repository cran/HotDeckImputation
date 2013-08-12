#
#    Do not delete!
#  File name   	HotDeckImputation.R
#  Part of:		   HotDeckImputation (GNU R contributed package)
#  Author:			Dieter William Joenssen
#  Copyright:		Dieter William Joenssen
#  Email:			Dieter.Joenssen@TU-Ilmenau.de
#  Created:		   14 May 2013
#  Last Update: 	11 August 2013
#  Description:	R code for Package HotDeckImputation. Implemented functionions include following:
#                 Topics:
#                 -impute.NN_HD           ~ An extremely comprehensive implementation Nearest-Neighbor hot deck algorithms
#                 -impute.mean            ~ Imputes the mean value of the complete cases in a variable for the missing values
#                 -reweight.data          ~ Reweights data with given wieghts, preprocessing for distance calculation
#				  -match.d_r_vam          ~ performs the Vogel's approximation method matching for donors-recipients


##Description of impute.NN_HD
{
# data is data, requires numeric matrix
# distance is either:
#     numeric matrix       ---   predefined donors x recipients distance matrix
#                                (row/colnames must correspond to object number in data)
#     string length=1      ---   defining dinstance metric to be used ("man","eukl","tscheb","mahal")
#     numeric length=1       ---   defining minkovski parameter
# weights is either:
#     string length=1        ---   defining reweighting type ("range", "sd", "var", "none")
#     numeric length=1       ---   defining one numeric weight for all variables
#     string vector        ---   like option 1, only different type for each variable
#     numeric vector       ---   like option 2, only different numeric weight for each variable
#     list                 ---   mixture of options 3 and 4
# comp is:
#     string length=1     ---   defining compensation of missing values for distance calculation
#                                ("rw_dist","mean","rseq","rsim")
# donor_limit is:
#     numeric length=1      ---   Interpretaion of value is range dependent:
#                                - if value is in (0,1) then it is interpreted as a dynamic donor limit 
#                                  i.e. .5 means any 1 donor may serve 50% of all recipients
#                                  in case of fractional results, this is rounded up
#                                - if value is 1 or larger then it is interpreted as a static donor limit
#                                  i.e. 2 means any 1 donor may serve up to 2 recipients
#                                - Inf is unlimited donor usage
#                                - Fractional parts for numers larger than 1 are discarded
# optimal_donor:
#  string length=1        ---   defines how the optimal donor is found when a donor limit is used: (problem explained extensively elsewhere)
#                                - "no"      recipients are allocated their donor by order they apear in the data
#                                - "rand"    recipients are allocated their donor by random order
#                                - "mmin"    matrix minmum method is used to allocate donors to the recipients
#                                - "vam"     Vogel's approximation method is used to allocate donors to the recipients
#                                - "modi"    Modified Stepingstone algorithm is used
#                                ("no", "rand","mmin","vam","modi")
}

impute.NN_HD<-function(data=NULL,distance="man",weights="range",attributes="sim", comp="rw_dist",donor_limit=Inf,optimal_donor="no",
                list_donors_recipients = NULL)
{
   original_data<-data
   n<-dim(data)[1]
   m<-dim(data)[2]
   
   #test what type of weights and create appropriate vector
   if(length(weights)==1)
   {
      if(is.character(weights))
      {
         weights = match.arg(arg=weights,choices=c("none","range","sd","var"),several.ok=FALSE)
         if(weights=="none")
         {weights <-rep(1,m)}
         else
         {
            if(weights=="range")
            {weights <-1/apply(apply(data,MARGIN=2,range,finite=TRUE),MARGIN=2,diff)}
            else
            {
               if(weights =="sd")
               {weights <-1/apply(data,MARGIN=2,sd,na.rm=TRUE)}
               else
               {#must be == "var"
                  weights <-1/apply(data,MARGIN=2,var,na.rm=TRUE)
               }
            }
         }
      }
      else
      {
         ## weights must be 1 numeric
         if(is.numeric(weights))
         {weights<-rep(weights,m)}
         else
         {stop(paste("uniplemented weights used: ",weights))}
      }
   }
   else
   {
      #more than one value for weights is given
      if(!is.numeric(weights))
      {#if it is not numeric it must be either a character vector or list
         if(is.character(weights))
         {
            if(length(weights)==m)
            {stop("unimplemented feature, option 3 for weights\n")}
            else
            {
               stop("inproper usage of option 3 for weights, number of weights must equal number of attributes\n")
            }
            weights<-rep(1,m)
         }
         else
         {
            #must be list
            if(length(weights)==m)
            {stop("unimplemented feature, option 5 for weights\n")}
            else
            {
               stop("inproper usage of option 5 for weights, number of weights must equal number of attributes\n")
            }
            weights<-rep(1,m)
         }
      }
   }
   
   weights[weights==0]<-1
   
   if(attributes=="sim")
   {
      
      #if no distance matrix is specified, then donors and recipients need to be found
      if(is.null(list_donors_recipients))
      {
         MV_Indicator_M<-matrix(.C("c_test_NA",data=as.double(data),
                                   min=as.integer(ncol(data)),nin=as.integer(nrow(data)),
                                   NAOK=TRUE)$data,nrow=nrow(data))
         NA_counts <- rowSums(MV_Indicator_M)
         if(attributes == "sim")
         {
            donors <- which(NA_counts==0)
            list_donors_recipients<-list(donors=donors, recipients=(1:dim(MV_Indicator_M)[1])[-donors])
            rm(list=c("donors","MV_Indicator_M","NA_counts"))
         }
         else
         {
            stop("Unimplemented way of handling attributes used: \"seq \"")
            donors <- which(NA_counts<dim(MV_Indicator_M)[2])
            recipients <- which(rowSums(MV_Indicator_M)>0)
         }
      }
      
      #if not matrix, matrix must be computed
      if(!is.matrix(distance))
      {
         if(is.character(distance))
         {
            distance = match.arg(arg=distance,choices=c("man","eukl","tscheb","mahal"),several.ok=FALSE)
            if(distance == "man"){distance <- 1}
            if(distance == "eukl"){distance <- 2}
            if(distance == "tscheb"){stop("unimplemented feature, distance = \"tscheb\" \n")}
            if(distance == "mahal"){stop("unimplemented feature, distance = \"mahal\" \n")}

         }
         
         #reweight data
         if(!all(weights==1))
         {
            data<- reweight.data(data,weights,distance)
         }
         
         #print(list_donors_recipients)
         
         comp = match.arg(arg=comp,choices=c("rw_dist","mean","rseq","rsim"),several.ok=FALSE)
            #if(comp=="rw_dist")
            #{#no action needed}
            if(comp=="mean")
            {data<-impute.mean(data)}
            if(comp=="rseq")
            {stop("unimplemented feature, comp = \"rseq\" \n")}
            if(comp=="rsim")
            {stop("unimplemented feature, comp = \"rsim\" \n")}

         
         out <- .C("c_dist_hot_deck",
                   distances = as.double(rep(0,length(list_donors_recipients$donors)*length(list_donors_recipients$recipients))),
                   data = as.double(data),
                   n = as.integer(dim(data)[1]),
                   m = as.integer(dim(data)[2]),
                   donors = as.integer(list_donors_recipients$donors-1),
                   recipients = as.integer(list_donors_recipients$recipients-1),
                   n_donors = as.integer(length(list_donors_recipients$donors)),
                   n_recipients = as.integer(length(list_donors_recipients$recipients)),
                   p = as.double(distance),
                   NAOK=TRUE)
         
            distance<-matrix(out[[1]],nrow=length(list_donors_recipients$donors),byrow=FALSE)
            rownames(distance)<-list_donors_recipients$donors
            colnames(distance)<-list_donors_recipients$recipients
            
         }
      

      
      #make donor_limit static 
      if(donor_limit >0)
      {
         if(donor_limit< 1)
         {
            #dynamic donor limit used
            donor_limit <- ceiling(donor_limit * length(list_donors_recipients$recipients))
         }else{
            #static donor limit used
            donor_limit<-ceiling(donor_limit)
         }
      }else{
         stop(paste("inproper usage of donor_limit: ",donor_limit))}
      if(is.infinite(donor_limit))
      {donor_limit<-length(list_donors_recipients$recipients)}
      
      if((donor_limit*length(list_donors_recipients$donors))<length(list_donors_recipients$recipients))
      {stop(paste("Donor limit set too stringent. Number of donors * donor-limit must be >= Number of recipients.\n Compare:",
                  donor_limit*length(list_donors_recipients$donors),"vs.",length(list_donors_recipients$recipients)))
      }
      optimal_donor = match.arg(arg=optimal_donor,choices=c("no", "rand", "mmin","modifvam","vam", "modi"),several.ok=FALSE)
      ##match.donor_recipient
      if(optimal_donor == "no")
      {
         out <- .C("c_match_d_r",
                   recipient_donor_list = as.integer(rep(-1,length(list_donors_recipients$recipients))),
                   distance = as.double(distance),
                   recipient_order = as.integer( (1:dim(distance)[2])-1),
                   donor_limit = as.double(rep(donor_limit,dim(distance)[1])),
                   n_donors = as.integer(dim(distance)[1]),
                   n_recipients = as.integer(dim(distance)[2]),
                   NAOK=TRUE)
         list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[out[[1]]])
         
      }
      if(optimal_donor == "rand")
      {
         out <- .C("c_match_d_r",
                   recipient_donor_list = as.integer(rep(-1,length(list_donors_recipients$recipients))),
                   distance = as.double(distance),
                   recipient_order = as.integer(sample(1:dim(distance)[2])-1),
                   donor_limit = as.double(rep(donor_limit,dim(distance)[1])),
                   n_donors = as.integer(dim(distance)[1]),
                   n_recipients = as.integer(dim(distance)[2]),
                   NAOK=TRUE)
         list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[out[[1]]])
         
      }
      if(optimal_donor == "mmin")
      {
         
            distribution <- .C("c_matmin",
                               distribution = as.integer(rep(0,length(distance))),
                               distance = as.double(distance),
                               supply = as.integer(rep(donor_limit,dim(distance)[1])),
                               demand = as.integer(rep(1,length(list_donors_recipients$recipients))),
                               nrow_s =as.integer(length(list_donors_recipients$donors)),
                               ncol_d =as.integer(length(list_donors_recipients$recipients)),
                               NAOK=TRUE)$distribution
            
            distribution<-matrix(distribution,nrow=length(list_donors_recipients$donors),byrow=FALSE)
            donor_index<-apply(distribution,2,which.max)
            list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[donor_index])
         
      }
      if(optimal_donor == "modifvam")
      {
         distribution <- .C("c_modifVAM",
                            distribution = as.integer(rep(0,length(distance))),
                            distance = as.double(distance),
                            supply = as.integer(rep(donor_limit,dim(distance)[1])),
                            demand = as.integer(rep(1,length(list_donors_recipients$recipients))),
                            nrow_s =as.integer(length(list_donors_recipients$donors)),
                            ncol_d =as.integer(length(list_donors_recipients$recipients)),
                            NAOK=TRUE)$distribution
            
         distribution<-matrix(distribution,nrow=length(list_donors_recipients$donors),byrow=FALSE)
         donor_index<-apply(distribution,2,which.max)
         list_recip_donor<-cbind(recipient=list_donors_recipients$recipients,donor=list_donors_recipients$donors[donor_index])
         
      }
      if(optimal_donor == "vam")
      {
         list_recip_donor<-match.d_r_vam(distance = distance,recipients=list_donors_recipients$recipients,
                                         donors=list_donors_recipients$donors,donor_limit=rep(donor_limit,dim(distance)[1]))
         
      }
      if(optimal_donor == "modi")
      {stop("unimplemented feature, optimal_donor = \"modi\" \n")}
      
      #IMPUTE!
         data <- .C("c_impute_NN_HD",
                   data = as.double(original_data),
                   n = as.integer(dim(data)[1]),
                   m = as.integer(dim(data)[2]),
                   recipients = as.integer(list_recip_donor[,1]-1),
                   n_recipients = as.integer(dim(list_recip_donor)[1]),
                   donors = as.integer(list_recip_donor[,2]-1),
                   NAOK=TRUE)$data
         data<-matrix(data,nrow=n,ncol=m)
}
   else
   {stop("Unimplemented way of handling attributes used: \"seq \"")}
   
   return(data)
}

impute.mean <-function(data=NULL)
{
   out <- .C("c_impute_mean",
             data = as.double(data),
             n = as.integer(dim(data)[1]),
             m = as.integer(dim(data)[2]),
             NAOK=TRUE
             )
   return(matrix(out[[1]],ncol=dim(data)[2]))
}

reweight.data<-function(data=NULL,weights=NULL,minkovski_factor=1)
{
   weights <- (weights)^(1/minkovski_factor)
   out <- .C("reweight_data",
             data = as.double(data),
             n = as.integer(dim(data)[1]),
             m = as.integer(dim(data)[2]),
             weights = as.double(weights),
             NAOK=TRUE
             )
   return(matrix(out[[1]],ncol=dim(data)[2]))
}

# matrix, vector, vector, vector
match.d_r_vam<-function(distance = NULL,recipients=NULL,donors=NULL,donor_limit=NULL)
{
   
   demand_fictive_recipient <- sum(donor_limit) - length(recipients);
   if(demand_fictive_recipient == 0)
   {
      distribution <- .C("c_VAM",
                         distribution = as.integer(rep(0,length(distance))),
                         distance = as.double(distance),
                         supply = as.integer(donor_limit),
                         demand = as.integer(rep(1,length(recipients))),
                         nrow_s =as.integer(length(donors)),
                         ncol_d =as.integer(length(recipients)),
                         NAOK=TRUE)$distribution
      
      distribution<-matrix(distribution,nrow=length(donors),byrow=FALSE)
   }
   if(demand_fictive_recipient>0)
   {
      distance<-cbind(distance, rep(0,length(donors)))
      
      distribution <- .C("c_VAM",
                         distribution = as.integer(rep(0,length(distance))),
                         distance = as.double(distance),
                         supply = as.integer(donor_limit),
                         demand = as.integer(c(rep(1,length(recipients)),demand_fictive_recipient)),
                         nrow_s =as.integer(length(donors)),
                         ncol_d =as.integer(length(recipients)+1),
                         NAOK=TRUE)$distribution
      
      distribution<-matrix(distribution,nrow=length(donors),byrow=FALSE)
      distribution<-distribution[,-(ncol(distribution))]
      
   }
   if(demand_fictive_recipient<0)
   {
      distance<-rbind(distance, rep(0,length(recipients)))
      
      distribution <- .C("c_VAM",
                         distribution = as.integer(rep(0,length(distance))),
                         distance = as.double(distance),
                         supply = as.integer(c(donor_limit,demand_fictive_recipient*-1)),
                         demand = as.integer(rep(1,length(recipients))),
                         nrow_s =as.integer(length(donors)+1),
                         ncol_d =as.integer(length(recipients)),
                         NAOK=TRUE)$distribution
      
      distribution<-matrix(distribution,nrow=(length(donors)+1),byrow=FALSE)
      distribution<-distribution[-(nrow(distribution)),]
      
   }
   donor_index<-apply(distribution,2,which.max)
   return(cbind(recipient=recipients,donor=donors[donor_index]))
}