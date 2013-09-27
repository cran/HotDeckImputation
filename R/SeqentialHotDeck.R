impute.SEQ_HD<-function(DATA=NULL,initialvalues=0, navalues=NA, modifyinplace = TRUE)
{
   imputeables<-seq_len(dim(DATA)[2])
   imputeables<-imputeables-1
   
   if(length(initialvalues)!=length(imputeables))
   {
      if(length(initialvalues)==1)
      {
         initialvalues <- rep(initialvalues,length(imputeables))
      }
      else
      {
         stop("length(initialvalues) must equal 1 or dim(DATA)[2]\n")
      }
   }
   if(length(navalues)!=length(imputeables))
   {
      if(length(navalues)==1)
      {
         navalues <- rep(navalues,length(imputeables))
      }
      else
      {
         stop("length(navalues) must equal 1 or dim(DATA)[2]\n")
      }
   }
   
   mode(initialvalues)<-"integer"
   mode(navalues)<-"integer"
   
   if(storage.mode(DATA)=="integer")
   {
      if(!modifyinplace)
      {
        DATA<- .deepcopy(DATA)
      }
      .Call("SeqHD_INTEGER", DATA,  initialvalues, navalues)
   }
   else
   {
      stop("currently only matrices are acceptable\n")
   }
   return(DATA)
}


impute.CPS_SEQ_HD<-function(DATA=NULL,covariates=NULL,initialvalues=0, navalues=NA, modifyinplace = TRUE)
{
   if(is.null(covariates)|(length(covariates)==0))
   {
      #normal SEQ_HD
      impute.SEQ_HD(DATA,initialvalues, navalues, modifyinplace = TRUE)
   }
   
   imputeables<-seq_len(dim(DATA)[2])[-covariates]
   covariates<-covariates-1
   imputeables<-imputeables-1
   
   if(length(initialvalues)!=length(imputeables))
   {
      if(length(initialvalues)==1)
      {
         initialvalues <- rep(initialvalues,length(imputeables))
      }
      else
      {
         stop("length(initialvalues) must equal 1 or length(imputeables)\n")
      }
   }
   if(length(navalues)!=length(imputeables))
   {
      if(length(navalues)==1)
      {
         navalues <- rep(navalues,length(imputeables))
      }
      else
      {
         stop("length(navalues) must equal 1 or length(imputeables)\n")
      }
   }
   
   mode(covariates)<-"integer"
   mode(imputeables)<-"integer"
   mode(initialvalues)<-"integer"
   mode(navalues)<-"integer"
   
   if(storage.mode(DATA)=="integer")
   {
     if(!modifyinplace)
     {
       DATA<- .deepcopy(DATA)
     }
      .Call("CPS_SeqHD_INTEGER", DATA, covariates, imputeables, initialvalues, navalues)
   }
   else
   {
      stop("currently only matricies are acceptable\n")
   }
   return(DATA)
}


.deepcopy<-function(DATA=NULL)
{
  DATA<-.Call("deepcopy_c", DATA)
  return(DATA)
}