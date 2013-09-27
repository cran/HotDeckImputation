\name{impute.CPS_SEQ_HD}
\alias{impute.CPS_SEQ_HD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CPS Sequential Hot-Deck Imputation
}
\description{
Resolves missing data by the CPS sequential Hot-Deck Imputation.
}
\usage{
impute.CPS_SEQ_HD(DATA = NULL, covariates = NULL, initialvalues = 0,
                  navalues = NA, modifyinplace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{Data containing missing values. Should be a matrix of integers.}
  \item{covariates}{Vector containing the covariates (columns that should be used to create the imputation classes). If \code{is.null(covariates) | length(covariates)==0} this function defaults to \code{impute.SEQ_HD}. See Section: Note for further Details.}
  \item{initialvalues}{The initial values for the start-up process of the imputation. Should be \code{"integer"} and \code{length(initialvalues)==1 | length(initialvalues)==dim(DATA)[2]}. The default of \code{0} is not normally a good value.}
  \item{navalues}{NA code for each variable that should be imputed. Should be \code{"integer"} and \code{length(initialvalues)==1 | length(initialvalues)==dim(DATA)[2]}. Default is R's NA value.}
  \item{modifyinplace}{Should \code{DATA} be modified in place? (See the Section: Warning.) If not, a copy is made.}
 }
 
 \section{Warning }{If \code{modifyinplace == FALSE} \code{DATA} or rather the variable supplied is edited directly! This is significantly faster if the data set is large.} 
\details{
This function imputes the missing values in any variable by creating imputation classes and then replicating the most recently observed value in the class and variable. Imputation classes are created by the adjustment cell method.
}
\value{
An imputed data matrix the same size as the input \code{DATA}.
}
\references{
Hanson, R.H. (1978) The Current Population Survey: Design and Methodology. \emph{Technical Paper No. 40 }. U.S. Bureau of the Census.

Joenssen, D.W. and Bankhofer, U. (2012) Donor Limited Hot Deck Imputation: Effects on Parameter Estimation. \emph{Journal of Theoretical and Applied Computer Science}. \bold{6}, 58--70. 

Joenssen, D.W. and Muellerleile, T. (2014) Fehlende Daten bei Data-Mining. \emph{HMD Praxis der Wirtschaftsinformatik}. \bold{51}, 458--468, 2014. doi: 10.1365/s40702-014-0038-8 [in German]
}
\author{Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}}

\note{
This is a very fast imputation method. Only one pass of the data is needed. With the use of proper covariates, data may be missing MAR.
Covariates should be complete (not missing data). If not, \code{NA} will be used for building classes. This may or may not be appropriate for the data. The presence of missing values in the covariates in not checked.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{impute.SEQ_HD}}, \code{\link{impute.mean}}, \code{\link{impute.NN_HD}}}
\examples{
#Set the random seed to an arbitrary number
set.seed(421)

n<-1000
m<-3
pmiss<-.1

#Generate matrix of random integers and 2 binary covariates
Y<-cbind(matrix(sample(0:1,replace=TRUE,size=n*2),nrow=n),
		 matrix(sample(0:9,replace=TRUE,size=n*m),nrow=n))

#generate missing values, MCAR, in all but the first two columns
Y[,-c(1,2)][sample(1:length(Y[,-c(1,2)]),
				   size=floor(pmiss*length(Y[,-c(1,2)])))]<-NA

#perform the sequential imputation Y within the 
#classes created by cross-classifying variables 1 and 2
impute.CPS_SEQ_HD(DATA=Y,covariates=c(1,2),initialvalues=0, navalues=NA, modifyinplace = FALSE)


####an example highlighting the modifyinplace option
#using cbind to show the results of the function and the intial data next to another
cbind(impute.CPS_SEQ_HD(DATA=Y,covariates=c(1,2),initialvalues=0,
                        navalues=NA, modifyinplace = FALSE),Y)
#notice that columns 8-10 (representing Y) still have missing data

#same procedure, except modifyinplace is set to TRUE
cbind(impute.CPS_SEQ_HD(DATA=Y,covariates=c(1,2),initialvalues=0,
                        navalues=NA, modifyinplace = TRUE),Y)
#notice that columns 8-10 (representing Y) are identical to columns 3-5, 
#Y has (and any Variables pointing to the same object have) been directly modified.
}

\keyword{ NA }
\keyword{ manip }
\keyword{ optimize }
\keyword{ multivariate }
\keyword{ impute }
