\name{impute.SEQ_HD}
\alias{impute.SEQ_HD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequential Hot-Deck Imputation
}
\description{
Resolves missing data by sequential Hot-Deck Imputation.
}
\usage{
impute.SEQ_HD(DATA = NULL, initialvalues = 0, navalues = NA,
              modifyinplace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{Data containing missing values. Should be a matrix of integers.}
  \item{initialvalues}{The initial values for the start-up process of the imputation. Should be \code{"integer"} and \code{length(initialvalues)==1 | length(initialvalues)==dim(DATA)[2]}. The default of \code{0} is not normally a good value.}
  \item{navalues}{NA code for each variable that should be imputed. Should be \code{"integer"} and \code{length(initialvalues)==1 | length(initialvalues)==dim(DATA)[2]}. Default is R's NA value.}
  \item{modifyinplace}{Should \code{DATA} be modified in place? (See the Section: Warning.) If not, a copy is made.}
 }
 
 \section{Warning }{If \code{modifyinplace == FALSE} \code{DATA} or rather the variable supplied is edited directly! This is significantly faster if the data set is large.} 
\details{
This function imputes the missing values in any variable by replicating the most recently observed value in that variable. 
}
\value{
An imputed data matrix the same size as the input \code{DATA}.
}
\references{
Hanson, R.H. (1978) The Current Population Survey: Design and Methodology. \emph{Technical Paper No. 40 }. U.S. Bureau of the Census.

Joenssen, D.W. (2015) Hot-Deck-Verfahren zur Imputation fehlender Daten -- Auswirkungen des Donor-Limits. Ilmenau: Ilmedia. [in German, Dissertation]

Joenssen, D.W. and Bankhofer, U. (2012) Donor Limited Hot Deck Imputation: Effects on Parameter Estimation. \emph{Journal of Theoretical and Applied Computer Science}. \bold{6}, 58--70. 

Joenssen, D.W. and Muellerleile, T. (2014) Fehlende Daten bei Data-Mining. \emph{HMD Praxis der Wirtschaftsinformatik}. \bold{51}, 458--468, 2014. doi: 10.1365/s40702-014-0038-8 [in German]
}
\author{Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}}

\note{
This is by far the fastest imputation method. Only one pass of the data is needed. However, no covariate information is used, thus only leads to good results if the data are missing MCAR.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{impute.CPS_SEQ_HD}}, \code{\link{impute.mean}}, \code{\link{impute.NN_HD}}}
\examples{
#Set the random seed to an arbitrary number
set.seed(421)

n<-1000
m<-5
pmiss<-.1


#Generate matrix of random integers
Y<-matrix(sample(0:9,replace=TRUE,size=n*m),nrow=n)

#generate 6 missing values, MCAR, in all but the first row
Y[-1,][sample(1:length(Y[-1,]),size=floor(pmiss*length(Y[-1,])))]<-NA

#perform the sequential imputation Y
impute.SEQ_HD(DATA=Y,initialvalues=0, navalues=NA, modifyinplace = FALSE)

####an example highlighting the modifyinplace option
#using cbind to show the results of the function and the intial data next to another
cbind(impute.SEQ_HD(DATA=Y,initialvalues=0, navalues=NA, modifyinplace = FALSE),Y)
#notice that columns 6-10 (representing Y) still have missing data

#same procedure, except modifyinplace is set to TRUE
cbind(impute.SEQ_HD(DATA=Y,initialvalues=0, navalues=NA, modifyinplace = TRUE),Y)
#notice that columns 6-10 (representing Y) are identical to columns 1-5, 
#Y has (and any Variables pointing to the same object have) been directly modified.

}

\keyword{ NA }
\keyword{ manip }
\keyword{ optimize }
\keyword{ multivariate }
\keyword{ impute }
