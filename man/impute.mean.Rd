\name{impute.mean}
\alias{impute.mean}

\title{Attribute Wise Mean Imputation}
\description{
This function imputes the column mean of the complete cases for the missing cases.  
Utilized by \code{impute.NN_HD} as a method for dealing with missing values in distance calculation.  
}
\usage{
impute.mean(data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data with missing values.}
}
%%\details{%%  ~~ If necessary, more details than the description above ~~}

\value{Returns an imputed data matrix with the same dimensions as \code{data}.}

\references{
Little, R.J.A and Rubin, D.B. (2002) \emph{Statistical Analysis with Missing Data}. Hoboken: Wiley.
}

\author{Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}}

%%\note{}

\seealso{\code{\link{impute.NN_HD}}}

\examples{
#Set the random seed to an arbitrary number
set.seed(421)

#Generate matrix of random integers
Y<-matrix(sample(0:9,replace=TRUE,size=6*3),nrow=6)

#generate 6 missing values, MCAR, in all but the first row
Y[-1,][sample(1:12,size=6)]<-NA

#Impute the colMeans of Y
impute.mean(data=Y)
}

\keyword{ NA }
\keyword{ manip }

