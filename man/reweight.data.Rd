\name{reweight.data}
\alias{reweight.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reweight Variables
}
\description{
Reweighting function to preprocess a data matrix prior to Minkovski distance calculation.
}
\usage{
reweight.data(DATA = NULL, weights = NULL, minkovski_factor = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
Data that should be reweighted.
}
  \item{weights}{
Numeric vector with length equal to the number of variables in \code{DATA}.
}
  \item{minkovski_factor}{
The desired Minkovski parameter that will be used for calculating the distances.
}
}
%%\details{%%  ~~ If necessary, more details than the description above ~~}
\value{
Returns a data matrix with the same dimensions as \code{DATA}.
}
%%\references{%% ~put references to the literature/web site here ~}
\author{Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}}

%%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{impute.NN_HD}}}

\examples{
#Set the random seed to an arbitrary number
set.seed(421)

#Generate matrix of random integers
Y<-matrix(sample(0:9,replace=TRUE,size=6*3),nrow=6)

#choose variable variances
Weights<-1/apply(X=Y,MARGIN=2,FUN=var)

#reweight data for faster Euclidean distance calculation
reweight.data(DATA = Y, weights = Weights, minkovski_factor = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

