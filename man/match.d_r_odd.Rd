\name{match.d_r_odd}
\alias{match.d_r_odd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Donor-Recipient Matching via Optimal Donor Distribution Method
}
\description{
A function that performs the optimal donor distribution method.  GLPK is used to find the optimal solution to the integer program.
}
\usage{
match.d_r_odd(distance = NULL, recipients=NULL, donors=NULL, donor_limit=NULL)
}
\arguments{
  \item{distance}{A distance matrix of dimensions length(donors) * length(recipients)}
  \item{recipients}{A vector of object (row) numbers from the original data matrix, indicating which objects require imputation.}
  \item{donors}{A vector of object (row) numbers from the original data matrix, indicating which objects may be used for imputation.}
  \item{donor_limit}{A vector of length(donors) detailing how often any one donor may be used.}
}

\value{A matrix of dimensions length(recipients) * 2.  The first column, named recipients, is equivalent to \code{recipients}.  The second column, named donors, contains the donor that is matched to any of the recipients.}
\references{
Domschke, W. (1995) \emph{Logistik: Transport}. Munich: Oldenbourg. [in German]

Joenssen, D.W. (2015) Donor Limited Hot Deck Imputation: A Constrained Optimization Problem. In: B. Lausen, S. Krolak-Schwerdt, and M. B\"ohmer (Eds.): \emph{Data Science, Learning by Latent Structures, and Knowledge Discovery. Studies in Classification, Data Analysis and Knowledge Organization}, pages 319--328. Berlin/Heidelberg: Springer.

Joenssen, D.W. (2015) Hot-Deck-Verfahren zur Imputation fehlender Daten -- Auswirkungen des Donor-Limits. Ilmenau: Ilmedia. [in German, Dissertation]

}
\author{Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}}

\seealso{\code{\link{impute.NN_HD}}, \code{\link{match.d_r_vam}}}


\examples{
#Set the random seed to an arbitrary number
set.seed(421)

#Set up a random distance matrix
   ndonor=20
   nrecip=20
   distance<-matrix(sample(1:100,replace=TRUE,size=ndonor*nrecip),nrow=ndonor,ncol=nrecip)
#Name donors and recipients
   donors<-1:nrow(distance)
   recipients<-(nrow(distance)+1):(nrow(distance)+ncol(distance))
   colnames(distance)<-recipients
   rownames(distance)<-donors
#Set up the donor limit variable
   donor_limit<-rep(1,length(donors))

#perform the matching
match.d_r_odd(distance=distance,recipients=recipients,donors=donors,
   donor_limit=donor_limit)
}
\keyword{ optimize }

