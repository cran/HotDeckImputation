\name{impute.NN_HD}
\alias{impute.NN_HD}

\title{
The Nearest Neighbor Hot Deck Algorithms
}
\description{
A comprehensive function that performs nearest neighbor hot deck imputation.  
Aspects such as variable weighting, distance types, and donor limiting are implemented.  
New concepts such as the optimal distribution of donors are also available.  
}

\usage{
impute.NN_HD(DATA = NULL, distance = "man", weights = "range", attributes = "sim",
 comp = "rw_dist", donor_limit = Inf, optimal_donor = "no",
 list_donors_recipients = NULL, diagnose = NULL)
}


\arguments{
  \item{DATA}{Data containing missing values. Must either be \code{data.frame}, then factors and strings will be recoded using \code{model.matrix} or Will be coerced by \code{data.matrix}.}
  \item{distance}{Distance type to use when searching for the nearest neighbor.  See the details section for options.}
  \item{weights}{Weights by which the variables should be scaled.  See the details section for options.}
  \item{attributes}{Determines how attributes should be handled.  Currently only "sim", meaning donor and recipient pools are disjoint, is implemented.}
  \item{comp}{Defines the compensation of missing values for distance calculation.  See the details section for options.}
  \item{donor_limit}{Limits how often a donor may function as such.  See the details section for options.}
  \item{optimal_donor}{Defines how the optimal donor is found when a donor limit is used.  See the details section for options.}
  \item{list_donors_recipients}{Option for manually specifying the donor and recipient pools via a list with components "donors" and "recipients".}
  \item{diagnose}{Option to recover the generated distances and donor-recipient matches. See details section for usage.}
}
\details{
\strong{argument:} \code{distance} can be defined as:
\itemize{
  \item \emph{numeric matrix}, donors x recipients distance matrix
  \item \emph{numeric length = 1}, Minkovski parameter
  \item \emph{string length = 1}, distance metric to be used:
  \itemize{
  \item \emph{"man"}, Manhattan distance
  \item \emph{"eukl"}, Euclidean distance
  \item \emph{"tscheb"}, Chebyshev distance
  \item \emph{"mahal"}, Mahalanobis distance (covariance matrix calculated after missing data compensation, incompatible with comp="rw_dist")
  }
}

\strong{argument:} \code{weights} can be defined as:
\itemize{
  \item \emph{string length = 1}, reweighting type "range", "sd", "var", or "none"
  \item \emph{numeric length = 1}, one numeric weight for all variables
  \item \emph{string vector}, like option 1, only different type for each variable (not implemented)
  \item \emph{numeric vector}, like option 2, only different numeric weight for each variable
  \item \emph{list}, mixture of options 3 and 4 (not implemented)
}

\strong{argument:} \code{comp} can be defined as:
\itemize{
  \item \emph{"rw_dist"}, total distance is reweighted by number of distances that may be computed
  \item \emph{"mean"}, mean imputation is performed before distance calculation
  \item \emph{"rseq"}, random hot deck imputation, each variable sequentially (uses impute.SEQ_HD)
  \item \emph{"rsim"}, random hot deck imputation, each variable simultaneously (not implemented)
}
\strong{argument:} \code{donor_limit} is a single number interpreted by its range:
\itemize{
  \item \emph{(0,1)}, dynamic donor limit, i.e., .5 means any 1 donor may serve up to 50\% of all recipients, rounded up if fractional
  \item \emph{[1,Inf)}, static donor limit, i.e., 2 means any 1 donor may serve up to 2 recipients, fractional parts discarded
  \item \emph{Inf}, no donor limit
}
\strong{argument:} \code{optimal_donor} is a single string interpreted by its value:
\itemize{
  \item \emph{"no"}, donor-recipient matching is performed in order by which the recipients appear in the data (fastest)
  \item \emph{"rand"}, donor-recipient matching is performed in a random recipient-order
  \item \emph{"mmin"}, donor-recipient matching is performed by the matrix minimum method (sequence independent)
  \item \emph{"modifvam"}, donor-recipient matching is performed by a modified (only columns considered) Vogel's approximation method (sequence independent)
  \item \emph{"vam"}, donor-recipient matching is performed by the Vogel's approximation method (sequence independent)
  \item \emph{"odd"}, donor-recipient matching is performed via the optimal donor distribution method (sequence independent, best results)
}
\strong{argument:} \code{diagnose} should be:
\itemize{
  \item \emph{\code{NULL}}, no diagnostics will be returned.
  \item \emph{character string}, desired variable name under which the diagnostics will be saved to \code{.GlobalEnv}. The following character strings will however default to \code{NULL} with a warning:
  "if", "else", "repeat", "while", "function", "for", "in", "next", "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_", "NA_complex_", "NA_character_", "c", "q", "s", "t", "C", "D", "F", "I", "T"
  \item \emph{anything else}, defaults to \code{NULL} with a warning.
  }
Should be a character string of the desired variable name which will be created in \code{.GlobalEnv}
}
\value{
An imputed data matrix the same size as the input \code{DATA}.
If the \code{diagnose} option is used correctly, a list containing the following components will be created in the workspace:
\item{distances }{the donor-recipient distance matrix used for matching}
\item{list_donors_recipients }{the resultant recipient-donor matches}
}

\references{
Andridge, R.R. and Little, R.J.A. (2010) A Review of Hot Deck Imputation for Survey Non-response. \emph{International Statistical Review}. \bold{78}, 40--64. 

Bankhofer, U. and Joenssen, D.W. (2014) On Limiting Donor Usage for Imputation of Missing Data via Hot Deck Methods. In: M. Spiliopoulou, L. Schmidt-Thieme, and R. Jannings (Eds.): \emph{Data Analysis, Machine Learning and Knowledge Discovery. Studies in Classification, Data Analysis and Knowledge Organization}, 3--11. Berlin/Heidelberg: Springer.

Domschke, W. (1995) \emph{Logistik: Transport}. Munich: Oldenbourg. [in German]

Ford, B. (1983) An Overview of Hot Deck Procedures. In: W. Madow, H. Nisselson and I. Olkin (Eds.): \emph{Incomplete Data in Sample Surveys}. New York: Academic Press, 185--207.

Joenssen, D.W. (2015) Donor Limited Hot Deck Imputation: A Constrained Optimization Problem. In: B. Lausen, S. Krolak-Schwerdt, and M. B\"ohmer (Eds.): \emph{Data Science, Learning by Latent Structures, and Knowledge Discovery. Studies in Classification, Data Analysis and Knowledge Organization}, pages 319--328. Berlin/Heidelberg: Springer.

Joenssen, D.W. (2015) Hot-Deck-Verfahren zur Imputation fehlender Daten -- Auswirkungen des Donor-Limits. Ilmenau: Ilmedia. [in German, Dissertation]

Joenssen, D.W. and Bankhofer, U. (2012) Donor Limited Hot Deck Imputation: Effects on Parameter Estimation. \emph{Journal of Theoretical and Applied Computer Science}. \bold{6}, 58--70. 

Kalton, G. and Kasprzyk, D. (1986) The Treatment of Missing Survey Data. \emph{Survey Methodology}. \bold{12}, 1--16.

Sande, I. (1983) Hot-Deck Imputation Procedures. In: W. Madow, H. Nisselson and I. Olkin (Eds.): \emph{Incomplete Data in Sample Surveys}. New York: Academic Press, 339--349.
}

\author{Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}}

%%\note{%%  ~~further notes~~}

\seealso{\code{\link{impute.mean}}, \code{\link{match.d_r_vam}}, \code{\link{reweight.data}}
}
\examples{
#Set the random seed to an arbitrary number
set.seed(421)

#Generate random integer matrix size 10x4
Y<-matrix(sample(x=1:100,size=10*4),nrow=10)

#remove 5 values, ensuring one complete covariate and 5 donors
Y[-c(1:5),-1][sample(1:15,size=5)]<-NA

#Impute using various different (arbitrarily chosen) settings
impute.NN_HD(DATA=Y,distance="man",weights="var")

impute.NN_HD(DATA=Y,distance=2,weights=rep(.5,4),donor_limit=2,optimal_donor="mmin")

impute.NN_HD(DATA=Y,distance="eukl",weights=.25,comp="mean",donor_limit=1,
 optimal_donor="odd")
 
#Recover some diagnostics
impute.NN_HD(DATA=Y,distance="eukl",weights=.25,comp="mean",donor_limit=1,
 optimal_donor="odd",diagnose = "diagnostics")
# look at the diagnostics
 diagnostics
}

\keyword{ NA }
\keyword{ manip }
\keyword{ optimize }
\keyword{ multivariate }

