\name{HotDeckImputation-package}
\alias{HotDeckImputation-package}
\alias{HotDeckImputation}
\docType{package}
\title{Hot Deck Imputation Methods for Missing Data}
\description{
This package provides hot deck imputation methods to resolve missing data.  
Methods provided are popular in survey methodology, mostly used in the context of large national statistics, but are also finding their way to data mining due to their computational simplicity.  
A key aspect of this package is the implementation of the commonly advocated donor-limit.
}
\details{
\tabular{ll}{
Package: \tab HotDeckImputation\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2014-09-27\cr
License: \tab GPL-3\cr
}
\code{HotDeckImputation} is the ever expanding implementation of hot deck imputation methods, such as the nearest neighbor, the CPS-sequential and random hot deck.  
The package aims to be comprehensive in the functionality provided, covering key aspects of hot deck imputation not found elsewhere.  
\cr
Currently implemented functions include:\cr
Nearest neighbor hot deck imputation.\cr
Sequential hot deck imputation.\cr
CPS sequential hot deck imputation.\cr
Development requests are welcome.
}
\author{
Dieter William Joenssen

Maintainer: Dieter William Joenssen <Dieter.Joenssen@googlemail.com>
}

\references{
Andridge, R.R. and Little, R.J.A. (2010) A Review of Hot Deck Imputation for Survey Non-response. \emph{International Statistical Review}. \bold{78}, 40--64. 

Bankhofer, U. and Joenssen, D.W. (2014) On Limiting Donor Usage for Imputation of Missing Data via Hot Deck Methods. In: M. Spiliopoulou, L. Schmidt-Thieme, and R. Jannings (Eds.): \emph{Data Analysis, Machine Learning and Knowledge Discovery. Studies in Classification, Data Analysis and Knowledge Organization}, 3--11. Berlin/Heidelberg: Springer.

Domschke, W. (1995) \emph{Logistik: Transport}. Munich: Oldenbourg. [in German]

Ford, B. (1983) An Overview of Hot Deck Procedures. In: W. Madow, H. Nisselson and I. Olkin (Eds.): \emph{Incomplete Data in Sample Surveys}. New York: Academic Press, 185--207.

Joenssen, D.W. (2015) Donor Limited Hot Deck Imputation: A Constrained Optimization Problem. In: B. Lausen, S. Krolak-Schwerdt, and M. B\"ohmer (Eds.): \emph{Data Science, Learning by Latent Structures, and Knowledge Discovery. Studies in Classification, Data Analysis and Knowledge Organization}, pages tba. Berlin/Heidelberg: Springer. (in press)

Joenssen, D.W. and Bankhofer, U. (2012) Donor Limited Hot Deck Imputation: Effects on Parameter Estimation. \emph{Journal of Theoretical and Applied Computer Science}. \bold{6}, 58--70. 

Joenssen, D.W. and Muellerleile, T. (2014) Fehlende Daten bei Data-Mining. \emph{HMD Praxis der Wirtschaftsinformatik}. \bold{51}, 458--468, 2014. doi: 10.1365/s40702-014-0038-8 [in German]

Kalton, G. and Kasprzyk, D. (1986) The Treatment of Missing Survey Data. \emph{Survey Methodology}. \bold{12}, 1--16.

Sande, I. (1983) Hot-Deck Imputation Procedures. In: W. Madow, H. Nisselson and I. Olkin (Eds.): \emph{Incomplete Data in Sample Surveys}. New York: Academic Press, 339--349.
}

\keyword{package}
\keyword{manip}
\keyword{NA}
\keyword{multivariate}

%%\examples{}
