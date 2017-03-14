\name{FMiC}
\alias{FMiC}
\title{
Fuzzy Micro-Cluster - Structure Definition
}
\description{
Defines object FMiC as the structure for a single Fuzzy Micro-Cluster (FMiC), creation, updating and auxiliary maintenance functions.
}
%\usage{
%}
%\arguments{
%}
\details{
A Fuzzy Micro-Cluster (FMiC) is defined as vector \eqn{(\overline{CF2^{x}}, \overline{CF1^{x}}, t, M)}{(CF2x, CF1x, t, M)}, where:

- \eqn{\overline{CF2^{x}}}{CF2x} is the quadratic sum of examples weighted by the membership values;

- \eqn{\overline{CF1^{x}}}{CF1x} is the linear sum of examples weighted by the membership values;

- \eqn{t}{t} is the timestamp for the newest example added;

- \eqn{M}{M} is the scalar cardinality of the micro-cluster (sum of membership values).

}
\value{
An object of class FMiC. Aside from general components described in \code{Details}, the centroid (c) for the FMiC is also available.
}
\references{
Lopes, P. A. and Camargo, H. A. 2017. FuzzStream: Fuzzy Data Stream Clustering Based on the Online-Offline Framework. In Proceedings of the 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE). Naples, Italy, 1--6.
}
\author{
Priscilla de Abreu Lopes
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{DSC_FMiC}}
}
\examples{
# set of points p
p <- data.frame(X = c(1,2,3), Y = c(4,5,6))

# creates a new fmic with first point in p and membership value of 1
fmic <- FMiC(p[1,], membership = 1)

# add the other points in p with 0.8 and 0.3 membership values respectively
fmic$addPoints(p[2:3,], membership = c(0.8,0.3))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
