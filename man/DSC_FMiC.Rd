\name{DSC_FMiC}
\alias{DSC_FMiC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fuzzy Micro-Cluster - Online Clustering
}
\description{
Class implements the online management of a set of Fuzzy Micro-Clusters (FMiC).
}
\usage{
DSC_FMiC(maxMiC = 100, m = 2, theta = 0.8, Theta = 0.8, description = NULL)
}
\arguments{
  \item{maxMiC}{
Maximum number of Fuzzy Micro-Clusters in the summary structure.
}
  \item{m}{
A number greater than 1 giving the degree of fuzzification.
}
  \item{theta}{
Membership threshold for a point to be added to the existing summary structure.
}
  \item{Theta}{
Memembership threshold for two FMiC to be merged.
}
  \item{description}{
Optional character string to describe the clustering method.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class DSC_FMiC (subclass of DSC, DSC_R, DSC_Micro)
}
\references{
Lopes, P. A. and Camargo, H. A. 2017. FuzzStream: Fuzzy Data Stream Clustering Based on the Online-Offline Framework. In Proceedings of the 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE). Naples, Italy, 1--6.
}
\author{
Priscilla de Abreu Lopes
}
\note{
\bold{Known Issues:} the clustering function has a time optimization problem. :-/
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{FMiC}}, \code{\link[stream]{DSC}}, \code{\link[stream]{DSC_Micro}}
}
\examples{
stream <- DSD_BarsAndGaussians(noise=.05)
plot(stream)

fuzzyMicro <- DSC_FMiC()
update(fuzzyMicro, stream, 1000)
fuzzyMicro

# micro-clusters
nclusters(fuzzyMicro)
head(get_centers(fuzzyMicro))

# plot
plot(fuzzyMicro, stream)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
