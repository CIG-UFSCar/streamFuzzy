\name{DSC_FCM}
\alias{DSC_FCM}
\title{
Fuzzy C-Means Macro-clusterer
}
\description{
Class implements the fuzzy c-means algorithm for reclustering a set of micro-clusters.
}
\usage{
DSC_FCM(centers, weighted = TRUE, iter.max = 100,
        method = c("cmeans", "ufcl"),
        dist = c("euclidean", "manhattan"),
        m = 2, min_weight = NULL, description = NULL)
}
\arguments{
  \item{centers}{
Number of clusters or initial values for cluster centres.
}
  \item{weighted}{
Use a weighted fuzzy c-means.
}
  \item{iter.max}{
Maximum number of iterations.
}
  \item{method}{
If \dQuote{cmeans}, then we have the c-means fuzzy clustering method, if \dQuote{ufcl} we have the on-line update. Abbreviations are also accepted.
}
  \item{dist}{
Must be one of the following: If \dQuote{euclidean}, the mean square error, if \dQuote{manhattan}, the mean absolute error is computed. Abbreviations are also accepted.
}
  \item{m}{
A number greater than 1 giving the degree of fuzzification.
}
  \item{min_weight}{
Micro-clusters with a weight less than this will be ignored for reclustering.
}
  \item{description}{
Optional character string to describe the clustering method.
}
}
\details{
Please refer to function cmeans in \strong{e1071} for more details on the algorithm.

Note that this clustering cannot be updated iteratively and every time it is used for (re)clustering, the old clustering is deleted.
}
\value{
An object of class DSC_FCM (subclass of DSC, DSC_R, DSC_Macro)
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Priscilla de Abreu Lopes
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[stream]{DSC}}, \code{\link[stream]{DSC_Macro}}, \code{\link[e1071]{cmeans}}
}
\examples{
stream <- DSD_Gaussians(k=3, noise=0)

# create micro-clusters via sampling
sample <- DSC_Sample(k=20)
update(sample, stream, 500)
sample

# recluster micro-clusters
fcm <- DSC_FCM(centers=3)
recluster(fcm, sample)
plot(fcm, stream, type="both")

# For comparison we use fuzzy c-means directly to cluster data
# Note: fuzzy c-means is not a data stream clustering algorithm
fcm <- DSC_FCM(centers=3)
update(fcm, stream, 500)
plot(fcm, stream)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
