#imports e1071

### creator
DSC_FCM <- function(centers,
                    weighted = TRUE,
                    iter.max = 100,
                    method = c("cmeans", "ufcl"),
                    dist = c("euclidean" , "manhattan"),
                    m = 2,
                    min_weight = NULL,
                    description=NULL)

{

  method <- match.arg(method)
  if(!is.null(description)) desc <- description
  else if(weighted) desc <- "Fuzzy C-Means (weighted)"
  else desc <-"Fuzzy C-Means"

  structure(list(description = desc,
                 RObj = fcm_refClass$new(
                   centers = centers, weighted = weighted, iter.max = iter.max,
                   method = method, dist = dist, m = m,
                   min_weight = min_weight)),
            class = c("DSC_FCM","DSC_Macro","DSC_R","DSC"))
}


fcm_refClass <- setRefClass("fcm",
                            fields = list(
                              centers	    = "numeric",
                              weighted = "logical",
                              iter.max    = "numeric",
                              method   = "character",
                              dist = "character",
                              m = "numeric",
                              assignment  = "numeric",
                              membership = "matrix",
                              data      = "data.frame",
                              weights	    = "numeric",
                              clusterCenters = "data.frame",
                              clusterWeights = "numeric",
                              details      = "ANY",
                              min_weight   = "numeric"
                            ),
                            methods = list(
                              initialize = function(
                                centers = 3,
                                weighted = TRUE,
                                iter.max = 100,
                                method = c("cmeans", "ufcl"),
                                dist = c("euclidean", "manhattan"),
                                m = 2,
                                min_weight = NULL
                              ) {
                                centers  	<<- centers
                                weighted <<- weighted
                                iter.max	<<- iter.max
                                method   <<- match.arg(method)
                                dist <<- match.arg(dist)
                                m <<- m
                                assignment	<<- numeric()
                                membership <<- matrix()
                                weights	<<- numeric()
                                clusterWeights <<- numeric()
                                clusterCenters <<- data.frame()
                                data	<<- data.frame()

                                if(is.null(min_weight)) min_weight <<- 0
                                else min_weight <<- min_weight

                                .self
                              }
                            ),
)

fcm_refClass$methods(
  cluster = function(x, weight = rep(1,nrow(x)), ...) {

    if(nrow(x)==1)
      warning("DSC_FCM does not support iterative updating! Old data is overwritten.")

    ### filter weak clusters
    if(min_weight>0) {
      x <- x[weight>min_weight,]
      weight <- weight[weight>min_weight]
    }


    weights <<- weight
    data <<- x

    if(nrow(data)>centers) {
      if(weighted) fcm <- e1071::cmeans(x=data,
                                        weight=weights,
                                        centers=centers,
                                        iter.max = iter.max,
                                        method = method,
                                        dist = dist,
                                        m = m)
      else fcm <- e1071::cmeans(x = data,
                                centers = centers,
                                iter.max = iter.max,
                                method = method,
                                dist = dist,
                                m = m)

      assignment <<- fcm$cluster
      membership <<- fcm$membership
      clusterCenters <<- data.frame(fcm$centers)
      details <<- fcm
    } else {
      assignment <<- 1:nrow(data)
      membership <<- diag(nrow(data))
      clusterCenters <<- x
      details <<- NULL
    }

    ## Cluster Weights = sum of (membership of example in cluster times weight of example)
    clusterWeights <<- sapply(1:nrow(clusterCenters), FUN =
                                function(i, j) sum(membership[j,i]*weights[j]), j = 1:nrow(data))

  },

  get_macroclusters = function(...) { clusterCenters },
  get_macroweights = function(...) { clusterWeights },

  get_microclusters = function(...) { data },
  get_microweights = function(x) { weights },

  get_membershipMatrix = function(...) { membership },

  microToMacro = function(micro=NULL, ...){
    if(is.null(micro)) micro <- 1:nrow(data)
    structure(assignment[micro], names=micro)
  }
)
