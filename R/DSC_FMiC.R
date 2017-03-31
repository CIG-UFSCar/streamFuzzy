## FMiC - Management
DSC_FMiC <- function(maxMiC = 100, m = 2, theta = 0.8, Theta = 0.8, description = NULL){

  if(!is.null(description)) desc <- description
  else desc <-"Fuzzy Micro-Cluster - Management"

  structure(list(description = desc,
                 RObj = fMicro_refClass$new(
                   maxMiC = maxMiC, m = m, theta = theta, Theta = Theta)),
            class = c("DSC_FMiC","DSC_Micro","DSC_R","DSC"))
}

fMicro_refClass <-
  setRefClass("fuzzyMicroCluster",
              fields = list(
                ## args
                maxMiC = "numeric",
                m = "numeric",
                theta = "numeric", # add to existing FMiC
                Theta = "numeric", # merge two FMiC
                ## store microcluster
                allFMiC = "list",
                microCenters = "data.frame", # needed for Macro-Clustering
                microWeights = "numeric" # needed for Macro-Clustering
              ),

              methods = list(
                initialize = function(
                  maxMiC = 100, m = 2, theta = 0.8, Theta = 0.8
                ) {
                  maxMiC <<- maxMiC
                  m <<- m
                  theta <<- theta
                  Theta <<- Theta

                  allFMiC <<- list()
                  microCenters <<- data.frame()
                  microWeights <<- numeric()

                  .self
                }
              )
  )

fMicro_refClass$methods(
  cluster = function(x, #weight = NULL,
                     ...){

    # Examples are not weighted

    #if(is.null(weight)){
    #  weights <- rep(1,nrow(x))
    #} else {
    #  if(length(weight)!=nrow(x)) stop("number of weights does not match number of points")
    #  weights <- weight
    #}

    numMiC <- nrow(microCenters)
    nX <- nrow(x)
    minNumFMiC <- 5

    while(nX != 0){

      tryCreate <- FALSE

      p <- x[1, ]
      #w <- weights[1]

      ## if there is at least a few micro cluster
      if(numMiC >= minNumFMiC){
        ## get membership of point to all micro clusters
        dxC <- as.matrix(dist(rbind(microCenters, p))) ### distance matrix for p and all FMiC centers
        dxC <- dxC[numMiC+1,-(numMiC+1)]  ### Only need the one line for distances considering p
        membership <- numeric()

        ## Calculating membership values for p in all availableFMiC
        for(j in 1:numMiC){
          s <- 0
          for(k in 1:numMiC) {
            s <- s + ((dxC[j]/dxC[k])^(2/(m-1)))
          }
          membership <- c(membership, 1/s)
        }

        maxMem <- max(membership)

        ## if maximum membership >= theta
        if(maxMem >= theta){

          ## Update all FMiC where membership is > 0.1
          aux <- membership > 0.1 # vector of booleans indicating FMiC to be updated

          ## add p pondered by membership
          foreach(isUpdated = aux, j = 1:numMiC, mem = membership) %do% {
            if(isUpdated) {
              allFMiC[[j]]$addPoints(x = p, membership = mem)
              microCenters[j,] <<- allFMiC[[j]]$get_center()
              microWeights[j] <<- allFMiC[[j]]$get_weight()
            }
          }

          ## remove p from x
          x <- x[-1,]
          #weights <- weights[-1]
          nX <- nX - 1

        } else {
          tryCreate <- TRUE
        }
      } else {
        tryCreate <- TRUE
      }

      if(tryCreate && numMiC == maxMiC){
        ## if trying to create FMiC, but it's already maxxed out
        ## remove older FMiC
        ts <- sapply(allFMiC, function(fmic) fmic$get_t())

        ## get minimum timestamp and compare with all micro clusters
        minTS <- min(ts)
        maintain <- ts != minTS

        ## maintain is a vector of bool with TRUE for micro clusters to be mantained
        allFMiC <<- allFMiC[maintain] ## remove old micro clusters

        microCenters <<- microCenters[maintain, ] ## remove centers for old micro clusters

        microWeights <<- microWeights[maintain] ## remove weights for old micro clusters

        numMiC <- nrow(microCenters)

      }

      if(tryCreate){
        ## create new micro cluster
        fmic <- FMiC(x = p, membership = 1)
        allFMiC[[numMiC+1]] <<- fmic
        microCenters <<- rbind(microCenters, fmic$get_center())
        microWeights <<- c(microWeights, fmic$get_weight())

        ## So I don't loose attribute names
        names(microCenters) <<- names(p)
        rownames(microCenters) <<- NULL

        numMiC <- nrow(microCenters)

        ## remove p from x
        x <- x[-1,]
        #weights <- weights[-1]
        nX <- nX - 1
      }

      ## merging step
      ## there should be at least minNumFMiC + 1 points to try and merge
      if(numMiC > minNumFMiC+1) {

        ## get distance matrix
        dCC <- as.matrix(dist(microCenters))

        ## calculate membership for minNumFMiC closest points
        aux <- as.numeric(rownames(microCenters))
        mem <-
          foreach(r = aux) %do% {

            distR <- order(dCC[-r,r])

            distR <- head(distR, minNumFMiC)
            distR <- distR + (distR >= r)

            membership <- numeric()

            for(j in distR){
              s <- 0
              for(k in distR) {
                s <- s + ((dCC[j,r]/dCC[k,r])^(2/(m-1)))
              }
              s <- 1/s
              membership <- c(membership, s)
            }

            names(membership) <- distR

            membership
          }

        ## get max membership
        maxMem <- sapply(mem, function(x) max(x))

        ## get mergeable MiC
        merge <- c(1:numMiC)
        merge <- merge[maxMem >= Theta]

        ## if there are MiC to be merged
        if(length(merge) > 0){

          mergedMiC <- list()
          mergedCenters <- numeric()
          mergedWeights <- numeric()
          numMerge <- 0
          indCenters <- c(1:numMiC)
          maintain <- rep(TRUE, times = numMiC)

          while(length(merge) > 0){

            ## get indexes for merging clusters
            c1 <- merge[[1]]
            c2 <- as.numeric(names(mem[[c1]])[1])

            if(maintain[c2]){

              maintain <- maintain & (indCenters != c1 & indCenters != c2)

              fmic1 <- allFMiC[[c1]]
              fmic2 <- allFMiC[[c2]]

              fmic1$merge(fmic2)

              numMerge <- numMerge + 1

              mergedMiC[[numMerge]] <- fmic1
              mergedCenters <- rbind(mergedCenters, fmic1$get_center())
              mergedWeights <- c(mergedWeights, fmic1$get_weight())

              merge <- merge[merge!= c1 & merge != c2]

            } else {
              merge <- merge[merge != c1]
            }
          }

          ## remove merged FCM
          allFMiC <<- allFMiC[maintain] ## remove c1 and c2
          microCenters <<- microCenters[maintain, ] ## remove centers for c1 and c2
          microWeights <<- microWeights[maintain] ## remove weights for c1 and c2

          allFMiC <<- append(allFMiC, mergedMiC)
          microCenters <<- rbind(microCenters, mergedCenters)
          microWeights <<- c(microWeights, mergedWeights)

          numMiC <- nrow(microCenters)

          names(microCenters) <<- names(p)
          rownames(microCenters) <<- NULL

        }

      }
    }
  },

  get_microclusters = function(...) { microCenters },
  get_microweights = function(...) { microWeights }
)
