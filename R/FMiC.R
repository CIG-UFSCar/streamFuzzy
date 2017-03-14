## Def structure for a Fuzzy Micro-Cluster (FMiC), based on TCF (Zhou, 2008)
## FMiC is defined as vector (CF2x, CF1x, t, M)
##  - CF2x: quadratic sum of examples weighted by the membership values
##  - CF1x: is the linear sum of examples weighted by the membership values
##  - t: timestamp for the newest example added
##  - M: scalar cardinality of the micro-cluster (sum of memeberships values)
FMiC <- function(M = 0, CF1x = NULL, CF2x = NULL, t = NULL){

  FMiC <- FMiC$new(M = M, CF1x = CF1x, CF2x = CF2x, t = t)

  desc <- "FMiC structure, based on TCF (Zhou, 2008)"

  structure(list(description = desc, RObj = FMiC), class = c("FMiC"))

}

FMiC <-
  setRefClass("FMiC",
              fields = list(
                ## args
                M = "numeric",
                CF1x = "numeric",
                CF2x = "numeric",
                t = "numeric",
                ## store microcluster
                c = "numeric", # cluster prototype
                N = "numeric", # number of associated examples
                # r = "numeric" # cluster radius
              ),

              methods = list(
                initialize = function(
                  x = NULL, #weight = NULL,
                  membership = NULL
                ) {

                  if(!is.null(x)){

                    ## Examples are weighted by their membership value!

                    # Membership can't be NULL
                    if(is.null(membership)) stop("Membership values can not be NULL")
                    else {
                      # else check if number of membership
                      if(length(membership)!=nrow(x)) stop("number of memberships does not match number of points")
                      # get the membership values passed as args
                      memberships <- membership
                    }

                    M <<- sum(memberships)
                    CF1x <<- mapply(function(x, m) sum(m*(x)), x, memberships)
                    CF2x <<- mapply(function(x, m) sum(m*(x^2)), x, memberships)
                    t <<- as.numeric(Sys.time())

                    N <<- nrow(x)
                    c <<- CF1x/M
                    ## what is the radius of a fuzzy group? :/
                    ## caluclating r as if crisp group
                    #r <<- sqrt((abs(CF2x)/M) - ((abs(CF1x)/M)^2))


                  }
                  else {
                    M <<- numeric()
                    CF1x <<- numeric()
                    CF2x <<- numeric()
                    t <<- numeric()

                    c <<- numeric()
                    N <<- 0
                    #r <<- numeric()
                  }

                  .self
                }
              )
  )

FMiC$methods(
  addPoints = function(x, membership = NULL, ...){

    ## Examples are weighted by their membership value!

    # Membership can't be NULL
    if(is.null(membership)) stop("Membership values can not be NULL")
    else {
      # else check if number of membership
      if(length(membership)!=nrow(x)) stop("number of memberships does not match number of points")
      # get the membership values passed as args
      memberships <- membership
    }

    M <<- M + sum(memberships)
    CF1x <<- CF1x + mapply(function(x, m) sum(m*(x)), x, memberships)
    CF2x <<- CF2x + mapply(function(x, m) sum(m*(x^2)), x, memberships)
    t <<- as.numeric(Sys.time())

    c <<- CF1x/M
    N <<- N + nrow(x)
    #r <<- sqrt((abs(CF2x)/w) - ((abs(CF1x)/w)^2))
  },

  merge = function(fmic, ...){

    M <<- M + fmic$get_M()
    CF1x <<- CF1x + fmic$get_CF1x()
    CF2x <<- CF2x + fmic$get_CF2x()
    t <<- max(t, fmic$get_t())

    c <<- CF1x/M
    N <<- N + fmic$get_N()
    #r <<- sqrt((abs(CF2x)/w) - ((abs(CF1x)/w)^2))

  },

  get_M = function(...) { M },
  get_CF1x = function(...) { CF1x },
  get_CF2x = function(...) { CF2x },
  get_t = function(...) { t },
  get_center = function(...) { c },
  get_N = function(...) {N},
  get_weight = function(...) { M } # M is the weight of the FMiC
  #  get_radius = function(...) { r }
)
