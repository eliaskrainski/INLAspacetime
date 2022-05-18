#' Illustrative code for building a mesh in 2d domain.
#'
#' Creates a mesh object.
#'
#' @param loc a two column matrix with location coordinates.
#' @param domain a two column matrix defining the domain.
#' @param max.edge the maximun edge length.
#' @param offset the length of the outer extension.
#' @param SP logical indicating if the output will include the SpatialPolygons.
#' @section Warning:
#'  This is just for illustration purpose and one should consider the
#'  efficient function available a the INLA package.
#' @return a mesh object.
#' @export
mesh2d <-
    function(loc, domain, max.edge, offset, SP=TRUE)
{
### arguments check
    if(missing(loc)) {
        if(missing(domain)) {
            stop("please provide 'loc' or domain'!")
        }
        xy <- domain
    } else {
        if(missing(domain)) {
            xy <- loc
        } else {
            xy <- rbind(loc, domain)
        }
    }
    if(!is.logical(SP))
        stop("'SP' must be logical")
### define retangle around
    xyl <- unique(apply(xy, 2, range, na.rm=TRUE))
    if(!missing(offset)) {
      ro <- which(offset<0)
      if (length(ro)>0)
        offset[ro] <- -(xyl[2, ro] - xyl[1,ro])*offset[ro]
      xyl[1,] <- xyl[1,] - offset
      xyl[2,] <- xyl[2,] + offset
    }
    ### define triangle coordinates
    xyr <- apply(xyl, 2, diff)
    transpose <- (xyr[1]<xyr[2])
    if(transpose) {
      xyl <- xyl[,2:1]
      xyr <- rev(xyr)
    }
    k0 <- trunc(xyr[1]/max.edge) + ((xyr[1]%%max.edge)>0)+1
    edgelen <- xyr[1]/(k0-1)
    h <- edgelen/sqrt(2)
    k <- trunc(xyr[2]/h) + ((xyr[2]%%h)>0)+1
    xu <- seq(xyl[1,1], xyl[2,1], edgelen)
    nxu <- length(xu)
    if(length(xu)!=k0)
      stop(paste0("length(xu)=", nxu, ", k0=", k0, " and k=", k))
    xx <- list(c(xu[1]-edgelen/2, (xu[-nxu]+xu[-1])/2,
                 xu[nxu]+edgelen/2), xu)

    k.a <- 0
    refine <- TRUE
    if (refine)
      xx <- xx[2:1]
    else  {
      while(length(xx[[k.a+2]])>1) {
        k.a <- k.a+1
        xx[[k.a+2]] <- xx[[k.a]][2:length(xx[[k.a+1]])]
      }
      xx <- xx[(k.a+2):1]
    }
    y0 <- ((-(k-1)):k.a)*h + xyl[2,2]
    y0 <- y0 + abs(min(y0)-xyl[1,2])/2
    k <- k + k.a
    if(k>length(xx))
      for(j in (length(xx)+1):k) {
        turn <- refine && ((j-1)>trunc(k/2))
        if (turn) {
          if (length(xx[[j-1]])>length(xx[[j-2]]))
            xj <- xx[[j-2]]
          else
            xj <- xx[[j-2]][-c(1, length(xx[[j-2]]))]
        }
        else
          xj <- c(xx[[j-2]][1]-edgelen, xx[[j-2]],
                  xx[[j-2]][length(xx[[j-2]])]+edgelen)
        xx[[j]] <- xj
      }
    xx <- xx[length(xx):1]
    triang <- list(loc=cbind(unlist(xx),
                             rep(y0, sapply(xx, length))))
    if(transpose)
      triang$loc <- triang$loc[,2:1]
    triang$manifold <- "R2"
    triang$meta$is.refined <- refine
    ### triangle identification
    tv <- list()
    a <- length(xx[[1]])
    i.b <- list()
    if(a==2)
      i.b[[1]] <- cbind(1,2)
    if(a>2)
      i.b[[1]] <- matrix(c(1, rep(2:(a-1), each=2), a),
                         ncol=2, byrow=TRUE)
    a0 <- 0
    for(j in 1:(length(xx)-1)) {
      a <- length(xx[[j]])
      turn <- refine && (a<length(xx[[j+1]]))
      if (turn) {
        tv[[j]] <- rbind(
          cbind(1:(a-1), 2:a, (2:a)+a),
          cbind(1:a, 1:a+a+1, 1:a+a)[a:1,])+a0
        i.b[[j+1]] <- rbind(c(1, a+1), c(a, 2*a+1))+a0
      } else {
        tv[[j]] <- cbind(1:(a-1), 2:a, 1:(a-1)+a)+a0
        if(nrow(tv[[j]])>1) {
          add <- cbind(2:(a-1), 2:(a-1)+a, 1:(a-2)+a)+a0
          tv[[j]] <- rbind(tv[[j]],
                           add[nrow(add):1, ,drop=FALSE])
        }
        i.b[[j+1]] <- rbind(c(1, a+1), c(a, 2*a-1))+a0
      }
      a0 <- a0 + a
    }
    if(refine) {
      a <- length(xx[[length(xx)]])
      if (a>1)
        i.b[[length(i.b)+1]] <- cbind(1:(a-1), 2:a)+a0
    }
    triang$graph <- list(tv=Reduce('rbind', tv))
    triang$n <- nrow(triang$loc)
    triang$segm <- list(bnd=list(loc=NULL,
                                 idx=Reduce('rbind', i.b)))
    triang$segm$bnd$grp <- matrix(0, nrow(triang$segm$bnd$idx), 1)
    triang$segm$bnd$is.bnd <- TRUE
    attr(triang$segm$bnd, 'class') <- 'inla.mesh.segment'
    if(SP) {
      triang$SP <- sp:::SpatialPolygons(
        lapply(1:nrow(triang$graph$tv), function(j) {
          jj <- triang$graph$tv[j, ]
          p <- sp:::Polygon(triang$loc[c(jj, jj[1]), ])
          sp:::Polygons(list(p), paste(j))
        }))
      triang$centroids <- sp:::coordinates(triang$SP)
    }
    triang$loc <- cbind(triang$loc, 0)
    attr(triang, 'class') <- 'inla.mesh'
    return(triang)
}

