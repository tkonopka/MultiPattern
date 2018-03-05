## Package: MultiPattern
##
## Functions for making custom plots 
##




##' Custom function to draw mds maps. Uses Rcssplot and some hard settings
##'
##' @param xylayout matrix or list. When matrix, should contain two columns x/y for positions of
##' points on 2d map. When list, should contain several matrices as just described.
##' @param color list with names. Should contain vectors of items in xylayout.
##' These items will be plotted with the css style determined by names(color)
##' @param highlight.points vector of rows in the xylayout that will be
##' plotted with the `highlight` style
##' @param xypadding numeric, determines empty space between points and box
##' @param main character, title for the plot
##' @param legend.separate logical, set TRUE to create a two-panel display. Highlighted
##' configurations will be labeled as A, B, C in main panel, then explained in the second panel.
##' @param squarexy logical, set TRUE to force the scales on x and y axis to be the same
##' (recommended so that euclidean distances on plot correspond to euclidean distances in
##' layout algos)
##' @param label logical, TRUE to print names of points
##' @param Rcss Rcss object
##' @param Rcssclass Rcss class vector
##' @param ... arguments passed on to plot (e.g. xlab, ylab)
##' 
##' @export
MPplotmap = function(xylayout, color=c(),
                     highlight.points=c(), 
                     xypadding=0.05, main="", legend.separate=FALSE,
                     squarexy=TRUE, label=FALSE,
                     Rcss="default", Rcssclass=c(), ...) {
  
  if (!class(xylayout) %in% c("matrix", "data.frame", "list")) {
    stop("xylayout must be a matrix, data.frame, or list\n")
  }
  
  RcssDefaultStyle = RcssGetDefaultStyle(Rcss)
  RcssCompulsoryClass = RcssGetCompulsoryClass(Rcssclass)
  RcssOverload()
  
  ## for ease-of-use, allow input to be a list of layouts
  if (class(xylayout)=="list") {
    for (i in 1:length(xylayout)) {           
      MPplotmap(xylayout[[i]], color=color, xypadding=xypadding,
                main=paste0(main, " - ", names(xylayout)[i]),
                squarexy=squarexy, label=label,
                highlight.points=highlight.points)
    }
    return();       
  }
  
  ## --- Determine which rows are plain, color, and highlighted   
  if (is.null(rownames(xylayout))) {
    rownames(xylayout) = paste0("row:", seq(1, nrow(xylayout)))
  }    
  highlight.rows = highlight.points
  
  normal.rows = rownames(xylayout)
  normal.rows = normal.rows[!(normal.rows %in% c(unlist(color), highlight.rows))]
  ## remove highlight rows from the color list
  color = lapply(color, function(x) { x[!x %in% highlight.rows] })
  
  ## figure out plotting ranges on axes
  xlim = range(xylayout[,1]);
  ylim = range(xylayout[,2]);
  if (squarexy) {
    xylim = MPsquarelim(xylayout[,1], xylayout[,2])
    xlim = xylim$xlim
    ylim = xylim$ylim        
  }
  xlim = xlim+((xlim[2]-xlim[1])*xypadding*c(-1, 1))
  ylim = ylim+((ylim[2]-ylim[1])*xypadding*c(-1, 1))
  
  ## create a two sided layout for the chart
  if (legend.separate) {
    par(mfrow=c(1,2))
  }
  
  ## add elements to the first plot
  plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim,
           xaxs="i", yaxs="i", ...)  
  rect(xlim[1], ylim[1], xlim[2], ylim[2], Rcssclass="box")    
  if (label) {        
    if (length(normal.rows)>0) {
      text(xylayout[normal.rows,1], xylayout[normal.rows,2],
               normal.rows)
    }
    for (nowg in names(color)) {
      color.rows = color[[nowg]]
      if (length(color.rows)>0) {
        text(xylayout[color.rows,1], xylayout[color.rows,2],
                 color.rows, Rcssclass=nowg)
      }
    }        
    if (length(highlight.rows)>0) {
      text(xylayout[highlight.rows,1], xylayout[highlight.rows,2],
               highlight.rows, Rcssclass="highlight")
    }            
  } else {
    if (length(normal.rows)>0) {
      points(xylayout[normal.rows,1], xylayout[normal.rows,2])
    }
    for (nowg in names(color)) {
      color.rows = color[[nowg]]
      if (length(color.rows)>0) {
        points(xylayout[color.rows,1], xylayout[color.rows,2],
                           Rcssclass=nowg)
      }
    }
    ## draw the highlight dots
    if (length(highlight.rows)>0) {
      points(xylayout[highlight.rows, 1],
                 xylayout[highlight.rows, 2],
                 Rcssclass="highlight")
      if (length(highlight.points)>0) {
        nowtext = highlight.points
        if (legend.separate) {
          nowtext = names(highlight.points);                    
        }
        text(xylayout[highlight.points, 1],
             xylayout[highlight.points, 2],
             names(highlight.points),
             Rcssclass="highlight")
      }
    }
  }
  mtext(main, side=3, Rcssclass="main")
  
  ## add elements to a second plot
  if (legend.separate) {
    plot(xlim, ylim, xlim=xlim, ylim=ylim,
         xaxs="i", yaxs="i", type="n", frame=FALSE)
    if (length(highlight.points)>0) {
      nnp = length(highlight.points)
      tempy = seq(ylim[2], ylim[1], length=nnp+2)[2:(nnp+1)] 
      points(rep(xlim[1], nnp), tempy,  Rcssclass="highlight") 
      text(rep(xlim[1], nnp), tempy,
           paste0(nowtext, " - ", highlight.points), pos=4,
           Rcssclass="highlight")            
    }
    if (length(color)>0) {
      nnp = length(color)
      tempy = seq(ylim[2], ylim[1], length=nnp+2)[2:(nnp+1)]
      for (i in seq_len(length(color))) {
        nowg = names(color)[i]
        points(mean(xlim), tempy[i], Rcssclass=nowg)           
        text(mean(xlim), tempy[i], nowg, pos=4, Rcssclass=nowg) 
      }
    }
  }
  
}




##' plot a 2D scatter diagram with clusters in various style
##'
##' @param coords matrix, first two columns will be interpreted as xy coordinates
##' @param clust character, mapping between items and cluster numbers
##' @param Gnames logical, set TRUE to force cluster group names into
##' G1, G2, etc. Set FALSE to use the existing group names in clust.
##' (This only affects Rcss styling)
##' @param main character, displayed as plot title
##' @param Rcss Rcss object
##' @param Rcssclass Rcss class
##' 
##' @export
MPplotScatterWithK = function(coords, clust, Gnames=TRUE, main="",
                              Rcss="default", Rcssclass=c()) {
  
  RcssDefaultStyle = RcssGetDefaultStyle(Rcss)
  RcssCompulsoryClass = RcssGetCompulsoryClass(Rcssclass)
  RcssOverload()
  
  ## use the defined cluster codes to split the points by cluster
  if (Gnames) {
    xycut = split(names(clust), as.integer(as.factor(clust)))
    names(xycut) = paste0("G", as.character(names(xycut)))
  } else {
    xycut = split(names(clust), clust)
  }
  
  ## split up the 2D layout into groups
  xy.coords = list()
  for (i in names(xycut)) {
    xy.coords[[i]] = coords[xycut[[i]],, drop=FALSE]
  }
  
  xylim = MPsquarelim(coords[,1], coords[,2])
  xlim = xylim$xlim
  ylim = xylim$ylim
  
  p.padding = RcssValue("ScatterWithK", "padding", default=0.02)
  xlim = xlim + ((xlim[2]-xlim[1])*p.padding*c(-1,1))
  ylim = ylim + ((ylim[2]-ylim[1])*p.padding*c(-1,1))
  
  ## create a plot area
  plot(xlim, ylim, xaxs="i", yaxs="i", type="n", axes=F, frame=F)
  ## create boxes for main plot and the projections
  rect(xlim[1], ylim[1], xlim[2], ylim[2], Rcssclass="main")
  
  ## add points to the main plot area and the projections
  for (i in names(xy.coords)) {
    nowxy = xy.coords[[i]]
    points(nowxy[,1], nowxy[,2], Rcssclass=i)
  }
  
  mtext(main, side=3)
  invisible(xycut)
}




##' plot a 2D scatter diagram with clusters in various style
##'
##' @param coords matrix, first two columns will be interpreted as xy coordinates
##' @param xyd distance object or matrix
##' @param seed character, sample for which to show neighbors
##' @param seed.links integer, number of links
##' @param main character, displayed as plot title
##' @param Rcss Rcss object
##' @param Rcssclass Rcss class
##' 
##' @export
MPplotScatterWithLinks = function(coords, xyd,
                                  seed = rownames(coords)[1],
                                  seed.links=round(nrow(coords)/4), main="", 
                                  Rcss="default", Rcssclass=c()) {
  
  RcssDefaultStyle = RcssGetDefaultStyle(Rcss)
  RcssCompulsoryClass = RcssGetCompulsoryClass(Rcssclass)
  RcssOverload()
  
  ## cluster and classify objects
  if (class(xyd)=="dist") {
    xyd = as.matrix(xyd)        
  }
  if (!is.null(seed)) {
    xyd = MPrankNeighbors(xyd)
    seed.neighbors = names(sort(xyd[seed,]))
    seed.neighbors = seed.neighbors[seed.neighbors!=seed]
  } else {
    seed.neighbors = c()
  }
  
  xylim = MPsquarelim(coords[,1], coords[,2])
  xlim = xylim$xlim
  ylim = xylim$ylim
  
  p.padding = RcssValue("ScatterWithK", "padding", default=0.02)
  xlim = xlim + ((xlim[2]-xlim[1])*p.padding*c(-1,1))
  ylim = ylim + ((ylim[2]-ylim[1])*p.padding*c(-1,1))
  
  ## create a plot area
  plot(xlim, ylim, xaxs="i", yaxs="i", type="n", axes=F, frame=F)
  ## create boxes for main plot and the projections
  rect(xlim[1], ylim[1], xlim[2], ylim[2], Rcssclass="main")
  
  ## add points to the main plot area and the projections
  notseed = !(rownames(coords) %in% c(seed, seed.neighbors[1:seed.links]))
  if (sum(notseed)>0) {
    points(coords[notseed,1], coords[notseed,2])
  }
  
  val2hex = function(x) {
    format(as.hexmode(round(x*255)), width=2)
  }
  
  ## draw lines from the seed to its nearest neighbors
  if (!is.null(seed)) {
    seedcol = RcssValue("lines", "col", default="#ff0000", Rcssclass="seed")
    for (i in 1:seed.links) {
      nowtrans = val2hex((seed.links-i+1)/(seed.links))
      nowneighbor = seed.neighbors[i]
      lines(coords[c(seed, nowneighbor), 1],
            coords[c(seed, nowneighbor), 2],
            col=paste0(seedcol, nowtrans), Rcssclass="seed")
      points(coords[nowneighbor, 1], coords[nowneighbor, 2],
             bg=paste0(seedcol, nowtrans), Rcssclass="neighbor")
    }
    ## draw the seed sample last (so appears on top)
    points(coords[seed,1], coords[seed,2], Rcssclass="seed")
  }
  
  ## finish up with title 
  mtext(main, side=3)
}


