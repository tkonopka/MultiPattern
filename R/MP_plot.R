


##' Custom function to draw TSNE/mds maps. Uses Rcssplot and some hard settings
##'
##' @param xylayout matrix or list. When matrix, should contain two columns x/y for positions of
##' points on 2d map. When list, should contain several matrices as just described.
##' @param color list with names. Should contain vectors of rownames in xylayout.
##' These items will be plotted with the css style determined by names(color)
##' @param highlight.points vector of rows in the xylayout that will be
##' plotted with the `highlight` style
##' @param highlight.links list with pairs of rows in the xylayout that will be
##' connected by a line with "highlight" style
##' @param xypadding numeric, determines empty space between points and box
##' @param main character, title for the plot
##' @param legend.separate logical, set TRUE to create a two-panel display. Highlighted
##' configurations will be labeled as A, B, C in main panel, then explained in the second panel.
##' @param squarexy logical, set TRUE to force the scales on x and y axis to be the same
##' (recommended so that euclidean distances on plot correspond to euclidean distances in
##' layout algos)
##' @param label logical, TRUE to print names of points
##' @param RC Rcss object
##' @param RCC Rcss class vector
##' @param ... arguments passed on to plot (e.g. xlab, ylab)
##' 
##' @export
MPplotmap = function(xylayout, color=c(),
    highlight.points=c(), highlight.links=list(),
    xypadding=0.05, main="", legend.separate=FALSE,
    squarexy=TRUE, label=FALSE, RC="default", RCC=c(), ...) {
    
    if (!class(xylayout) %in% c("matrix", "data.frame", "list")) {
        stop("xylayout must be a matrix, data.frame, or list\n")
    }
    
    ## for ease-of-use, allow input to be a list of layouts
    if (class(xylayout)=="list") {
        for (i in 1:length(xylayout)) {           
            MPplotmap(xylayout[[i]], color=color, xypadding=xypadding,
                      main=paste0(main, " - ", names(xylayout)[i]),
                      squarexy=squarexy, label=label,
                      highlight.points=highlight.points,
                      highlight.links=highlight.links,
                      RC=RC, RCC=RCC)          
        }
        return();       
    }
    
    ## --- Determine which rows are plain, color, and highlighted   
    if (is.null(rownames(xylayout))) {
        rownames(xylayout) = paste0("row:", seq(1, nrow(xylayout)))
    }    
    highlight.rows = c(c(highlight.points, unlist(highlight.links)))
    
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
        Rcsspar(mfrow=c(1,2), Rcss=RC, Rcssclass=RCC)
    }
    
    ## add elements to the first plot
    Rcssplot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", Rcss=RC, Rcssclass=RCC, ...)
    Rcssrect(xlim[1], ylim[1], xlim[2], ylim[2], Rcss=RC, Rcssclass=c(RCC, "box"))    
    if (label) {        
        if (length(normal.rows)>0) {
            Rcsstext(xylayout[normal.rows,1], xylayout[normal.rows,2], normal.rows,
                     Rcss=RC, Rcssclass=c(RCC))
        }
        for (nowg in names(color)) {
            color.rows = color[[nowg]]
            if (length(color.rows)>0) {
                Rcsstext(xylayout[color.rows,1], xylayout[color.rows,2], color.rows,
                         Rcss=RC, Rcssclass=c(RCC, nowg))
            }
        }        
        if (length(highlight.rows)>0) {
            Rcsstext(xylayout[highlight.rows,1], xylayout[highlight.rows,2], highlight.rows,
                     Rcss=RC, Rcssclass=c(RCC, "highlight"))
        }            
    } else {
        if (length(normal.rows)>0) {
            Rcsspoints(xylayout[normal.rows,1], xylayout[normal.rows,2], Rcss=RC, Rcssclass=c(RCC))
        }
        for (nowg in names(color)) {
            color.rows = color[[nowg]]
            if (length(color.rows)>0) {
                Rcsspoints(xylayout[color.rows,1], xylayout[color.rows,2],
                           Rcss=RC, Rcssclass=c(RCC, nowg))
            }
        }
        ## draw the highlight links
        if (length(highlight.links)>0) {
            for (i in 1:length(highlight.links)) {
                nowpair = highlight.links[[i]][1:2]
                Rcsslines(xylayout[nowpair, 1], xylayout[nowpair, 2],
                          Rcss=RC, Rcssclass=c(RCC, "highlight"))
                nowname = names(highlight.links)[i]
                if (!is.null(nowname)) {
                    Rcsstext(mean(xylayout[nowpair, 1]), mean(xylayout[nowpair, 2]), nowname,
                             Rcss=RC, Rcssclass=c(RCC, "highlight"))
                }
            }
        }
        ## draw the highlight dots
        if (length(highlight.rows)>0) {
            Rcsspoints(xylayout[highlight.rows, 1], xylayout[highlight.rows, 2],
                       Rcss=RC, Rcssclass=c(RCC, "highlight"))
            if (length(highlight.points)>0) {
                nowtext = highlight.points
                if (legend.separate) {
                    nowtext = names(highlight.points);                    
                }
                Rcsstext(xylayout[highlight.points, 1], xylayout[highlight.points, 2],
                         names(highlight.points), Rcss=RC, Rcssclass=c(RCC, "highlight"))
            }
        }
    }
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    
    ## add elements to a second plot
    if (legend.separate) {
        Rcssplot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n",
                 frame=FALSE, Rcss=RC, Rcssclass=RCC)
        if (length(highlight.points)>0) {
            nnp = length(highlight.points)
            tempy = seq(ylim[2], ylim[1], length=nnp+2)[2:(nnp+1)] 
            Rcsspoints(rep(xlim[1], nnp), tempy,  Rcss=RC, Rcssclass=c(RCC, "highlight"))           
            Rcsstext(rep(xlim[1], nnp), tempy, paste0(nowtext, " - ", highlight.points), pos=4,
                     Rcss=RC, Rcssclass=c(RCC, "highlight"))            
        }
        if (length(color)>0) {
            nnp = length(color)
            tempy = seq(ylim[2], ylim[1], length=nnp+2)[2:(nnp+1)]
            for (i in seq_len(length(color))) {
                nowg = names(color)[i]
                Rcsspoints(mean(xlim), tempy[i],  Rcss=RC, Rcssclass=c(RCC, nowg))           
                Rcsstext(mean(xlim), tempy[i], nowg, pos=4, Rcss=RC, Rcssclass=c(RCC, nowg)) 
            }
        }
    }
    
}





## get a matrix indicating the ranks of samples relative to a given seed sample
##
## MPds - list of dissimilarity objects
## samplenames - names of columns/rows in dissimilarity objects
## seedsample - name of sample to consider
## returns a matrix showing how far away all the samples are from the seed
MPgetDistanceFromSeed = function(MPds, samplenames, seedsample) {
    
    ## create a summary object holding distances from the input objects
    ans = matrix(0, ncol=length(samplenames), nrow=length(MPds))
    colnames(ans) = samplenames
    rownames(ans) = names(MPds)

    if (length(MPds)<1) {
        return(ans)
    }
    
    ## fill in the summary table with data from MPds
    for (i in 1:length(MPds)) {
        temp = MPds[[i]]
        if (class(temp)=="dist") {
            temp = as.matrix(temp)
        } else if (class(temp)=="Rle") {
            temp = matrix(as.numeric(temp), ncol=length(samplenames))
            rownames(temp) = colnames(temp) = samplenames
        }
        ans[i, samplenames] = temp[seedsample, samplenames]        
    }
    
    return(ans)    
}





## plot a heatmap of nearest neighbors
MPplotNeighbors = function(MPnei, samplenames, seedsample, consensus=TRUE, 
    Rcss="default", Rcssclass=c()) {
    
    ## vectorize the function on seedsample
    if (length(seedsample)>1) {
        for (nowseed in seedsample) {
            MPplotNeighbors(MPnei, samplenames, nowseed, consensus=consensus,
                            Rcss=Rcss, Rcssclass=Rcssclass)
        }
        return()
    }
    
    ## interpret MPnei either as distance object list, or as a matrix with precomputed
    ## neighbors list
    if (class(MPnei)=="list") {
        MPnei = MPgetDistanceFromSeed(MPnei, samplenames, seedsample)
    }

    RC = Rcss
    RCC = Rcssclass

    val2hex = function(x) {
        format(as.hexmode(round(x*255)), width=2)
    }
    
    ## get some plot info from the Rcss
    nowmai = RcssGetPropertyValueOrDefault(Rcss, "par", "mai", default=c(1,1,1,1),
        Rcssclass=Rcssclass)
    nowcol = RcssGetPropertyValueOrDefault(Rcss, "MPneighbors", "col", default="#ff0000",
        Rcssclass=Rcssclass)
    nowlabspace = RcssGetPropertyValueOrDefault(Rcss, "MPneighbors", "labspace", default=1,
        Rcssclass=Rcssclass)
    
    ## create a table of colors
    MPneicol = matrix("", ncol=ncol(MPnei), nrow=nrow(MPnei))
    MPnei.max = ncol(MPnei)
    MPneicol = paste0(nowcol, val2hex(1-(MPnei/ncol(MPnei))))
    MPneicol = matrix(MPneicol, ncol=ncol(MPnei), nrow=nrow(MPnei))
    colnames(MPneicol) = colnames(MPnei)
    
    xlim = c(0, ncol(MPnei))
    ylim = c(0, nrow(MPnei))
    
    Rcsspar(mai=nowmai, Rcss=RC, Rcssclass=RCC)
    plot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n",
         frame=FALSE, axes=FALSE, ylab="", xlab="")
        
    ## write out the labels on the heatmap
    Rcsstext(rep(-nowlabspace, nrow(MPnei)), seq(ylim[1]+0.5, ylim[2]), rev(rownames(MPnei)),
             Rcss=RC, Rcssclass=c(RCC, "ylab"))
    othersamples = samplenames; ##[!samplenames %in% seedsample]
    Rcsstext(seq(0.5, xlim[2]), rep(ylim[2]+nowlabspace, length(othersamples)), othersamples,
             Rcss=RC, Rcssclass=c(RCC, "xlab"))
    
    ## add color rectangles to the heatmap
    for (i in 1:length(othersamples)) {
        is = othersamples[i]
        Rcssrect(i-1, seq(0, ylim[2]-1), i, seq(1, ylim[2]), col=rev(MPneicol[,is]),
                 Rcss=Rcss, Rcssclass=RCC)
    }
    ## add rectangle around whole thing
    Rcssrect(0, 0, xlim[2], ylim[2], Rcss=RC, Rcssclass=c(RCC, "box"))

    ## add title to the heatmap
    Rcssmtext(seedsample, side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    
    ## create a summary/consensus rank rank
    if (consensus) {
        consensus = apply(MPnei, 2, median)
        consensus.col = paste0(nowcol, val2hex(1-(consensus/ncol(MPnei))))
        names(consensus.col) = colnames(MPnei)
        
        Rcsstext(-nowlabspace, -1.5, "Consensus", Rcss=RC, Rcssclass=c(RCC, "ylab"))
        for (i in 1:length(othersamples)) {
            is = othersamples[i]
            Rcssrect(i-1, -2, i, -1, col=consensus.col[is], Rcss=RC, Rcssclass=RCC)
        }
        Rcssrect(0, -2, xlim[2], -1, Rcss=RC, Rcssclass=c(RCC, "box"))
    }

    ## compute a consensus rank for each of the neighbors
    consensus = median(apply(MPnei, 2, median))    
    return(consensus)
}





##' plot pairwise overlaps summaries
##'
##' e.g. set A overlaps with B by 50%, with C by 30%, with D by 70%, etc.
##'
##' @param ds list of distance matrices
##' @param samplenames vector with samples contained in distance objects
##' @param seedsample one or more sample name (produces one plot for each seed)
##' @param neighborhood integer. Size of neighborhood around seed sample
##' @param ds.names vector with labels to replace ds original names
##' @param bg vector with labels from ds.names that are treated as class 'bg'
##' @param plot logical. Set TRUE to get figures. Set FALSE to get comparison data frame.
##' @param order logical. Set TRUE to arrange multiple seedsample by sum JI
##' @param density.n integer. Passed on to density to compute smoothing fo JI when seedsample=NA
##' @param density.adjust integer
##' @param Rcss Rcss object for Rcssplot
##' @param Rcssclass class for Rcssplot
##'
##' @export
MPplotPairwiseSets = function(ds, samplenames, seedsample, neighborhood=10,
    ds.names=NULL, bg=NULL,
    plot=TRUE, order=TRUE, density.n=128, density.adjust=2, 
    Rcss="default", Rcssclass=c()) {
    
    ## vectorize the function on seedsample
    if (length(seedsample)>1) {
        ## perhaps obtain an ordering of the seed samples
        if (order) {
            seedmetric = list()
            for (nowseed in seedsample) {
                seedmetric[[nowseed]] = MPplotPairwiseSets(ds, samplenames, nowseed,
                              neighborhood=neighborhood, ds.names=ds.names, bg=bg,
                              plot=FALSE, order=FALSE)
            }
            seedmetric = sapply(seedmetric, function(x) {sum(x[,"JI"])})
            seedsample = seedsample[order(seedmetric)]
        }

        ## here, seedsample is either the original seedsample or reordered
        ## to display extreme sumJI values on either end
        for (nowseed in seedsample) {
            MPplotPairwiseSets(ds, samplenames, nowseed, neighborhood=neighborhood,
                               ds.names=ds.names, bg=bg,
                               plot=plot, order=FALSE, 
                               Rcss=Rcss, Rcssclass=Rcssclass)
        }
        return()
    }

    ## helper to compute Jaccard index 
    JI = function(a, b) {
        a.and.b = sum(a%in%b)
        a.or.b = length(unique(c(a, b)))
        return(a.and.b/a.or.b)
    }    
    
    ## ---------------- Compute set comparisons ---------------

    seedsamplemain = seedsample
    if (is.na(seedsample)) {
        seedsample = samplenames
        seedsamplemain = "Average"
    }

    ## create an object holiding JIs for all pairs of similarity objects
    numpairs = length(ds)*(length(ds)-1)/2;
    neicompare = data.frame(Set1=rep("", numpairs), Set2="",
        matrix(0, ncol=length(seedsample), nrow=numpairs), stringsAsFactors=F)
    colnames(neicompare) = c("Set1", "Set2", paste0("JI.", seedsample))
    
    for (nowseed in seedsample) {
        ## get distances from seedsample to all other samples
        MPnei = MPgetDistanceFromSeed(ds, samplenames, nowseed)
        
        ## now MPnei has a matrix indicating the samples nearest to the seed sample
        ## compute sets of sample neighborhoods
        neisets = list()
        MPneicols = samplenames
        for (nowrow in rownames(MPnei)) {
            neisets[[nowrow]] = MPneicols[rank(MPnei[nowrow, ], ties.method="min")<=neighborhood]
        }
        neisets = lapply(neisets, function(x) {x[x!=nowseed]})

        cc = 1
        for (i in 1:(length(neisets)-1)) {
            for (j in (i+1):length(neisets)) {
                neicompare[cc, c("Set1", "Set2")] = names(neisets)[c(i,j)]
                neicompare[cc, paste0("JI.", nowseed)] = JI(neisets[[i]], neisets[[j]])
                cc=cc+1;
            }
        }
        
    }

    ## remove any non-finite values, set to zero
    for (nows in seedsample) {
        temp = paste0("JI.", nows)
        neicompare[!is.finite(neicompare[, temp]), temp] = 0
        rm(temp)
    }
    
    ## order the neicompare rows by JI
    temp = apply(neicompare[, paste0("JI.", seedsample), drop=FALSE], 1, mean)
    neicompare = neicompare[order(temp, decreasing=T), ]
    
    if (!plot) {
        return(neicompare)
    }
    
    ## ---------------- PLOT ---------------
    
    ## get some inputs from the Rcss
    RC = Rcss
    RCC = Rcssclass
    nowmai = RcssGetPropertyValueOrDefault(RC, "par", "mai", default=c(1,1,1,1),
        Rcssclass=RCC)
    nowcol = RcssGetPropertyValueOrDefault(RC, "MPneighbors", "col", default="#ff0000",
        Rcssclass=RCC)
    nowlabspace = RcssGetPropertyValueOrDefault(RC, "MPneighbors", "labspace", default=1,
        Rcssclass=RCC)
    nowcellsize = RcssGetPropertyValueOrDefault(RC, "MPneighbors", "cellsize", default=c(0.05,0.05),
        Rcssclass=RCC)    
    barlength = RcssGetPropertyValueOrDefault(RC, "MPneighbors", "barlength", default=6, 
        Rcssclass=RCC)
    barwidth = RcssGetPropertyValueOrDefault(RC, "MPneighbors", "barwidth", default=0.8, 
        Rcssclass=RCC)
    barwidth2 = barwidth/2
    ## compute the effective margin
    effmai = nowmai
    effmai[1] = nowmai[1]+(nowcellsize[2]*nrow(MPnei))
    ## compute the effective class
    neiclass = rep("good", nrow(neicompare))
    neiclass[neicompare[, "Set1"] %in% bg | neicompare[, "Set2"] %in% bg] = "bg"

    
    ## draw the bar chart with Jaccard Index
    xlim = c(-0.5+barwidth2, nrow(neicompare))
    ylim = c(0, 1)
    allsets = rownames(MPnei)
    vscale = (ylim[2]-ylim[1])/barlength;
    
    Rcsspar(mai=effmai, Rcss=RC, Rcssclass=RCC)
    plot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n",
         frame=FALSE, axes=FALSE, ylab="", xlab="")
    
    ## ---------------------------------------
    
    ## plot the individual bars
    if (ncol(neicompare)==3) {
        for (i in 1:nrow(neicompare)) {
            Rcssrect(i-0.5-barwidth2, 0, i-0.5+barwidth2, neicompare[i, paste0("JI.", seedsample)],
                     Rcss=RC, Rcssclass=c(RCC, "barplot", neiclass[i]))
        }
    } else {
        for (i in 1:nrow(neicompare)) {
            temp = as.numeric(neicompare[i, 3:ncol(neicompare)])
            tempdensity = density(temp, n=density.n, adjust=density.adjust)
            tempmax = max(tempdensity$y)            
            tempx = c(tempdensity$x, rev(tempdensity$x))
            tempx[tempx<0]=0
            tempx[tempx>1]=1
            tempy = c(tempdensity$y/tempmax, -rev(tempdensity$y/tempmax))
            Rcsspolygon(i-0.5+(tempy*barwidth2), tempx, Rcss=RC, Rcssclass=c(RCC, neiclass[i]))
            rm(tempx, tempy, tempmax, temp, tempdensity)
        }
    }
    ## plot the axis for the bar chart
    Rcssaxis(side=2, labels=NA, Rcss=RC, Rcssclass=RCC)
    Rcssaxis(side=2, lwd=0, Rcss=RC, Rcssclass=c(RCC, "ylab"))
    Rcssaxis(side=1, at=xlim, labels=NA, tck=0, Rcss=RC, Rcssclass=RCC)
    
    ## ---------------------------------------    
    ## write out the labels on the setmap
    allsets = rownames(MPnei)
    vscale = (ylim[2]-ylim[1])/barlength;
    setv = (seq(-0.5, -length(allsets))-nowlabspace)*vscale
    names(setv) = allsets
    if (!is.null(ds.names)) {
        rownames(MPnei) = ds.names
    }
    Rcsstext(rep(-nowlabspace, nrow(MPnei)), setv, rownames(MPnei),
             Rcss=RC, Rcssclass=c(RCC, "ylab"))
    
    ## add connectors on the setmap
    for (i in 1:nrow(neicompare)) {
        nsets = as.character(neicompare[i, c("Set1", "Set2")])
        nbg = allsets[!(allsets %in% nsets)]
        Rcsspoints(rep(i-0.5, length(nbg)), setv[nbg], Rcss=RC, Rcssclass=c(RCC, "bg"))
        Rcsslines(rep(i-0.5, 2), setv[nsets], Rcss=RC, Rcssclass=c(RCC, "highlight", neiclass[i]))
        Rcsspoints(rep(i-0.5, 2), setv[nsets], Rcss=RC, Rcssclass=c(RCC, "highlight", neiclass[i]))
    }
    
    ## add title to the heatmap
    Rcssmtext(paste0(seedsamplemain, " - ", neighborhood , " neighbors"),
              side=3, Rcss=RC, Rcssclass=c(RCC, "main"))
    Rcssmtext("Jaccard Index", side=2, Rcss=RC, Rcssclass=c(RCC, "ylab"))
    
}






##' plot a 2D scatter diagram with clusters in various style
##'
##' @param coords matrix, first two columns will be interpreted as xy coordinates
##' @param xyd distance object or matrix
##' @param clust.method character, set as either method for hclust, or "pam", or "none". When "none",
##' the xyd object is interpreted as a list with items
##' @param clust.k integer, number of clusters to color. Clusters will be called G1, G2, G3, etc.
##' Use these codes to tune appearance of these classes in the Rcss. 
##' @param main character, displayed as plot title
##' @param RC Rcss object
##' @param RCC Rcss class
##' 
##' 
MPplotScatterWithK.old = function(coords, xyd, clust.method="single", clust.k=2, 
    main="", 
    RC="default", RCC=c()) {
    
    ## cluster and classify objects
    if (class(xyd)=="matrix") {
        xyd = as.dist(xyd)        
    }
    if (clust.method=="none") {
        xycut = xyd
    } else {
        if (clust.method=="pam") {
            xycut = pam(xyd, k=clust.k, cluster.only=TRUE)
        } else {
            xycut = MPcutree(hclust(xyd, method=clust.method), k=clust.k)        
        }
        xycut = split(names(xycut), xycut)        
        names(xycut) = paste0("G", as.character(names(xycut)))
    }
    
    ## split up the 2D layout into groups
    xy.coords = list()
    for (i in names(xycut)) {
        xy.coords[[i]] = coords[xycut[[i]],, drop=FALSE]
    }
    
    xylim = MPsquarelim(coords[,1], coords[,2])
    xlim = xylim$xlim
    ylim = xylim$ylim
    
    p.padding = RcssGetPropertyValueOrDefault(RC, "ScatterWithK", "padding",
        default=0.02, Rcssclass=RCC)
    xlim = xlim + ((xlim[2]-xlim[1])*p.padding*c(-1,1))
    ylim = ylim + ((ylim[2]-ylim[1])*p.padding*c(-1,1))
    
    ## create a plot area
    Rcssplot(xlim, ylim, xaxs="i", yaxs="i", type="n", Rcss=RC, Rcssclass=RCC, axes=F, frame=F)
    ## create boxes for main plot and the projections
    Rcssrect(xlim[1], ylim[1], xlim[2], ylim[2], Rcss=RC, Rcssclass=c(RCC, "main"))
    
    ## add points to the main plot area and the projections
    for (i in names(xy.coords)) {
        nowxy = xy.coords[[i]]
        Rcsspoints(nowxy[,1], nowxy[,2], Rcss=RC, Rcssclass=c(RCC, i))
    }
    
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=RCC)
    invisible(xycut)
}



##' plot a 2D scatter diagram with clusters in various style
##'
##' @param coords matrix, first two columns will be interpreted as xy coordinates
##' @param clust character, mapping between items and cluster numbers
##' @param Gnames logical, set TRUE to force cluster group names into G1, G2, etc.
##' set FALSE to use the existing group names in clust
##' @param main character, displayed as plot title
##' @param RC Rcss object
##' @param RCC Rcss class
##' 
##' @export
MPplotScatterWithK = function(coords, clust, Gnames=TRUE, main="", RC="default", RCC=c()) {
    
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
    
    p.padding = RcssGetPropertyValueOrDefault(RC, "ScatterWithK", "padding",
        default=0.02, Rcssclass=RCC)
    xlim = xlim + ((xlim[2]-xlim[1])*p.padding*c(-1,1))
    ylim = ylim + ((ylim[2]-ylim[1])*p.padding*c(-1,1))
    
    ## create a plot area
    Rcssplot(xlim, ylim, xaxs="i", yaxs="i", type="n", Rcss=RC, Rcssclass=RCC, axes=F, frame=F)
    ## create boxes for main plot and the projections
    Rcssrect(xlim[1], ylim[1], xlim[2], ylim[2], Rcss=RC, Rcssclass=c(RCC, "main"))
    
    ## add points to the main plot area and the projections
    for (i in names(xy.coords)) {
        nowxy = xy.coords[[i]]
        Rcsspoints(nowxy[,1], nowxy[,2], Rcss=RC, Rcssclass=c(RCC, i))
    }
    
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=RCC)
    invisible(xycut)
}






##' plot a 2D scatter diagram with clusters in various style
##'
##' @param coords matrix, first two columns will be interpreted as xy coordinates
##' @param xyd distance object or matrix
##' @param seed character, sample for which to show neighbors
##' @param seed.links integer, number of links
##' @param main character, displayed as plot title
##' @param RC Rcss object
##' @param RCC Rcss class
##' 
##' @export
MPplotScatterWithLinks = function(coords, xyd,
    seed = rownames(coords)[1], seed.links=round(nrow(coords)/4),
    main="", 
    RC="default", RCC=c()) {
    
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
    
    p.padding = RcssGetPropertyValueOrDefault(RC, "ScatterWithK", "padding",
        default=0.02, Rcssclass=RCC)
    xlim = xlim + ((xlim[2]-xlim[1])*p.padding*c(-1,1))
    ylim = ylim + ((ylim[2]-ylim[1])*p.padding*c(-1,1))
    
    ## create a plot area
    Rcssplot(xlim, ylim, xaxs="i", yaxs="i", type="n", Rcss=RC, Rcssclass=RCC, axes=F, frame=F)
    ## create boxes for main plot and the projections
    Rcssrect(xlim[1], ylim[1], xlim[2], ylim[2], Rcss=RC, Rcssclass=c(RCC, "main"))
    
    ## add points to the main plot area and the projections
    notseed = !(rownames(coords) %in% c(seed, seed.neighbors[1:seed.links]))
    if (sum(notseed)>0) {
        Rcsspoints(coords[notseed,1], coords[notseed,2], Rcss=RC, Rcssclass=RCC)
    }
    
    val2hex = function(x) {
        format(as.hexmode(round(x*255)), width=2)
    }
    
    ## draw lines from the seed to its nearest neighbors
    if (!is.null(seed)) {
        seedcol = RcssGetPropertyValueOrDefault(RC, "lines", "col",
            default="#ff0000", Rcssclass=c(RCC, "seed"))
        for (i in 1:seed.links) {
            nowtrans = val2hex((seed.links-i+1)/(seed.links))
            nowneighbor = seed.neighbors[i]
            Rcsslines(coords[c(seed, nowneighbor), 1],
                      coords[c(seed, nowneighbor), 2],
                      col=paste0(seedcol, nowtrans), Rcss=RC, Rcssclass=c(RCC, "seed"))
            Rcsspoints(coords[nowneighbor, 1], coords[nowneighbor, 2],
                       bg=paste0(seedcol, nowtrans), Rcss=RC, Rcssclass=c(RCC, "neighbor"))
        }
        ## draw the seed sample last (so appears on top)
        Rcsspoints(coords[seed,1], coords[seed,2], Rcss=RC, Rcssclass=c(RCC, "seed"))
    }
    
    ## finish up with title 
    Rcssmtext(main, side=3, Rcss=RC, Rcssclass=RCC)
}




## #######################################################################################
## Plot functions borrowed from package ExpCube (Tomasz Konopka)



##' Draw a vertical recrangle with a color scale and labels
##'
##' This function is an (edited) copy of E3drawScaleLegend from ExpCube (Tomasz Konopka)
##' 
##' @param legend - named vector of colors. Colors defining a color scale. All
##' names associated with vector will be draw onto to the plot. (Set some names to ""
##' to avoid labeling every single color in a color gradient)
##' @param xlim - numeric vector of two elements. Gives range of x axis.
##' @param ylim - numeric vector of two elements. Gives range of y axis.
##' @param main - text two write above the legend
##' @param RC - Rcss object. Style to use for plotting, uses package Rcssplot
##' @param RCC - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
MPdrawScaleLegend = function(legend=NULL, xlim=c(0,1), ylim=c(0,1), main="",   
    RC="default", RCC=c()) {
    
    ## get features of the legend from the RC
    relposx =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relposx",
        default=-0.1, Rcssclass=RCC)
    relwidth =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relwidth",
        default=0.05, Rcssclass=RCC)
    relposy =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relposy",
        default=0.75, Rcssclass=RCC)    
    relheight =  RcssGetPropertyValueOrDefault(RC, "scalelegend", "relheight",
        default=0.15, Rcssclass=RCC)
    
    ## find distances on x/y axes
    xd = xlim[2]-xlim[1]
    yd = ylim[2]-ylim[1]
    ## define coordinates of legend corners
    bl = c(xlim[1]+(xd*relposx), ylim[1]+(yd*relposy))
    br = c(bl[1]+(xd*relwidth), bl[2])
    tl = c(bl[1], bl[2]+(yd*relheight))
    tr = c(br[1], tl[2])
    
    ysteps = seq(bl[2], tl[2], length=length(legend)+1)
    absv = yd*RcssGetPropertyValueOrDefault(RC, "scalelegend", "relv",
        default=0.05, Rcssclass=RCC)
    
    ## draw a legend in the corner
    if (!is.null(legend) & length(legend)>0) {        
        for (j in 1:length(legend)) {
            Rcssrect(bl[1], ysteps[j], tr[1], ysteps[j+1], col=legend[j], Rcss=RC,
                     Rcssclass=c(RCC,"scalelegend"))
        }
        Rcsslines(c(bl[1], br[1], tr[1], tl[1], bl[1]), c(bl[2], br[2], tr[2], tl[2], bl[2]), 
                  Rcss=RC, Rcssclass=c(RCC, "scalelegend"))
        Rcsstext(br[1], seq(bl[2], tl[2], length=length(legend)), names(legend),
                 Rcss=RC, Rcssclass=c(RCC, "scalelegend"))

        Rcsstext(tl[1], tl[2]+absv, main,
                 Rcss=RC, Rcssclass=c(RCC, "scalelegend", "main"))        
    }
    
    
}

##' Plot a basic heat map
##'
##' This function is an (edited) copy of E3PlotBasicHeat from ExpCube (Tomasz Konopka)
##'
##' @param dat - matrix of colors.
##' @param legend - vector of colors. Legend shading. Warning: the user bares responsibility
##' to make sure the legend matches with the dat matrix. 
##' @param xlabels - logical. Toggle display of labels on x axis.
##' @param ylabels - logical. Toggle display of labels on y axis.
##' @param dividers - logical. Toggle vertical dotted lines between x items.
##' @param xlab - character string. Text to display below x axis
##' @param ylab - character string. Text to display below y axis.
##' @param main - character string. Text to display as title, above heatmap.
##' @param Rcss - Rcss object. Used to style the heatmap with Rcssplot.
##' @param Rcssclass - character vector. Classes to tune Rcssplot formatting.
##' 
##' @export
MPplotBasicHeat = function (dat, legend = NULL, ylabels = FALSE, xlabels = TRUE, 
    dividers = TRUE, xlab = "", ylab = "", main = "", Rcss = "default", 
    Rcssclass = c()) {
    
    RC = Rcss
    RCC = Rcssclass
    xlim = c(0, ncol(dat))
    ylim = c(0, nrow(dat))
    xdown = RcssGetPropertyValueOrDefault(RC, "basicheat", "xdown", 
        default = -0.1, Rcssclass = RCC)
    absdown = RcssGetPropertyValueOrDefault(RC, "basicheat", 
        "absdown", default = 0, Rcssclass = RCC)
    Rcsspar(Rcss = RC, Rcssclass = RCC)
    plot(xlim, ylim, type = "n", xlim = xlim, ylim = ylim, xaxs = "i", 
        yaxs = "i", xlab = "", ylab = "", frame = F, axes = F)
    for (i in 1:ncol(dat)) {
        jj = 1:nrow(dat)
        ii = rep(i, nrow(dat))
        Rcssrect(ii - 1, jj - 1, ii, jj, col = dat[, i], Rcss = RC, 
            Rcssclass = RCC)
        dat.na = is.na(dat[, i])
        if (sum(dat.na) > 0) {
            Rcssrect((ii - 1)[dat.na], (jj - 1)[dat.na], ii[dat.na], 
                jj[dat.na], Rcss = RC, Rcssclass = c(RCC, "NA"))
        }
    }
    if (xlabels) {
        if (absdown > 0) {
            xlabpos = xdown
        }
        else {
            xlabpos = xdown * ylim[2]
        }
        Rcsstext(seq(1, ncol(dat)) - 0.5, xlabpos, colnames(dat), 
            Rcss = RC, Rcssclass = c(RCC, "x"))
    }
    if (ylabels) {
        Rcsstext(0, seq(1, nrow(dat)) - 0.5, rownames(dat), pos = 2, 
            Rcss = RC, Rcssclass = c(RCC, "y"))
    }
    if (dividers) {
        for (i in 1:(ncol(dat) - 1)) {
            Rcsslines(rep(i, 2), ylim, Rcss = RC, Rcssclass = c(RCC, 
                "divider"))
        }
    }
    Rcsslines(c(xlim, rev(xlim), xlim[1]), c(rep(ylim, each = 2), 
        ylim[1]), Rcss = RC, Rcssclass = c(RCC, "box"))
    MPdrawScaleLegend(legend = legend, xlim = xlim, ylim = ylim, 
        RC = RC, RCC = RCC)
    Rcssmtext(main, side = 3, Rcss = RC, Rcssclass = c(RCC, "main"))
    Rcssmtext(ylab, side = 2, Rcss = RC, Rcssclass = c(RCC, "y"))
    Rcssmtext(xlab, side = 1, Rcss = RC, Rcssclass = c(RCC, "x"))
}





##' Convert a matrix of values (range [-Inf,Inf] to colors using transparency)
##'
##' This funciton is a copyt of E3Val2Col from package ExpCube (Tomasz Konopka)
##' 
##' @param x - numeric matrix.
##' @param col - vector of two colors in #XXXXXX format. First element determines
##' color associated with negative values. Second element determines color associated
##' with positive values
##' @param maxval - numeric. Value for which color reaches saturation
##'
##' @export
MPvalMat2ColMat = function(x, col=c("#0000ff", "#ff0000"), maxval=1) {

    ## helper function to convert a number [0,1] into a hex transparency
    val2hex = function (x) {
        ans = as.character(as.hexmode(round(x * 255)))
        shortans = nchar(ans) < 2
        if (sum(shortans, na.rm=TRUE) > 0) {
            ans[shortans] = paste("0", ans[shortans], sep = "")
        }
        ans
    }

    ## keep track of all NAs (will restore at the end)
    x.na = is.na(x)
    
    xpos = x>0 & !x.na;
    xneg = x<0 & !x.na;        
    ## force absolute values into range [0,1]
    y = abs(x);
    y[y>maxval]=maxval;
    yvals = as.vector(y)/maxval;
    ## make a matrix of colors
    temp = paste0(col[1], val2hex(yvals));
    if (sum(xpos, na.rm=TRUE)>0) {
        temp[xpos] = paste0(col[2], val2hex(yvals[xpos]));
    }
    ans = matrix(temp, ncol=ncol(x), nrow=nrow(x));
    rownames(ans) = rownames(x);
    colnames(ans) = colnames(x);
    ## replace full transparency with #ffffff
    ans[x==0] = "#ffffff";
    ans[x.na] = NA
    ans
}




##' Convert a number into a color using  transparency
##'
##' This is an (edited) copy of E3Val2Col (Tomasz Konopka)
##' 
##' @param x - numeric matrix.
##' @param col - vector of two colors in #XXXXXX format. First element determines
##' color associated with negative values. Second element determines color associated
##' with positive values
##' @param maxval - numeric. Value for which color reaches saturation
##'
##' @export
MPval2Col = function(x, col=c("#0000ff", "#ff0000"), maxval=1) {

    ## helper function to convert a number [0,1] into a hex transparency
    val2hex = function (x) {
        ans = as.character(as.hexmode(round(x * 255)))
        shortans = nchar(ans) < 2
        if (sum(shortans, na.rm=TRUE) > 0) {
            ans[shortans] = paste("0", ans[shortans], sep = "")
        }
        ans
    }    

    ## keep track of all NAs (will restore at the end)
    x.na = is.na(x)
    
    ## create a boolean vector of colors based on sign of x
    ans = rep(col[1], length(x))
    ans[x>0 & !x.na] = col[2]    
    ## make sure values are within [-1, 1] range
    x[x>maxval & !x.na] = maxval;
    x[x<(-maxval) & !x.na] = -maxval;
    x = x/maxval    
    ## append a transparency value and that's it
    ans = paste0(ans, val2hex(abs(x)))
    names(ans) = names(x)
    ans[x.na] = NA
    ans
}

