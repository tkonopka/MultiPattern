## Some temporary functions, include in package or not?
##
##


##' Computes simple Fisher-based enrichment calculations using a set of hits, a set of
##' background items, and concept definitions
##'
##' @param MPconcepts list linking items to concepts
##' @param hits character vector with items
##' @param background character vector with items (background set)
##' 
##' @export
MPenrichment = function(MPconcepts, hits, background) {
    
    allconcepts = unique(unlist(MPconcepts))
    
    ## ensure that background is independent of hits
    nothits = background[!(background %in% hits)]
   
    ##
    ans = rep(0, length(allconcepts))
    names(ans) = allconcepts
    
    ## turn the concepts list around (from item->concepts to concepts->item)
    c2i = setNames(vector("list", length(allconcepts)), allconcepts)
    for (nowitem in names(MPconcepts)) {
        c2i[[nowitem]] = data.frame(Item=nowitem, Concept=MPconcepts[[nowitem]],
               stringsAsFactors=F)
    }
    c2i = data.frame(rbindlist(c2i), stringsAsFactors=F)
    c2i = split(c2i[, "Item"], c2i[, "Concept"])
    
    ## create a matrix with summary stats, for each concept, membershipt in hits and nothits
    temp = matrix(0, ncol=4, nrow=length(allconcepts))
    rownames(temp) = allconcepts
    colnames(temp) = c("hits.1", "hits.0", "bg.1", "bg.0")
    for (nowconcept in allconcepts) {
        hits.1 = sum(hits %in% c2i[[nowconcept]])
        nothits.1 = sum(nothits %in% c2i[[nowconcept]])
        temp[nowconcept, ] = c(hits.1, length(hits)-hits.1, nothits.1, length(nothits)-nothits.1)
    }
    ## apply simple fisher test on the summary stats
    ans = apply(temp, 1, function(x) {fisher.test(matrix(x, nrow=2))$p.value })
    
    return(list(details=temp, fisher=ans))
}





##' Create a new list with concepts
##'
##' @param items character vector
##' 
##' @export
MPnewconcepts = function(items) {
    ans = setNames(vector("list", length(items)), items)
    class(ans) = "MPconcept"
    return(ans)
}



##' print 
##'
##' @param x aa
##' @param ... other arguments, not used
##'
##' @export
print.MPconcept = function(x, ...) {
    cat("MPconcept object:\n")
    cat("  -", length(x), "items\n")
    cat("  -", length(unique(unlist(x))), "concepts\n")
}





##' Manipulates an MPconcepts list
##'
##' @param MPc and MPconcepts object
##' @param item2concept list with concept membership,
##' e.g. list(concept1=c(item1, item2), etc.)
##' @param concept2item list with concept membership,
##' e.g. list(item1=c(concept1, concept2), etc.)
##' @param conceptmatrix matrix or data frame with membership.
##' Rows are trated as items, column names as features/concepts.
##' 
##' @export
MPaddconcepts = function(MPc, item2concept=NULL, concept2item=NULL, conceptmatrix=NULL) {
    
    ## remember how the function was called
    captureMPc = deparse(substitute(MPc))
    
    ## add item2concept information
    if (!is.null(item2concept)) {
        if (class(item2concept)!="list") {
            stop("item2concept must be a list\n")
        }
        for (nowitem in item2concept) {
            MPc[[nowitem]] = c(MPc[[nowitem]], item2concept[[nowitem]])
        }
    }

    ## add from concept2item
    if (!is.null(concept2item)) {
        if (class(concept2item)!="list") {
            stop("concept2item must be a list with names and items\n")
        }        
        ## copy data from the newconcept list into the MPc object
        for (nowconcept in names(concept2item)) {
            nowitems = concept2item[[nowconcept]]
            for (nowitem in nowitems) {
                MPc[[nowitem]] = c(MPc[[nowitem]], nowconcept)
            }
        }
    }
    
    ## add from conceptmatrix
    if (!is.null(conceptmatrix)) {
        if (!class(conceptmatrix) %in% c("matrix", "data.frame")) {
            stop("conceptmatrix must be a matrix or dataframe\n")
        }
        for (nowrow in rownames(conceptmatrix)) {
            nowconcepts = colnames(conceptmatrix)[which(as.logical(conceptmatrix[nowrow,]))]
            if (length(nowconcepts)>0) {
                MPc[[nowrow]] = c(MPc[[nowrow]], nowconcepts)
            }
        }        
    }
    
    ## ensure that duplicate concepts are eliminates
    MPc = lapply(MPc, unique)
    
    ## return the new concepts list silently
    assign(captureMPc, MPc, parent.frame())
    invisible(MPc)          
}



