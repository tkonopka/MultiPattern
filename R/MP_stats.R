## functions for statistical tests 
##

##'
##' Compute two p-values from a list of numbers emphasizing contrasts
##'
##' the first p-value is an ranked anova p
##' the second p-value is a high-contrast wilcoxon p
##' for a list with several elements, the wilcoxon groups will be: examples:
##' (1), 2, (3)
##' (1), 2, 3, (4)
##' (1, 2), 3, (4, 5)
##' (1, 2, 3, 4), 5, 6, (7, 8, 9, 10)
##' The groupings can be changed by increasing the contrast, e.g. contrast=2
##' (1), 2, 3, (4)
##' (1), 2, 3, 4, (5)
##' In these examples, the groups witin parenthese are merged together, then
##' the two parentheses are compared using wilcoxon rank sum test
##' 
##' @param li ordered list with numeric values, length must be at least 2
##' @param contrast integer, number of middle groups that will be omitted
##' (when total number of groups is odd, the groups skipped will be contrast+1)
##' 
##' @export
contrast.test = function(li, contrast=1) {
    
    if (class(li)!="list") {
        stop("input must be a list\n")
    }        
    if (length(li)<2) {
        stop("input list must have at least two elements")
    }
    
    ## compute anova from a list
    temp = list()
    for (i in 1:length(li)) {
        if (length(li[[i]])>0) {
            temp[[i]] = data.frame(val=li[[i]], fa=letters[i], stringsAsFactors=F)
        }
    }
    temp = do.call(rbind, temp)
    temp[, "fa"] = as.factor(temp[, "fa"])
    ans.anova = anova(lm(rank(temp[, "val"])~temp[, "fa"]))
    ans.anova = (ans.anova$"Pr(>F)")[1]
    rm(temp)
    
    ## compute wilcoxon from extremes of the list
    halflen = max(1, floor((length(li)-contrast)/2))
    linames = names(li)
    if (is.null(linames)) {
        linames = 1:length(li)
    }
    contrast.groups = list(A=linames[1:halflen], B=rev(linames)[1:halflen])
    contrast.groups$Other = linames[!linames %in% unlist(contrast.groups)]
    
    ans.wilcoxon = list(A=as.numeric(unlist(li[1:halflen])),
        B=as.numeric(unlist(rev(li)[1:halflen])))
    temp=min(sapply(ans.wilcoxon, function(x) {sum(!is.na(x))}))
    if (temp==0) {
        ans.wilcoxon = NA
        ans.diff = NA
    } else {        
        ans.diff = median(ans.wilcoxon$A, na.rm=TRUE)-median(ans.wilcoxon$B, na.rm=TRUE)
        ans.wilcoxon = wilcox.test(ans.wilcoxon$A, ans.wilcoxon$B)
        ans.wilcoxon = ans.wilcoxon$p.value        
    }
    rm(temp)
    
    ans = list(data=li, contrast.groups=contrast.groups,
        anova.p=ans.anova,
        contrast.p=ans.wilcoxon, contrast.diff=ans.diff)
    class(ans) = "contrast"
    
    ans
}




##' pretty print of information from a contrast test
##'
##' @param x object of class contrast
##' @param ... additional parameters (ignored)
##'
##' @export
print.contrast = function(x, ...) {
    
    if (class(x)!="contrast") {
        stop("input must be of class contrast\n")
    }
    
    cat("\nRank-based test on high-contrast groups\n\n")
    cat("data:\t", length(x$data), " groups\n")
    cat("groups:\t", paste0("(",
                           paste(x$contrast.groups$A, collapse=", "), "), ",
                           paste(x$contrast.groups$Other, collapse=", "),
                           ifelse(length(x$contrast.groups$Other)>0, ", ", ""),
                           "(",
                           paste(x$contrast.groups$B, collapse=", "), ")"), "\n\n")
    cat("rank anova p-value:\t", signif(x$anova, 4), "\n")
    cat("contrast p-value:\t", signif(x$contrast.p, 4), "\n\n")
    
}


