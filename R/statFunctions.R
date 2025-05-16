utils::globalVariables(c("Log2FC", "sdata"))


#' @title Differential expression analysis
#'
#' @description This function identifies differentially expressed repertoire levels between groups of samples.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param colGrp a vector of character specifying the column names in the mData slot corresponding to the experimental condition to be analyzed.
#' @param level a character specifying the level of the repertoire to be compared. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param group a vector of character indicating the column name in the mData slot, as well as the two groups to be compared.
#' @details  This function uses the \href{https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeq}{DESeq2} package.
#' Briefly, it estimates the size factors using the poscounts method which deals with zero counts.
#' It then performs a default analysis by estimating the dispersion using a local regression of log dispersions over log base mean.
#' Finally, a generalized linear model is fitted using a Negative Binomial distribution and Wald statistics.
#' @return a data.frame with 6 columns: the repertoire level in rownames, and the baseMean, log2FoldChange, lfcSE, stat, pvalue and padj in columns.
#' The table is ordered by adjusted p-values.
#' @references Hill, M.O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology 54, 427–473.
#' Kindt, R., Van Damme, P., Simons, A.J. (2006). Tree diversity in western Kenya: using profiles to characterise richness and evenness. Biodiversity and Conservation 15, 1253–1270.
#' @export
#' @examples
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#' dds1 <- diffExpGroup(x = RepSeqData,
#'                     colGrp = "cell_subset" ,
#'                     level = "V",
#'                     group = c("cell_subset", "amTreg", "nTreg"))
#'}
#'
diffExpGroup <- function(x, colGrp,
                         level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                         group) {
    if (any(grepl(paste(c("\\+","\\-"), collapse="|"),group[-1]))) stop("Subgroups should not contain (-) or (+).")

    dds <- .toDESeq2(x, colGrp = colGrp, level = level)
    dds <- DESeq2::estimateSizeFactors(dds, type="poscounts")
    dds <- DESeq2::DESeq(dds, fitType = 'local')
    res <- DESeq2::results(dds, contrast = group)
    res <- as.data.frame(res[order(res$padj),])

    return(res)
}

#' @title Differential expression analysis
#'
#' @description Differential expression analysis using DESeq2
#'
#' @details function compares the expression level of a gene segment or a clone using the DESEq2 package.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param colGrp a vector of character, specify the groups to be compared.
#' @param level character, the level of the repertoire to estimate. Should be one of "aaClone", "V", "J", "VJ" or "aaCDR3".
#'
#' @export
#' @keywords internal
#' @return a DESeqDataSet object
#'
.toDESeq2 <- function(x, colGrp, level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3")) {

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The DESeq2 package is required but not installed. Please install it using BiocManager::install('DESeq2').")
  }
  
  if (missing(x)) stop("x is missing. An object of class RepSeqExperiment is epxected.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (missing(colGrp)) stop("user musts provide at least 1 condition.")

  coldat <- mData(x)[, colGrp, drop=FALSE]
  coldat <- apply(coldat, 2, function(x) gsub("\ ", ".", x))
  coldat <- apply(coldat, 2, function(x) gsub("\\+", "p", x))
  coldat <- apply(coldat, 2, function(x) gsub("\\-", "m", x))
  levelChoice <- match.arg(level)
  cts <- countFeatures(x, level=levelChoice)
  if (length(colGrp) > 1) colGrp <- paste(colGrp, collapse="+")
  rownames(coldat) <- gsub("-", ".", rownames(coldat))
  cts <- data.frame(cts, row.names=1, check.names = FALSE)
  colnames(cts) <- gsub("-", ".", colnames(cts))
  cts <- cts[, match(rownames(coldat),colnames(cts))]
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = coldat, design = as.formula(paste0("~" , colGrp)))

  return(dds)
}

# estimate size factor
#
# function
#
# @param x an object of class RepSeqExperiment
# @param level level of repertoire to analyze \code{aaClone, VJ, V, J, aaCDR3}.
# @param method normaliztion method used for size factor computation.
# @param UsePseudoRef a boolean indicating if Chao indices will be computed according to a reference repertoire (geometric mean repertoire across all samples).
# @return a vector of normalized size factors
# @export
# @examples
# estimateSF <- function(x, level=c("aaClone", "aaCDR3"), method=c("Chao", "iChao", "worChao", "Chao.gmmean", "Chao.median"), UsePseudoRef=TRUE) {
#     if (missing(x)) stop("x is missing.")
#     if (!is.RepSeqExperiment(x)) stop("a RepSeqCount object is expected.")
#     chaoest <- function(x) {
#         return(.chao1(x)$chao.est)
#     }
#     gm_mean <- function(x, na.rm=TRUE) {
#         exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
#     }
#     levelChoice <- match.arg(level)
#     meth <- match.arg(method)
#     dat <- data.table::copy(assay(x))
#     chao1 <- ichao <- chaowor <- sj <- s0 <- NULL
#
#     if (meth == "Chao") {
#         chao <- dat[, .(count=sum(count)), by=c("sample_id", levelChoice)][, .(s0=sum(count>0), sj=sum(count), chao1=chaoest(count)), by="sample_id"]
#         if (UsePseudoRef) {
#             ref <- dat[, .(gmmean=gm_mean(count)), by=levelChoice]
#             chao.ref <- chaoest(ref$gmmean)
#             out <- chao[, (chao.ref/s0) * sj]
#             } else out <- chao[, (chao1/s0) * sj]
#         names(out) <- chao$sample_id
#         }
#     if (meth == "Chao.gmmean") {
#         chao <- dat[, .(count=sum(count)), by=c("sample_id", levelChoice)][, .(chao1=chaoest(count)), by="sample_id"]
#         out <- dat[, .(s0=sum(count>0), sj=sum(count)), by="sample_id"][, (gm_mean(chao$chao1)/s0)*sj]
#         names(out) <- tmp$sample_id
#         }
#     if (meth == "Chao.median") {
#         chao <- dat[, .(count=sum(count)), by=c("sample_id", levelChoice)][, .(chao1=chaoest(count)), by="sample_id"]
#         tmp <- dat[, .(s0=sum(count>0), sj=sum(count)), by="sample_id"][, stats::median(chao$chao1[chao$chao1>0], na.rm=TRUE)*sj/s0, by="sample_id"]
#         out <- tmp$V1
#         names(out) <- tmp$sample_id
#         }
#     if (method == "iChao") {
#         chao <- dat[, .(count=sum(count)), by=c("sample_id", levelChoice)][, .(s0=sum(count>0), sj=sum(count), ichao=iChao(count)), by="sample_id"]
#         if (UsePseudoRef) {
#             ref <- dat[, .(gmmean=gm_mean(count)), by=levelChoice]
#             chao.ref <- iChao(ref$gmmean)
#             out <- chao[, (chao.ref/s0) * sj]
#             } else out <- chao[, (ichao/s0) * sj]
#         names(out) <- chao$sample_id
#         }
#     if (method == "worChao") {
#
#         chao <- dat[, .(count=sum(count)), by=c("sample_id", levelChoice)][, .(s0=sum(count>0), sj=sum(count), chaowor=Chaowor(count)), by="sample_id"]
#         if (UsePseudoRef) {
#             ref <- dat[, .(gmmean=gm_mean(count)), by=levelChoice]
#             chao.ref <- Chaowor(ref$gmmean)
#             out <- chao[, (chao.ref/s0) * sj]
#             } else out <- chao[, (chaowor/s0) * sj]
#         names(out) <- chao$sample_id
#         }
#     return(out)
# }


#' @title  Contrast matrix creation
#'
#' @description  create contrast matrix from factor level names
#'
#' @details function creates contrast matrix from factor level names.
#'
#' @param  x a factor
#' @param ... other parameters will be passed to function \code{contr.sum}
#'
#' @return a contrast design matrix
#' @keywords internal
#' @export
named.contr.sum <- function(x, ...) {
    if (is.factor(x)) {
        x <- levels(x)
    } else if (is.numeric(x) & length(x)==1L) {
        stop("cannot create names with integer value. Pass factor levels")
    }
    x <- stats::contr.sum(x, ...)
    colnames(x) <- apply(x,2,function(x)
         paste(names(x[x>0]), names(x[x<0]), sep="-")
    )
    x
}

# Get normalized count
#
# function computes the estimated size factor according to the choice and level of the repertoire, and the
#
# @param x a RepSeqExperiment.
# @param method method used for normalization.
# @param UsePseudoRef a boolen indicatif whether a reference repertoire will be used for normalizaition.
# @return an object of class RepSeqExperiment with normalized counts.
# @export
# @examples
# Discussion https://support.bioconductor.org/p/66067/

# normalizeCounts <- function(x, method=c("Chao", "iChao", "worChao", "Chao.gmmean", "Chao.median"), UsePseudoRef=TRUE) {
#     if (missing(x)) stop("x is missing.")
#     if (!is.RepSeqExperiment(x)) stop("a RepSeqCount object is expected.")
#     choice <- match.arg(method)
#     dat <- data.table::copy(assay(x))
#     sampleinfo <- sData(x)
#     sf <- estimateSF(x, level="aaClone", method=choice, UsePseudoRef=UsePseudoRef)
#     dat[, freq:=count/rep(sf, table(dat$sample_id))]
#     sampleinfo$sf <- sf
#     x.hist <- data.frame(rbind(History(x), history = paste0("normalizedCounts; x=", deparse(substitute(x)), "; method=", choice, "; UsePseudoRef=", UsePseudoRef)), stringsAsFactors = FALSE)
#     out <- methods::new("RepSeqExperiment", assayData=dat, sampleData=sampleinfo, metaData=mData(x), History=x.hist)
#     out@assayData$count<- NULL
#    colnames(out@assayData)[9]<-"count"
#    return(out)
# }

#' @title Calculation of the perturbation score
#'
#' @description This function computes the perturbation scores of the aaCDR3 length distribution within each V gene.
#'
#' Scores are calculated as a distance between each repertoire and the mean of the control group using the ISEApeaks method (Colette and Six., 2002).
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param ctrl.names a vector of characters indicating the sample_ids to be used as controls.
#' @param distance a character specifying the distance method to be used in the calculation of the perturbation scores. Should be one of the following: "manhattan", "euclidean", "canberra", "minkowski" or "maximum".
#' @param p an integer indicating the power of the Minkowski distance. Default is 2.
#' @return a data frame containing the perturbation scores for each V-gene.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' pert <- perturbationScore(RepSeqData,
#'                           ctrl.names = c("tripod-30-813",
#'                                          "tripod-30-815",
#'                                          "tripod-31-846"),
#'                           distance = "manhattan")
#'
perturbationScore <- function(x, ctrl.names, distance = c("manhattan", "euclidean", "canberra", "minkowski" ,"maximum"), p = 2) {
    if (missing(x)) stop("x is missing, an object of class RepSeqExperiment is expected.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    snames <- rownames(mData(x))
    if (missing(ctrl.names)) stop("ctrl.names is missing. A vector of characters containing the names of control samples is expected.")
    if (!any(ctrl.names %in% snames)) stop("all sample names in ctrl.names are not found in sampleData.")
    cts <- data.table::copy(assay(x))
    cts[, aaCDR3.length:=nchar(aaCDR3)]
    spectratype <- cts[, .(count=sum(count)), by=.(sample_id, V, aaCDR3.length)]
    spectratype[,pct:=prop.table(count), by=.(sample_id, V)]
    spectratypew <- data.table::dcast(spectratype, V+aaCDR3.length~sample_id, value.var="pct", fill=0)
    setkey(spectratypew, V, aaCDR3.length)
    spectratypew[, ctrl.mean:=rowMeans(.SD), .SDcols=ctrl.names]
    spectratypew[, ID:=paste0(V, "_", aaCDR3.length)]
    spectratypem <- data.table::melt(spectratypew, id.vars=c("ID", "V", "aaCDR3.length", "ctrl.mean"), variable.name = "sample_id", value.name = "pct")

    d <- tolower(match.arg(distance))
    ctrl.dist <- switch(d,
    manhattan = {
        spectratypem[, .(perturb=sum(abs(pct-ctrl.mean))),by=.(sample_id, V)]
    },
    euclidean = {
        spectratypem[, .(perturb=sqrt(sum((pct-ctrl.mean)^2))),by=.(sample_id, V)]
    },
    maximum = {
        spectratypem[, .(perturb=abs(max(pct-ctrl.mean))),by=.(sample_id, V)]
    },
    minkowski = {
        spectratypem[, .(perturb=(sum((pct-ctrl.mean)^p))^(1/p)),by=.(sample_id, V)]
    },
    canberra = {
        spectratypem[, .(perturb=sum(abs(pct-ctrl.mean)/(abs(pct) + abs(ctrl.mean)))),by=.(sample_id, V)]
    })
    ctrl.distw <- dcast(ctrl.dist, V~sample_id, value.var="perturb")
    out <- data.frame(ctrl.distw, row.names=1, check.names=FALSE)
    return(out)
}


#' @title Differential expression analysis between two samples
#'
#' @description This function compares the expression of a repertoire level between two samples.
#' Log-ratios are calculated on the occurrence of a selected repertoire level between two compared samples, and differentially expressed genes/sequences are identified based on the user-defined log-ratio threshold.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleNames a vector of character with the sample_ids of the repertoires to drop from the RepSeqExperiment object.
#' @param level a character specifying the repertoire level to be analyzed. Should be one of "aaClone","ntClone","aaCDR3","ntCDR3","V", "J",or "VJ".
#' @param scale a character specifying the type of occurrence to take into account: "count" or "frequency".
#' @param th the lof2FC threshold to be used
#' @param remove.zeros a boolean indicating whether or not repertoire levels that are completely absent in one of the two compared samples should be to take into account in the calculation.
#' @return a data frame containing the list of differentially expressed repertoire levels, their occurrence in both samples and their calculated log2FC.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' diff_expression <- diffExpInd(RepSeqData, 
#'                               level="V", 
#'                               scale="frequency", 
#'                               sampleNames = c("tripod-30-813","tripod-30-815"), 
#'                               remove.zeros = FALSE)
#'
diffExpInd <- function(x,  sampleNames = NULL, level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                    scale = c("frequency", "count"), th=1.5, remove.zeros=TRUE) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (length(sampleNames) == 1) stop("Two sampleNames are required.")

  
  if (is.null(sampleNames)) sampleNames <- rownames(sdata)[seq_len(2)]
  levelChoice <- match.arg(level)
  scaleChoice <- match.arg(scale)
  cols <- c("sample_id", levelChoice, "count")
  cts <- data.table::copy(assay(x))
  counts <- cts[sample_id %in% sampleNames, ..cols]
  data2plot <- data.table::dcast(counts, paste(levelChoice, "~sample_id"), value.var="count", fun.aggregate = sum)
  freqs<- function(x){
    frequency= x/sum(x)
    return(frequency)
  }

  if(remove.zeros==TRUE){
    data2plot=data2plot[apply(data2plot!=0, 1, all),]
  }

  if(scaleChoice == "frequency"){
  add_freq<- apply(data2plot[,-1], 2, freqs)
  data2plot<- cbind(data2plot[,1], as.data.frame(add_freq))
  }

  data2plot$Log2FC <- log2(data2plot[,2]/data2plot[,3])
  data2plot<-  data2plot %>%
                  filter(Log2FC >=th | Log2FC <= -th)
  data2plot<- setDT(data2plot)

  return(data2plot)

}
