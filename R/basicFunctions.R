utils::globalVariables(c("J", ".", "..sNames", "sdata", "metaData", ".SD", 
                        ".SDcols", "key", ".N", "count", "..keep.cols","ntCDR3",
                        "sampleNames", "aaCDR3.length", "aaCDR3", "ntClone",
                        "pct", "ctrl.mean", "ID","prop", "nSequences", "ranks",
                        "cumulative_frequency"))


#' @title Computing the occurrence of any repertoire level
#'
#' @description  This function calculates the occurrence of a selected
#' repertoire level and returns the calculated values for all the samples
#' within a RepSeqExperiment object.
#' It takes into account the weight of the studied level, i.e. the number of
#' sequences expressing a certain gene segment, or the count of a sequence in
#' a sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the repertoire level to be analyzed.
#' Should be one of "aaClone","ntClone","aaCDR3","ntCDR3","V", "J",or "VJ".
#' @param scale a character specifying the type of occurrence to return:
#' "count" or "frequency".
#' @param group a vector of character indicating the group column name in
#' mData and one experimental group within this column. Samples belonging to
#' the chosen group will be analyzed. The column must be of class factor.
#' Default is NULL, values are calculated in all the samples within the dataset.
#'
#' @return a data.table summarizing the count or frequency of the analyzed
#' level. In this table, rows correspond to the repertoire level, and columns
#' correspond to the sample_ids.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' level_statistics <- countFeatures(x = RepSeqData,
#'                                   level = "V",
#'                                   group=c("cell_subset", "amTreg"),
#'                                   scale="frequency")
#'
#' level_statistics <- countFeatures(x = RepSeqData,
#'                                   level = "J",
#'                                   group=c("sex", "F"),
#'                                   scale="count")
#'
#'
countFeatures <- function(x,
                          level = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "ntClone", "aaClone"),
                          scale = c("count", "frequency"), group = NULL) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment 
                                      is expected.")
  
  metaData <- mData(x)
  if(!is.null(group) & nrow(metaData[metaData[, group[1]] == group[2],])==0) stop("The chosen group is not present in the metaData")
  
  levelChoice <- match.arg(level)
  cts <- data.table::copy(assay(x))
  scl <- match.arg(scale)
  
  if (!is.null(group)) {
    grp <- metaData[, group[1]]
    grp.name <- group[2]
    sampleNames <- rownames(metaData[grp %in% grp.name, ])
    cts <- cts[sampleNames, on = "sample_id"]
  }
  if (scl == "count"){
    out <- data.table::dcast(cts, as.formula(paste0(levelChoice,"~ sample_id")), value.var="count", fun=sum)
  } else {
    cts <- cts[, .(count = sum(count)), by=c("sample_id", levelChoice)][, prop := count/sum(count), by = "sample_id"]
    out <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~sample_id")), value.var = "prop", fill = 0)
  }
    return(out)
}

# @title Gene and clone usage
#
# @description compute gene segment or clone usage
#
# @details function computes segment or clone frequency in each sample. The function takes into account the occurence of the analyzed level.
#
# @param x an object of class \code{\linkS4class{RepSeqExperiment}}
# @param level character, the level of the repertoire to estimate. Should be one of "aaClone", "aaCDR3","V", "J" or "VJ"
# @return a data.table


# segmentUsage <- function(x, level=c("aaClone", "V", "J", "VJ","aaCDR3")) {
#     if (missing(x)) stop("x is missing.")
#     if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#     levelChoice <- match.arg(level)
#     cts <- data.table::copy(assay(x))
#     cts <- cts[, .(count = sum(count)), by=c("sample_id", levelChoice)][, prop := count/sum(count), by = "sample_id"]
#     prop <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~sample_id")), value.var = "prop", fill = 0)
#     return(prop)
# }

#' @title Filtering based on sequence count
#'
#' @description This function filters out, in each sample, sequences having a
#' count equal to or below a specified threshold. It can be applied on all the
#' samples within a RepSeqExperiment object, or a group of samples belonging to
#' a specific experimental group.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the type of sequences on which the
#' filtering will be applied. Should be one of "aaClone","ntClone","aaCDR3" or
#' "ntCDR3". For instance, for level="aaCDR3", counts will first be recalculated
#' based on this column. Then, aaCDR3 sequences with counts equal or below the
#' n threshold will be excluded.
#' @param n an integer specifying the count threshold below which sequences will
#' be filtered out. For instance, for n=2, sequences with a count of 1 and 2
#' are filtered out.
#' @param group a vector of character indicating the group column name in mData
#' and one experimental group within this column. Samples belonging to the
#' chosen group will be analyzed. The column must be of class factor. Default
#' is NULL, values are calculated in all the samples within the dataset.

#' @return an object of class \code{RepSeqExperiment}
#' @export
#' @examples
#'
#' filterdata <- filterCount(x = RepSeqData,
#'                           n = 1,
#'                           level = "aaCDR3",
#'                           group = c("sex", "M"))
#'
#'
filterCount <- function(x, level=c("aaClone","ntClone","aaCDR3","ntCDR3"), n=1, group=NULL) {
    V1 <- NULL
    chao1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if (!is(n, 'numeric')) stop("n should be an integer")
    
    
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    metaData <- mData(x)
    if(!is.null(group) & nrow(metaData[metaData[, group[1]] == group[2],])==0) stop("The chosen group is not present in the metaData")
    

    if (!is.null(group)) {
      grp <- metaData[, group[1]]
      grp.name <- group[2]
      sampleNames <- rownames(metaData[grp %in% grp.name, ])
      keep <- cts[, sum(count) <= n, by =  c(levelChoice,"sample_id")]
      keep <- keep[ V1==FALSE | !sample_id %in% sampleNames,]
      res <- cts[keep, on = c(levelChoice,"sample_id")][, V1:=NULL]
      setkey(res, sample_id)
      filterout<- cts[, sum(count) <= n, by =  c(levelChoice,"sample_id")][V1 == TRUE & sample_id %in% sampleNames, ]
    }  else {
    keep <- cts[, sum(count) <= n, by =  c(levelChoice,"sample_id")][V1 == FALSE, ]
    res <- cts[keep, on = c(levelChoice,"sample_id")][, V1:=NULL]
    setkey(res, sample_id)
    filterout<- cts[, sum(count) <= n, by =  c(levelChoice,"sample_id")][V1 == TRUE, ]
    }

    nfilter <- nrow(cts) - nrow(res)

    stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
    
    pastek <- function(y) list( chao1 = .chao1(y)$chao.est, iChao = iChao(y))
    out <- data.frame(res[, .(count = sum(count)), by = c("sample_id", "ntClone")][, pastek(count), by = "sample_id"], row.names = 1)
    stats <- data.frame(base::merge(stats, out, by = 0, sort=FALSE), row.names = 1, stringsAsFactors = TRUE)

    metaData<- metaData %>%
      dplyr::select(-c(nSequences,ntCDR3,aaCDR3,V,J,VJ,aaClone,ntClone, chao1, iChao))
    sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                        row.names=1,  stringsAsFactors = TRUE)
    
    sdata <- sdata[order(match(rownames(sdata), rownames(stats))), ]

    out <- new("RepSeqExperiment",
                    assayData = res,
                    metaData = sdata,
                    otherData = oData(x),
                    History = data.frame(rbind(History(x),
                                        data.frame(history = paste(nfilter,levelChoice,"were filtered using filterCount: n=",n, "and group=",group))),
                                        stringsAsFactors=FALSE))
    oData(out) <- c(oData(out), filterCount=list(  cts[filterout, on = c(levelChoice,"sample_id")][, V1:=NULL]))

    rm(cts, keep)

    return(out)
}


#' @title Extract private sequences
#'
#' @description This function extract private sequences from the whole dataset,
#'  i.e. sequences found exclusively in one sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken
#' into account. Should be one of aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param singletons a boolean indicating whether or not private sequences with
#' a count of 1 should be extracted. Default is FALSE.
#' @return an object of class \code{\linkS4class{RepSeqExperiment}} composed
#' exclusively of private sequences.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' privateclones <- getPrivate(RepSeqData,
#'                             level = "ntClone",
#'                             singletons = FALSE)
#'
getPrivate <- function(x,  level=c("ntCDR3", "aaCDR3", "ntClone", "aaClone"), singletons=FALSE) {
    V1 <- NULL
    chao1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x))
    metaData <- mData(x)
    levelChoice <- match.arg(level)
    keep <- cts[, sum(count>0) == 1, by = levelChoice][V1 == TRUE]
    res <- cts[keep, on = levelChoice][, V1 := NULL]

    if (singletons){
      res<- res[res$count == 1,]
    }

    setkey(res, sample_id)
    nfilter <- nrow(cts) - nrow(res)

    stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
    metaData<- metaData %>%
      dplyr::select(-c(nSequences,ntCDR3,aaCDR3,V,J,VJ,aaClone,ntClone, chao1, iChao))
    sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                        row.names=1,  stringsAsFactors = TRUE)
    sdata <- sdata[match(rownames(metaData), rownames(stats)), ]

    out <- new("RepSeqExperiment",
                    assayData = res,
                    metaData = sdata,
                    otherData = oData(x),
                    History = data.frame(rbind(History(x),
                                        data.frame(history=paste(nfilter, "private", levelChoice, "were extracted with singletons set to",singletons))),
                                        stringsAsFactors=FALSE))
    rm(cts, keep)

    return(out)
}

#' @title Extract public sequences
#'
#' @description This function allows to subset a RepSeqExperiment object in
#' order to keep sequences that are shared by at least two samples:
#'
#' - belonging to a specified group
#'
#' - within the whole dataset
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken
#' into account. Should be one of "aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param group a vector of character indicating the group column name in the
#' mData slot and one experimental group within this column.
#'
#' Samples belonging to the chosen experimental group will be analyzed.
#' The column must be of class factor.
#'
#' Default is NULL, the analysis is performed on the whole dataset.
#'
#'
#' @return an object of class \code{\linkS4class{RepSeqExperiment}} composed
#' exclusively of shared sequences between the specified samples.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' publicSeq <- getPublic(x = RepSeqData,
#'                              level = "aaClone",
#'                              group = c("cell_subset", "amTreg"))
#'
getPublic <- function(x, level=c("ntCDR3", "aaCDR3", "ntClone", "aaClone"),
                        group = NULL) {
  V1 <- NULL
  chao1 <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  cts <- data.table::copy(assay(x))
  metaData <- mData(x)
  if(!is.null(group) & nrow(metaData[metaData[, group[1]] == group[2],])==0) stop("The chosen group is not present in the metaData")
  
  levelChoice <- match.arg(level)
  if (!is.null(group)){
    grp <- metaData[, group[1]]
    grp.name <- group[2]
    sampleNames <- rownames(metaData[grp %in% grp.name,])
    cts <- cts[sample_id %in% sampleNames,]
    keep <- cts[, sum(count > 0) >= 2, by = levelChoice][V1 == TRUE,]
    x.hist <- paste0("getPublic", "; group:", group[1], "; subgroup:", grp.name)
  } else {
    n <- nrow(metaData)
    keep <- cts[, sum(count > 0) >= 2, by = levelChoice][V1 == TRUE,]
    x.hist <- paste0("getPublic", "; group:", group)
  }
  res <- cts[keep, on = levelChoice][, V1 := NULL]
  setkey(res, sample_id)
  rm(cts, keep)

  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,ntCDR3,aaCDR3,V,J,VJ,aaClone,ntClone, chao1, iChao))
  sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                      row.names=1,  stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)), ]

  out <- new("RepSeqExperiment",
             assayData = res,
             metaData = sdata,
             otherData = oData(x),
             History = data.frame(rbind(History(x),
                                  data.frame(history = x.hist)),
                                  stringsAsFactors=FALSE))
  return(out)
}

#' @title Extract the most frequent sequences
#'
#' @description This function extracts the top n sequences within all samples
#' or a group of samples.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken
#' into account. Should be one of aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param prop a numeric between 0 and 1 indicating the proportion of top
#' sequences to extract.
#' @param group character, column name in mData indicating the group on which
#' the extraction will be applied. Must be of class factor. Default is NULL.
#' @return an object of class \code{\linkS4class{RepSeqExperiment}}
#' @export
#' @examples
#'
#' data(RepSeqData)
#' topClones <- getTopSequences(x = RepSeqData,
#'                           level = "aaClone",
#'                           group = c("cell_subset", "amTreg"), prop = 0.1)
#'
#'
getTopSequences <- function(x, level=c("aaClone","ntClone","aaCDR3","ntCDR3"),
                         group = NULL, prop=0.01) {
  V1 <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (prop>1) stop("prop should be a proportion between 0 and 1")
  cts <- data.table::copy(assay(x))
  sdata <- mData(x)
  if(!is.null(group) & nrow(sdata[sdata[, group[1]] == group[2],])==0) stop("The chosen group is not present in the metaData")
  
  levelChoice <- match.arg(level)
  if (!is.null(group)){
    grp <- sdata[, group[1]]
    grp.name <- group[2]
    sampleNames <- rownames(sdata[grp %in% grp.name,])
    cts <- cts[sample_id %in% sampleNames,]
    }

  cts_b <- cts[, .(count = sum(count)), by=c("sample_id", levelChoice)][, ranks := lapply(.SD, frankv, ties.method = "first", order = -1L), by = "sample_id", .SDcols = "count"]
  keep <- cts_b[cts_b[, .I[ranks %in% seq_len(ceiling(.N*prop))], by = "sample_id" ]$V1]

  res <- cts[keep, on = c(levelChoice, "sample_id")][, c("i.count", "ranks") := NULL]

  setkey(res, sample_id)
  nfilter <- nrow(cts) - nrow(res)
  rm(cts, keep)

  sdata <- sdata[sdata$sample_id %in% unique(res$sample_id), ]
  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("aaClone", "V", "J", "VJ", "aaCDR3"), by = "sample_id"], row.names = 1)
  sdata <- setDT(sdata)[, c("nSequences" ,"aaClone" , "V" ,"J" ,"VJ","aaCDR3" ) := NULL]
  sdata<- data.frame(sdata,row.names =sdata$sample_id )
  sdata <- data.frame(base::merge(sdata, stats, by = 0, sort=FALSE), row.names = 1, stringsAsFactors = TRUE)
  # sdata <- sdata[match(rotiwnames(sdata), rownames(stats)),]
  sdata <- sdata[order(match(sdata$sample_id, rownames(stats))),]
  
  out <- new("RepSeqExperiment",
             assayData = res,
             metaData = sdata,
             otherData = oData(x),
             History = data.frame(rbind(History(x),
                                        data.frame(history = paste0("Top", levelChoice,"were extracted with the group parameter set to",group))),
                                  stringsAsFactors=FALSE))
  return(out)
}


#' @title Give type of file
#'
#' @description get the type of file
#'
#' @details This function allows to get the type of a file
#'
#' @param path path to a file
#' @return a string indicating file type
#' @export
#' @keywords internal
#' @examples
#'
#' l <- list.files(system.file(file.path('extdata/MiAIRR'),
#'                      package = 'AnalyzAIRR'),
#'                      full.names = TRUE)
#' path <- l[1]
#'
#' file_type <- filetype(path)
#'
filetype <- function(path) {
    f <- file(path)
    ext <- summary(f)$class
    close.connection(f)
    return(ext)
}



#' @title Merge RepSeqExperiment objects
#'
#' @description This function allows to merge of two RepSeqExperiement objects.
#'
#' The two objects must however contain unique and different sample_ids.
#'
#' This function can be useful in case users wish to analyze alignment files
#' from different formats, or biological/technical replicates.
#' @param a the first  \code{\linkS4class{RepSeqExperiment}} object.
#' @param b the second  \code{\linkS4class{RepSeqExperiment}} object.
#' @return a  \code{\linkS4class{RepSeqExperiment}} object containing all
#' information from the two merged objects.
#' @export
# l <- list.files(system.file(file.path('extdata/MiAIRR'),
#                      package = 'AnalyzAIRR'),
#                      full.names = TRUE)
# 
# metaData <- read.table(system.file(file.path('extdata/sampledata.txt'),
#                          package='AnalyzAIRR'),
#                          sep = "\t",
#                          row.names = 1, header = TRUE)
# 
# dataset1 <- readAIRRSet(fileList = l[c(1,5)],
#                        cores=3,
#                        fileFormat = "MiAIRR",
#                        chain = "TRA",
#                        sampleinfo = metaData[c(1,5),],
#                        filter.singletons = FALSE,
#                        outFiltered = FALSE)
# 
# dataset2 <- readAIRRSet(fileList = l[7:8],
#                        cores=3,
#                        fileFormat = "MiAIRR",
#                        chain = "TRA",
#                        sampleinfo = metaData[7:8,],
#                        filter.singletons = FALSE,
#                        outFiltered = FALSE)
# 
# dataset <- mergeRepSeq(a = dataset1, b = dataset2)
#'
#'
mergeRepSeq <- function(a, b) {
    if (missing(a) | missing (b)) stop("Two RepSeqExperiment objects are required.")
    if (!is.RepSeqExperiment(a)) stop("a is not an object of class RepSeqExperiment.")
    if (!is.RepSeqExperiment(b)) stop("b is not an object of class RepSeqExperiment.")
    if (any(rownames(mData(a)) %in% rownames(mData(b)))) stop("Duplicates in sample names are not allowed, please use the function names()<- to rename them.")
    cts <- rbind(assay(a), assay(b))
    cts[, sample_id := as.character(sample_id)]
    sampleinfo <- rbind(mData(a), mData(b))
    a.history <- paste0(deparse(substitute(a)), ":", History(a)$history)
    b.history <- paste0(deparse(substitute(b)), ":", History(b)$history)
    concat.history <- paste(date(),"- concatenation of", deparse(substitute(a)), "and", deparse(substitute(b)), "using the function concateRepSeq")
    all.history <- data.frame(history=c(a.history, b.history, concat.history))
    metainfo <- ifelse(length(oData(a)) > 0 | length(oData(b)) > 0, c(oData(a), oData(b)), list())
    
    out <- new("RepSeqExperiment",
            assayData = cts,
            metaData = sampleinfo,
            otherData = metainfo,
            History = all.history)
    
    oData(out) <- c(oData(out), label_colors=list(plotColors(x = out)))
    
    return(out)
}

#' @title Exclude repertoires from a RepSeqExperiment object.
#'
#' @description This function allows the dropping of one or several samples
#' from a RepSeqExperiment object by specifying their corresponding sample ids.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleNames a vector of character with the sample_ids of the
#' repertoires to drop from the RepSeqExperiment object.
#' @return a RepSeqExperiment object.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' dropRepSeqData<- dropSamples(x = RepSeqData,
#'                              sampleNames=c("tripod-30-813", "tripod-30-815"))
#'
dropSamples <- function(x, sampleNames) {
    if (missing(x)) stop("A RepSeqExperiment object is required.")
    if (!is.RepSeqExperiment(x)) stop("x is not an object of class RepSeqExperiment.")
    sampleinfo <- mData(x)
    if (is.numeric(sampleNames)) {
        index <- sampleNames
        snames <- rownames(sampleinfo)[index]
    }
    if (is.character(sampleNames)) {
        if (!all(sampleNames %in% rownames(sampleinfo))) stop("Sample names not found in x.")
        index <- which(rownames(sampleinfo) %in% sampleNames)
        snames <- sampleNames
    }
    cts <- copy(assay(x))
    cts <- cts[!(sample_id %in% snames)]
    cts[, sample_id := as.character(sample_id)]
    sampleinfo <- droplevels(sampleinfo[-c(index), , drop=FALSE])
    x.history <- data.frame(rbind(History(x), data.frame(history = paste("sample",paste0(snames, collapse=", "), "was excluded from the RepSeqExperiment."))))
    metainfo <- oData(x)
    out <- new("RepSeqExperiment",
            assayData = cts,
            metaData = sampleinfo,
            otherData = metainfo,
            History = x.history)
    return(out)
}


#' @title Exclude a sequence from a RepSeqExperiment object.
#'
#' @description This function allows to drop of one or several sequences from a
#' RepSeqExperiment object in all samples or a specified group of samples.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken
#' into account. Should be one of "aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param name a vector of character specifying the name(s) of the sequence(s)
#' to filter out from the RepSeqExperiment object.
#' @param group a vector of character indicating the group column name in the
#' mData slot and one experimental group within this column.
#' @return a RepSeqExperiment object.
#' @export
#' @examples
#'
#' RepSeqData<- filterSequence(x = RepSeqData,
#'                             level="aaClone",
#'                             name="TRAV13-2 CAETQSLQRALIGNLQSPISRF TRAJ50",
#'                             group=c("cell_subset","nTreg"))
#'
filterSequence <- function(x, level=c("aaClone","ntClone","aaCDR3","ntCDR3"), 
                           name, 
                           group=NULL) {
  V1 <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  
  levelChoice <- match.arg(level)
  cts <- data.table::copy(assay(x))
  metaData <- mData(x)

  if(! name %in% cts[, get(levelChoice)]) stop("a valid sequence must be specified.")
  
  if (!is.null(group)) {
    if(!is.null(group) & nrow(metaData[metaData[, group[1]] == group[2],])==0) stop("The chosen group is not present in the metaData")
    
    sampleNames <- rownames(metaData[metaData[, group[1]] == group[2],])
    # 
    res <- cts[!(get(levelChoice) == name & sample_id %in% sampleNames)]
    # keep <- cts[, get(levelChoice) == name, by =  c(levelChoice, "sample_id")]
    # keep <- keep[ V1==FALSE | !sample_id %in% sampleNames,]
    # res <- cts[keep, on = c(levelChoice,"sample_id")][, V1:=NULL]
    # setkey(res, sample_id)
    filterout<- cts[get(levelChoice) == name & sample_id %in% sampleNames]
  }  else {
    res <- cts[!get(levelChoice) == name]
    filterout<- cts[get(levelChoice) == name ]
    # keep <- cts[, get(levelChoice) == name, by =  c(levelChoice,"sample_id")][V1 == FALSE, ]
    # res <- cts[keep, on = c(levelChoice,"sample_id")][, V1:=NULL]
    # setkey(res, sample_id)
    # filterout<- cts[, get(levelChoice) == name, by =  c(levelChoice,"sample_id")][V1 == TRUE, ]
  }
  
  nfilter <- nrow(cts) - nrow(res)
  
  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
  
  sdata <- data.frame(base::merge(metaData %>% dplyr::select(-c(nSequences,ntCDR3,aaCDR3,V,J,VJ,aaClone,ntClone)),
                                  stats, by = 0, sort = FALSE),
                               row.names=1,  stringsAsFactors = TRUE)
  sdata <- sdata[order(match(rownames(sdata), rownames(stats))), ]
  
  out <- new("RepSeqExperiment",
             assayData = res,
             metaData = sdata,
             otherData = oData(x),
             History = data.frame(rbind(History(x),
                                        data.frame(history = paste(nfilter,levelChoice,"were filtered using filterSequence: name=",name, "group1=",group[1],"and group2=",group[2]))),
                                        stringsAsFactors=FALSE))
  oData(out) <- c(oData(out), filterSequence=list(filterout))

  
  return(out)
}



#' @title Extract productive sequences
#'
#' @description Extract productive sequences from a RepSeqExperiment object by
#' filtering out unproductive ones. Filtered sequences include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of
#' nucleotides in the ntCDR3 column.
#'
#' - sequences containing stop codons: aaCDR3s with a "*" or "~" symbols.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}.
#' @return a filtered \code{\linkS4class{RepSeqExperiment}} object exclusively
#' containing productive sequences.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' productiveData <- getProductive(x = RepSeqData)
#'
getProductive <- function(x) {
  ntCDR3 <- NULL
  chao1 <- NULL
  metaData <- mData(x)
  if (missing(x)) stop("A RepSeqExperiment object is required.")
  if (!is.RepSeqExperiment(x)) stop("x is not an object of class RepSeqExperiment.")
  cts <- data.table::copy(assay(x))

  indx <- !cts[, nchar(ntCDR3) %% 3 > 0 | grepl("\\*", aaCDR3) | grepl("\\~", aaCDR3)]
  res <- cts[indx, ]

  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols =c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
  pastek <- function(y) list( chao1 = .chao1(y)$chao.est, iChao = iChao(y))
  out <- data.frame(res[, .(count = sum(count)), by = c("sample_id", "ntClone")][, pastek(count), by = "sample_id"], row.names = 1)
  stats <- data.frame(base::merge(stats, out, by = 0, sort=FALSE), row.names = 1, stringsAsFactors = TRUE)
  
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,ntCDR3,aaCDR3,V,J,VJ,aaClone,ntClone, chao1, iChao))
  sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                      row.names=1,  stringsAsFactors = TRUE)
  
  sdata <- sdata[order(match(rownames(sdata), rownames(stats))), ]
  
  out <- new("RepSeqExperiment",
             assayData = res,
             metaData = sdata,
             otherData = oData(x),
             History = data.frame(rbind(History(x),
                                        data.frame(history = "Unproductive sequences were filtered out from the RepSeqExperiment object")),
                                  stringsAsFactors=FALSE))

  return(out)
}


#' @title Extract unproductive sequences
#'
#' @description Extract unproductive sequences from a RepSeqExperiment object,
#' which include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of
#' nucleotides in the ntCDR3 column
#'
#' - sequences containing stop codons: aaCDR3s with a "*" or "~" symbols.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}.
#' @return a filtered \code{\linkS4class{RepSeqExperiment}} object exclusively
#' containing unproductive sequences.
#' @export
#' @examples
#' 
#' data(RepSeqData)
#'
#' unproductiveData <- getUnproductive(x = RepSeqData)
#' 
getUnproductive <- function(x) {
  ntCDR3 <- NULL
  chao1 <- NULL
  if (missing(x)) stop("A RepSeqExperiment object is required.")
  if (!is.RepSeqExperiment(x)) stop("x is not an object of class RepSeqExperiment.")
  cts <- data.table::copy(assay(x))
  metaData <- mData(x)
  indx <- cts[, nchar(ntCDR3) %% 3 > 0 | grepl("\\*", aaCDR3) | grepl("\\~", aaCDR3) ]
  res <- cts[indx, ]

  if(nrow(res) == 0){
    return(res)
    stop("No unproductive sequences were found in the RepSeqExperiment object.")

  } else {
  
  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols =c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,ntCDR3,aaCDR3,V,J,VJ,aaClone,ntClone, chao1, iChao))
  sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                      row.names=1,  stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)), ]

  out <- new("RepSeqExperiment",
             assayData = res,
             metaData = sdata,
             otherData = oData(x),
             History = data.frame(rbind(History(x),
                                        data.frame(history = "Productive sequences were filtered out from the RepSeqExperiment object")),
                                  stringsAsFactors=FALSE))


  return(out)
  }
}


#' @title Color palette
#'
#' @description This function allows an automatic color assigning, chosen from
#' color-blinded friendly palettes, to every sample and experimental group in
#' the metadata.
#' As such, each sample/group will be assigned the same color in all the
#' visualization functions.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @return a list of distinct colors assigned to each sample/group in the mData.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' colors <- plotColors(x = RepSeqData)
#'

plotColors<- function(x){

  col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category  == "qual" ,]
  mycolors = as.vector(unlist(mapply(RColorBrewer::brewer.pal, col_pals$maxcolors, rownames(col_pals))))
  
  mycolors<- mycolors[-which(base::duplicated(mycolors))]
  names=unique(as.vector(mData(x)[, unlist(lapply(mData(x), is.factor)), drop = FALSE]) %>% names())
  
  ann_colors<-vector("list")
  l<- length(unique(mData(x)[,names[1]]))
  
  if(l<=70){
    mycolors_b<- mycolors[seq_len(l)]
    names(mycolors_b) <- unique(mData(x)[,names[1]])
    ann_colors[["sample_id"]]<- mycolors_b
    
  }else{
    new_l<- ceiling(l/length(mycolors))
    colors_b=  rep(mycolors,new_l) 
    mycolors_b<- colors_b[seq_len(l)]
    
    names(mycolors_b) <-unique(mData(x)[,names[1]])
    ann_colors[["sample_id"]]<- mycolors_b
  }
  
  len<- sum(apply(mData(x)[,names[-1]], 2, dplyr::n_distinct))
  if(len>70) stop ("A maximum of 70 colors can be assigned. 
                    The number of different subgroups is higher than 70.")
  
  for (i in unique(names)[-1]) {
    l <- length(unique(mData(x)[, i]))
    mycolors_b <- mycolors[seq_len(l)]
    names(mycolors_b) <- unique(mData(x)[i][[i]])
    ann_colors[[i]] <- mycolors_b
    mycolors <- mycolors[!mycolors %in% mycolors_b]
  }
  
  return(ann_colors)
  
}

# ggName -> changes a string so it is enclosed in back-ticks.
#   This can be used to make column names that have spaces (blanks)
#   or non-letter characters acceptable to ggplot2.
#   This version of the function is vectorized with sapply.
ggname <- function(x) {
    if (!is(x, "character")) {
        return(x)
    }
    y <- vapply(x, function(s) {
    if (!grepl("^`", s)) {
        s <- paste("`", s, sep="", collapse="")
    }
    if (!grepl("`$", s)) {
        s <- paste(s, "`", sep="", collapse="")
    }}, FUN.VALUE = character(1))
    y
}


#' @title A theme for the visualization plots
#'
#' @description A unique plot theme used in all the visualization plots.
#' As this function is generated using the ggplot2 package, users can easily
#' customize all of the generated plots using ggplot2 arguments.
#'
#' @import ggplot2
#' @export
#'
theme_RepSeq<- function(){
    ggplot2::theme_light()+
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size=12),
                   axis.title = ggplot2::element_text(size=11),
                   axis.text = ggplot2::element_text(size=10),
                   legend.key.size = ggplot2::unit(.2, 'cm'),
                   strip.background = ggplot2::element_rect(fill="white", color="gray"),
                   strip.text.x = ggplot2::element_text(size = 10, color = "black"),
                   strip.text.y = ggplot2::element_text(size = 10, color = "black"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "gray89",linetype="dashed",size=0.1))
}


#' @title Calculation of the clonal distribution per interval
#'
#' @description This function calculates the clonal distribution per a set of intervals in all the samples within the dataset. 
#'
#' The column titled "Distribution" calculates the proportion of each interval in the whole repertoire, whereas the one titled "Cumulative frequency" shows the cumulative frequency of the sequences within each interval.
#'
#' This could allow a global view of the repertoire fraction contributing the most to the repertoire. For instance, top sequences belonging to the highest interval often constitute a low fraction of the whole repertoire but contribute more significantly in terms of cumulative frequency in view of their high occurrence.
#'
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the clonal distribution. Should be one of "aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param interval_scale whether intervals should be determined in count or frequency
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' Intervals(x = RepSeqData, level="aaCDR3",  interval_scale="count")
#' 
Intervals <- function(x, 
                      level = c("aaClone","ntClone", "ntCDR3","aaCDR3"),
                      interval_scale=c("count", "frequency")){
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  
  interval=percent <- NULL
  sdata<-mData(x)
  levelChoice <- match.arg(level)
  
  data2plot <- data.table::copy(assay(x))
  data2plot<- data2plot[, lapply(.SD, sum), by = c(levelChoice, "sample_id"), .SDcols = "count"][,frequency := prop.table(count),by="sample_id"]
  
  if(interval_scale=="count"){
 
    f <- function(x){
      if(x == 1) "1"
      else if(x <= 10) "]1, 10]"
      else if(x <= 100) "]10, 100]"
      else if(x <= 1000) "]100, 1000]"
      else if(x <= 10000) "]1000, 10000]"
      else "]10000, Inf]"
    }
  } else if(interval_scale=="frequency"){
   
    f <- function(x){
      if(x <= 0.000001) "]0, 0.000001]"
      else if(x <= 0.00001) "]0.000001, 0.00001]"
      else if(x <= 0.0001) "]0.00001, 0.0001]"
      else if(x <= 0.001) "]0.0001, 0.001]"
      else if(x <= 0.01) "]0.001, 0.01]"
      else "]0.01, 1]"
    }
  }
  
  data2plot <- data2plot[, lapply(.SD, sum), by = c(levelChoice, "sample_id"), .SDcols = interval_scale][, `:=`(interval, unlist(lapply(get(interval_scale), 
                                                                                                                                   f))), by = "sample_id"]
  data2plot_b <- data2plot %>% dplyr::group_by(interval, sample_id) %>% dplyr::summarize(sum=dplyr::n())
  data2plot_b <- data2plot_b %>% dplyr::group_by( sample_id) %>% dplyr::mutate(distribution=sum/sum(sum))
  
  data2plot <- data2plot[, lapply(.SD, sum), by = c("interval",  "sample_id"), .SDcols = interval_scale][, `:=`(cumulative_frequency, prop.table(get(interval_scale))), by = "sample_id"]
  data2plot <- merge(data2plot[,c(1,2,4)], data2plot_b[,c(1,2,4)], by=c("interval", "sample_id"))
  return(data2plot)
}


#' @title Calculation of repertoire dissimilarities
#'
#' @description This function assesses pairwise repertoire dissimilarities using a specific dissimilarity method.
#'
#' It calculates a list of dissimilarity indices, each taking into account different parameters. The proposed methods include:
#' 
#'  The Jaccard similarity: a measure of similarity between sample sets defined as the size of the intersection divided by the size of the union of the sample sets.
#'
#'  The Morisita-Horn similarity: a measure of similarity that tends to be over-sensitive to abundant species.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the indices are computed. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param method a character specifying the distance method to be computed. Should be one of the following: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis."
#' @param binary a boolean indicating whether or not to transform the data into a presence/absence data. Default is FALSE
#'
#' @details Details on the calculated indices can be found in the vegan package: 
#' https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist
#' @export
#' @examples
#'
#'
#' CalcDissimilarity(x = RepSeqData, 
#'                   level = "V", method = "jaccard" )
#'


CalcDissimilarity <- function(x, 
                              level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                              method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
                                         "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
                                         "binomial", "chao", "cao", "mahalanobis"),
                              binary = FALSE) {
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(method)) stop("a distance method is expected.")
  
  variable <- NULL
  levelChoice <- match.arg(level)
  methodChoice <- match.arg(method)
  cols <- c("sample_id", levelChoice, "count")
  tmp <- data.table::copy(assay(x))[, ..cols]
  sdata <- mData(x)
  sNames <- rownames(sdata)

  dat <- data.table::dcast(data = tmp, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)
  simmat <- as.matrix(dat[, vegan::vegdist(t(.SD), method = methodChoice, diag = TRUE, upper = TRUE, binary = binary), .SDcols=sNames])
  return(simmat)
  
  }

# extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(vapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# calculate log scale
log10_points <- function(data, mapping, ..., alpha) {
  ggally_points(data, mapping, ..., alpha=0.5) + scale_x_log10() + scale_y_log10()
}

