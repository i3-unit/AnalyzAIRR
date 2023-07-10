utils::globalVariables(c("J", ".", "..sNames", "sdata", "metaData", ".SD", ".SDcols", "key", ".N", "count", "..keep.cols", "sampleNames", "CDR3aa.length", "CDR3aa", "CDR3nt","pct", "ctrl.mean", "ID","prop", "nSequences", "clonotype", "ranks"))


#' @title Computing the occurrence of any repertoire level
#'
#' @description  This function calculates the occurrence of a selected repertoire level and returns the calculated values for all the samples within a RepSeqExperiment object.
#' It takes into account the weight of the studied level, i.e. the number of sequences expressing a certain gene segment, or the count of a sequence in a sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the repertoire level to be analyzed. Should be one of "clone","clonotype","CDR3aa","CDR3nt","V", "J",or "VJ".
#' @param scale a character specifying the type of occurrence to return: "count" or "frequency".
#' @param group a vector of character indicating the group column name in mData and one experimental group within this column. Samples belonging to the chosen experimental group will be analyzed. The column must be of class factor. Default is NULL, values are calculated in all the samples within the dataset.
#'
#' @return a data.table summarizing the count or frequency of the analyzed level. In this table, rows correspond to the repertoire level, and columns correspond to the sample_ids.
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
#'
countFeatures <- function(x,
                          level=c("clone", "V", "J", "VJ", "CDR3aa", "clonotype","CDR3nt"),
                          scale=c("count","frequency"), group=NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    scl <- match.arg(scale)
    metaData <- mData(x)
    if (!is.null(group)) {
      grp <- metaData[, group[1]]
      grp.name <- group[2]
      sampleNames <- rownames(metaData[grp %in% grp.name, ])
      cts <- cts[sampleNames, on = "sample_id"]
    }
    if (scl=="count"){
    out <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~ sample_id")), value.var="count", fun=sum)
    } else {
    cts <- cts[, .(count = sum(count)), by=c("sample_id", levelChoice)][, prop := count/sum(count), by = "sample_id"]
    out <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~sample_id")), value.var = "prop", fill = 0)
    }
    return(out)
}

# @title Gene and clonotype usage
#
# @description compute gene segment or clonotype usage
#
# @details function computes segment or clonotype frequency in each sample. The function takes into account the occurence of the analyzed level.
#
# @param x an object of class \code{\linkS4class{RepSeqExperiment}}
# @param level character, the level of the repertoire to estimate. Should be one of "clone", "CDR3aa","V", "J" or "VJ"
# @return a data.table


# segmentUsage <- function(x, level=c("clone", "V", "J", "VJ","CDR3aa")) {
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
#' @description This function filters out, in each sample, sequences having a count equal to or below a specified threshold.
#' It can be applied on all the samples within a RepSeqExperiment object, or a group of samples belonging to a specific experimental group.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the type of sequences on which the filtering will be applied. Should be one of "clone","clonotype","CDR3aa" or "CDR3nt". For instance, for level="CDR3aa", counts will first be recalculated based on this column. Then, CDR3aa sequences with counts equal or below the n threshold will be excluded.
#' @param n an integer specifying the count threshold below which  sequences will be filtered out. For instance, for n=2, sequences with a count of 1 and 2 will be filtered out.
#' @param group a vector of character indicating the group column name in mData and one experimental group within this column. Samples belonging to the chosen experimental group will be analyzed. The column must be of class factor. Default is NULL, values are calculated in all the samples within the dataset.

#' @return an object of class \code{RepSeqExperiment}
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' filterdata <- filterCount(x=RepSeqData, n=1, level = "clone", group=c("cell_subset", "amTreg"))
#'
#' filterdata <- filterCount(x=RepSeqData, n=1, level = "clonotype", group=NULL)
#'
filterCount <- function(x, level=c("clone","clonotype","CDR3aa","CDR3nt"), n=1, group=NULL) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    metaData <- mData(x)

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

    stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("CDR3nt", "CDR3aa", "V", "J", "VJ","clone","clonotype"), by = "sample_id"], row.names = 1)
    metaData<- metaData %>%
      dplyr::select(-c(nSequences,CDR3nt,CDR3aa,V,J,VJ,clone,clonotype))
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
#' @description This function extract private sequences within the whole dataset, i.e. sequences found exclusively in one sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account. Should be one of clone","clonotype", "CDR3nt" or "CDR3aa".
#' @param singletons a boolean indicating whether or not private sequences with a count of 1 should be extracted. Default is FALSE.
#' @return an object of class \code{\linkS4class{RepSeqExperiment}} composed exclusively of private sequences.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' privateclonotypes <- getPrivate(RepSeqData, level = "clonotype", singletons = FALSE)
#'
getPrivate <- function(x,  level=c("clone","clonotype","CDR3aa","CDR3nt"), singletons=FALSE) {
    V1 <- NULL
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

    stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("CDR3nt", "CDR3aa", "V", "J", "VJ","clone","clonotype"), by = "sample_id"], row.names = 1)
    metaData<- metaData %>%
      dplyr::select(-c(nSequences,CDR3nt,CDR3aa,V,J,VJ,clone,clonotype))
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
#' @description This function allows to subset a RepSeqExperiment object in order to keep sequences that are shared by at least two samples:
#'
#' - belonging to a specified group
#'
#' - within the whole dataset
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account. Should be one of "clone","clonotype", "CDR3nt" or "CDR3aa".
#' @param group a vector of character indicating the group column name in the mData slot and one experimental group within this column.
#'
#' Samples belonging to the chosen experimental group will be analyzed. The column must be of class factor.
#'
#' Default is NULL, the analysis is performed on the whole dataset.
#'
#'
#' @return an object of class \code{\linkS4class{RepSeqExperiment}} composed exclusively of shared sequences between the specified samples.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' publicSeq <- getPublic(x = RepSeqData,
#'                              level = "clone",
#'                              group = c("cell_subset", "amTreg"))
#'
getPublic <- function(x, level=c("clone","clonotype","CDR3aa","CDR3nt"),
                        group = NULL) {
  V1 <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  cts <- data.table::copy(assay(x))
  metaData <- mData(x)
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

  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("CDR3nt", "CDR3aa", "V", "J", "VJ","clone","clonotype"), by = "sample_id"], row.names = 1)
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,CDR3nt,CDR3aa,V,J,VJ,clone,clonotype))
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
#' @description This function extracts the top n sequences within all samples or a group of samples.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account. Should be one of clone","clonotype", "CDR3nt" or "CDR3aa".
#' @param prop a numeric between 0 and 1 indicating the proportion of top sequences to extract.
#' @param group character, column name in sData. Must be of class factor. Default is NULL, the threshold is calculated across all the samples within the dataset.
#' @return an object of class \code{\linkS4class{RepSeqExperiment}}
#' @export
#' @examples
#'
#' data(RepSeqData)
#' topClones <- getTopSequences(x = RepSeqData,
#'                           level = "clone",
#'                           group = c("cell_subset", "amTreg"), prop = 0.1)
#'
#'
getTopSequences <- function(x, level=c("clone","clonotype","CDR3aa","CDR3nt"),
                         group = NULL, prop=0.01) {
  V1 <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  cts <- data.table::copy(assay(x))
  sdata <- mData(x)
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
  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("clone", "V", "J", "VJ", "CDR3aa"), by = "sample_id"], row.names = 1)
  sdata <- setDT(sdata)[, c("nSequences" ,"clone" , "V" ,"J" ,"VJ","CDR3aa" ) := NULL]
  sdata<- data.frame(sdata,row.names =sdata$sample_id )
  sdata <- data.frame(base::merge(sdata, stats, by = 0, sort=FALSE), row.names = 1, stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)),]

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
#' l <- list.files(system.file(file.path('extdata/mixcr'),
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
#' @description This function allows the merging of two RepSeqExperiement objects.
#'
#' The two objects must however contain completely different sample_ids, as duplicates are not supported.
#'
#' This function can be useful in case users wish to analyze alignment files from different formats, or biological/technical replicates.
#' @param a the first  \code{\linkS4class{RepSeqExperiment}} object.
#' @param b the second  \code{\linkS4class{RepSeqExperiment}} object.
#' @return a  \code{\linkS4class{RepSeqExperiment}} object containing all information from the two merged objects.
#' @export
#' @examples
#' l <- list.files(system.file(file.path('extdata/mixcr'),
#'                      package = 'AnalyzAIRR'),
#'                      full.names = TRUE)
#'
#' metaData <- read.table(system.file(file.path('extdata/sampledata.txt'),
#'                          package='AnalyzAIRR'),
#'                          sep = "\t",
#'                          row.names = 1, header = TRUE)
#'                          
#' dataset1 <- readAIRRSet(fileList = l[c(1:3)],
#'                        cores=1L,
#'                        fileFormat = "MiXCR",
#'                        chain = "TRA",
#'                        sampleinfo = metaData[1:3,],
#'                        filter.singletons = FALSE,
#'                        outFiltered = FALSE,
#'                        raretab = FALSE)
#'
#' dataset2 <- readAIRRSet(fileList = l[c(4:8)],
#'                        cores=1L,
#'                        fileFormat = "MiXCR",
#'                        chain = "TRA",
#'                        sampleinfo = metaData[4:8,],
#'                        filter.singletons = FALSE,
#'                        outFiltered = FALSE,
#'                        raretab = FALSE)
#'
#' dataset <- mergeRepSeq(a = dataset1, b = dataset2)
#'
#'
mergeRepSeq <- function(a, b) {
    if (missing(a) | missing (b)) stop("Two RepSeqExperiment objects are required.")
    if (!is.RepSeqExperiment(a)) stop("a is not an object of class RepSeqExperiment.")
    if (!is.RepSeqExperiment(b)) stop("b is not an object of class RepSeqExperiment.")
    if (any(rownames(mData(a)) == rownames(mData(b)))) stop("Duplicates in sample names are not allowed, please use the function names()<- to rename them.")
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
    
    return(out)
}

#' @title Exclude repertoires from a RepSeqExperiment object.
#'
#' @description This function allows the dropping of one or several samples from a RepSeqExperiment object by specifying their corresponding sample ids.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleNames a vector of character with the sample_ids of the repertoires to drop from the RepSeqExperiment object.
#' @return a RepSeqExperiment object.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' dropRepSeqData<- dropSamples(x = RepSeqData,
#'             sampleNames=c("tripod-30-813", "tripod-30-815"))
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
#' @description This function allows the dropping of one or several sequences from a RepSeqExperiment object in all samples or a specified group of samples.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account. Should be one of "clone","clonotype", "CDR3nt" or "CDR3aa".
#' @param name a vector of character specifying the name(s) of the sequence(s) to filter out from the RepSeqExperiment object.
#' @param group a vector of character indicating the group column name in the mData slot and one experimental group within this column.
#' @return a RepSeqExperiment object.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' RepSeqData<- filterSequence(x = RepSeqData,
#'                             level="clone",
#'                             name="TRAV11 CVVGDRGSALGRLHF TRAJ18",
#'                             group=c("cell_subset","Teff"))
#'
filterSequence <- function(x, level=c("clone","clonotype","CDR3aa","CDR3nt"), name, group=NULL) {
  V1 <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  
  levelChoice <- match.arg(level)
  cts <- data.table::copy(assay(x))
  metaData <- mData(x)
  
  if(! name %in% cts[, get(levelChoice)]) stop("a valid sequence must be specified.")
  
  if (!is.null(group)) {
    grp <- metaData[, group[1]]
    grp.name <- group[2]
    sampleNames <- rownames(metaData[grp %in% grp.name, ])
    keep <- cts[, get(levelChoice) == name, by =  c(levelChoice,"sample_id")]
    keep <- keep[ V1==FALSE | !sample_id %in% sampleNames,]
    res <- cts[keep, on = c(levelChoice,"sample_id")][, V1:=NULL]
    setkey(res, sample_id)
    filterout<- cts[, get(levelChoice) == name, by =  c(levelChoice,"sample_id")][V1 == TRUE & sample_id %in% sampleNames, ]
  }  else {
    keep <- cts[, get(levelChoice) == name, by =  c(levelChoice,"sample_id")][V1 == FALSE, ]
    res <- cts[keep, on = c(levelChoice,"sample_id")][, V1:=NULL]
    setkey(res, sample_id)
    filterout<- cts[, get(levelChoice) == name, by =  c(levelChoice,"sample_id")][V1 == TRUE, ]
  }
  
  nfilter <- nrow(cts) - nrow(res)
  
  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("CDR3nt", "CDR3aa", "V", "J", "VJ","clone","clonotype"), by = "sample_id"], row.names = 1)
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,CDR3nt,CDR3aa,V,J,VJ,clone,clonotype))
  sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                      row.names=1,  stringsAsFactors = TRUE)
  sdata <- sdata[order(match(rownames(sdata), rownames(stats))), ]
  
  out <- new("RepSeqExperiment",
             assayData = res,
             metaData = sdata,
             otherData = oData(x),
             History = data.frame(rbind(History(x),
                                        data.frame(history = paste(nfilter,levelChoice,"were filtered using filterSequence: name=",name, "group1=",group[1],"and group2=",group[2]))),
                                        stringsAsFactors=FALSE))
  oData(out) <- c(oData(out), filterSequence=list(  cts[filterout, on = c(levelChoice,"sample_id")][, V1:=NULL]))
  
  rm(cts, keep)
  
  return(out)
}



#' @title Extract productive sequences
#'
#' @description Extract productive sequences from a RepSeqExperiment object by filtering out unproductive ones. Filtered sequences include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of nucleotides in the CDR3nt column.
#'
#' - sequences containing stop codons: CDR3aa sequences with a "*" or "~" symbols.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}.
#' @return a filtered \code{\linkS4class{RepSeqExperiment}} object exclusively containing productive sequences.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' productiveData <- getProductive(x = RepSeqData)
#'
getProductive <- function(x) {
  CDR3nt <- NULL
  metaData <- mData(x)
  if (missing(x)) stop("A RepSeqExperiment object is required.")
  if (!is.RepSeqExperiment(x)) stop("x is not an object of class RepSeqExperiment.")
  cts <- data.table::copy(assay(x))

  indx <- !cts[, nchar(CDR3nt) %% 3 > 0 | grepl("\\*", CDR3aa) | grepl("\\~", CDR3aa)]
  res <- cts[indx, ]

  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("CDR3nt", "CDR3aa", "V", "J", "VJ","clone","clonotype"), by = "sample_id"], row.names = 1)
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,CDR3nt,CDR3aa,V,J,VJ,clone,clonotype))
  sdata <- data.frame(base::merge(metaData, stats, by = 0, sort = FALSE),
                      row.names=1,  stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)), ]

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
#' @description Extract unproductive sequences from a RepSeqExperiment object, which include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of nucleotides in the CDR3nt column
#'
#' - sequences containing stop codons: CDR3aa sequences with a "*" or "~" symbols.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}.
#' @return a filtered \code{\linkS4class{RepSeqExperiment}} object exclusively containing unproductive sequences.
#' @export
#' @examples
#' 
#' \dontrun{
#' data(RepSeqData)
#'
#' unproductiveData <- getUnproductive(x = RepSeqData)
#' 
#' }
#' 
getUnproductive <- function(x) {
  CDR3nt <- NULL
  if (missing(x)) stop("A RepSeqExperiment object is required.")
  if (!is.RepSeqExperiment(x)) stop("x is not an object of class RepSeqExperiment.")
  cts <- data.table::copy(assay(x))
  metaData <- mData(x)
  indx <- cts[, nchar(CDR3nt) %% 3 > 0 | grepl("\\*", CDR3aa) | grepl("\\~", CDR3aa) ]
  res <- cts[indx, ]

  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("CDR3nt", "CDR3aa", "V", "J", "VJ","clone","clonotype"), by = "sample_id"], row.names = 1)
  metaData<- metaData %>%
    dplyr::select(-c(nSequences,CDR3nt,CDR3aa,V,J,VJ,clone,clonotype))
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


#' @title Color palette
#'
#' @description This function allows an automatic color assigning, chosen from color-blinded friendly palettes, to every sample and experimental group in the metadata.
#' As such, each sample/group will be assigned the same color in all the visualization functions.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @return a list of distinct colors assigned to each sample/group in the metadata.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' colors <- plotColors(x = RepSeqData)
#'
plotColors<- function(x, samplenames=TRUE){
  # qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$maxcolors  == 8 & RColorBrewer::brewer.pal.info$colorblind == 'TRUE',]
  # PAIRED = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[c(2,4,6,8,10,12,1,3,5,7,9,11)]
  # mycolors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(as.vector(as.matrix(mData(x)[, unlist(lapply(mData(x), is.factor)), drop = FALSE])) %>% unique() %>% length())
  # mycolors = as.vector(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

  col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category  == "qual" ,]
  mycolors = as.vector(unlist(mapply(RColorBrewer::brewer.pal, col_pals$maxcolors, rownames(col_pals))))
  
  names=as.vector(mData(x)[, unlist(lapply(mData(x), is.factor)), drop = FALSE]) %>% names()
  
  if(samplenames){
    ann_colors<-vector("list")
    l<- length(unique(mData(x)[,names[1]]))
  
    if(l<74){
    mycolors_b<- mycolors[seq_len(l)]
    names(mycolors_b) <-levels(mData(x)[[1]])
    ann_colors[["sample_id"]]<- mycolors_b
    }else{
    new_l<- ceiling(l/length(mycolors))
    colors_b=  rep(mycolors,new_l) 
    mycolors_b<- colors_b[seq_len(l)]
    
    names(mycolors_b) <-levels(mData(x)[[1]])
    ann_colors[["sample_id"]]<- mycolors_b
      
    }
  
  } else {
    len<- sum(apply(mData(x)[,names[-1]], 2, dplyr::n_distinct))
    if(len>74) stop ("A maximum of 74 colors can be assigned. The number of different subgroups is higher than 74.")
      
    ann_colors<-vector("list")
    for (i in unique(names)[-1]) {
      l <- length(unique(mData(x)[, i]))
      mycolors_b <- mycolors[seq_len(l)]
      names(mycolors_b) <- levels(mData(x)[i][[i]])
      ann_colors[[i]] <- mycolors_b
      mycolors <- mycolors[!mycolors %in% mycolors_b]
    }
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
#' As this function is generated using the ggplot2 package, users can easily customize all of the generated plots using ggplot2 arguments.
#'
#' @import ggplot2
#' @export
#'
theme_RepSeq<- function(){
    ggplot2::theme_light()+
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size=12),
                   axis.title = ggplot2::element_text(size=11),
                   axis.text = ggplot2::element_text(size=10),
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(.2, 'cm'),
                   strip.background = ggplot2::element_rect(fill="white", color="gray"),
                   strip.text.x = ggplot2::element_text(size = 10, color = "black"),
                   strip.text.y = ggplot2::element_text(size = 10, color = "black"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "gray89",linetype="dashed",size=0.1))
}

# extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
