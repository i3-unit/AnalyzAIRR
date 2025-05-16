utils::globalVariables(c("J", ".", "..sNames", ".SD", "variable2", "cdr3", "cdr3_nt", 
                         "v_gene", "j_gene", "barcode","raw_clonotype_id", ".SDcols",
                         "key", ".N", "count", "..keep.cols", "sampleNames", "aaCDR3.length",
                         "aaCDR3", "pct", "ctrl.mean", "ID","prop","clone"))

#' @title parse MiXCR output
#'
#' @description parse clonotype tables from MiXCR
#'
#' @details function imports clonotype tables produced by the MiXCR aligner.
#'
#' @param path full path to the aligned file. Files can be loaded as gzipped.
#' @param chain character, the TCR chain to be analyzed. One of \code{A} or \code{B}. Default value is \code{A}.
#' @return a data.table of 9 columns: \code{sample_id} name, \code{V} V gene, \code{J} J gene,\code{VJ} V-J gene combinations, \code{aaCDR3} amino acid CDR3 sequence, \code{ntCDR3} nucleotide CDR3 sequence, \code{ntClone} full clone sequence in nucleotide, \code{aaClone} full clone sequence in amino acid, \code{count} clone count. Clones are eliminated if a STOP codon (*) is detected in aaCDR3 chain, if the ntCDR3 length is not divisible by 3 or if the ntCDR3 sequence is ambiguous (contains a "N" base).
#' @export
#' @keywords internal

parseMiXCR <- function(path, chain = c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL")) {
  if (path == "" | missing(path))  stop("Empty file name.")
  tab=V=ntCDR3=aaCDR3=J=sample_id <- NULL
  tab <- data.table::fread(path)

  namelist <- colnames(tab)
  aacdr3 <- ifelse(length(grep("aaSeqShortCDR3", namelist)) > 0, "aaSeqShortCDR3", "aaSeqCDR3")
  ntcdr3 <- ifelse(length(grep("nSeqShortCDR3", namelist)) > 0, "nSeqShortCDR3", "nSeqCDR3")
  if("readCount" %in% namelist){
    count <- "readCount"
  } else if("cloneCount" %in% namelist){
    count <- "cloneCount"
  } else {
    stop("Can't find a column with clonotype counts")
  }
    
  if(length(grep("allVHitsWithScore", namelist)) > 0){
    vHits <- "allVHitsWithScore"
  } else if(length(grep("vHitsWithScore", namelist)) > 0){
    vHits <- "vHitsWithScore"
  } else if(length(grep("allvhits", namelist)) > 0){
    vHits <- "allvhits"
  } else if(length(grep("bestvhit", namelist)) > 0){
    vHits <- "bestvhit"
  } else if(length(grep("vhits", namelist)) > 0){
    vHits <- "vhits"
  } else if(length(grep("bestvgene", namelist)) > 0){
    vHits <- "bestvgene"
  } else if(length(grep("vHit", namelist)) > 0){
    vHits <- "vHit"
  } else if(length(grep("vGene", namelist)) > 0){
    vHits <- "vGene"
  } else if(length(grep("vHitScore", namelist)) > 0){
    vHits <- "vHitScore"
  } else if(length(grep("vGenes", namelist)) > 0){
    vHits <- "vGenes"
  } else {
    stop("Can't find a column with V genes")
  }

  if(length(grep("allJHitsWithScore", namelist)) > 0){
    jHits <- "allJHitsWithScore"
  } else if(length(grep("jHitsWithScore", namelist)) > 0){
    jHits <- "jHitsWithScore"
  } else if(length(grep("alljhits", namelist)) > 0){
    jHits <- "alljhits"
  } else if(length(grep("bestjhit", namelist)) > 0){
    jHits <- "bestjhit"
  } else if(length(grep("jhits", namelist)) > 0){
    jHits <- "jhits"
  } else if(length(grep("bestjgene", namelist)) > 0){
    jHits <- "bestjgene"
  } else {
    stop("Can't find a column with J genes")
  }

  keep.cols <- c(ntcdr3, aacdr3, vHits, jHits, count)

  tab <- tab[, ..keep.cols]
  data.table::setnames(tab, c("ntCDR3", "aaCDR3", "V", "J", "count"))
  tab[, V := {t1 = gsub("\\*..", "", V); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
  tab[, J := {t1 = gsub("\\*..", "", J); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
  tab[, aaCDR3 := gsub("\\_", "*", aaCDR3)]

  ch <- match.arg(chain)
  tab <- switch(ch,
                TRA = {
                  tab[grepl("TRA", V) & grepl("TRA", J), ]
                },
                TRB = {
                  tab[grepl("TRB", V) & grepl("TRB", J), ]
                },
                TRG = {
                  tab[grepl("TRG", V) & grepl("TRG", J), ]
                },
                TRD = {
                  tab[grepl("TRD", V) & grepl("TRD", J), ]
                },
                IGH = {
                  tab[grepl("IGH", V) & grepl("IGH", J), ]
                },
                IGK = {
                  tab[grepl("IGK", V) & grepl("IGK", J), ]
                },
                IGL = {
                  tab[grepl("IGL", V) & grepl("IGL", J), ]
                })
  sName <- gsub(".tsv|.txt|.gz|.zip|.tar|.csv", "", basename(path))
  tab$sample_id <- sName
  out <- tab[, c("VJ", "aaClone", "ntClone") := list(paste(V, J), paste(V, aaCDR3, J), paste( V, ntCDR3, J))][, c("sample_id", "V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone" ,"count")]
  return(out)
}

#' @title parse immunoseq output
#'
#' @description parse output tables from immunoseq
#'
#' @details function imports clonotype tables produced by the immunoseq aligner.
#'
#' @param path full path to the aligned file. Files can be loaded as gzipped.
#' @param chain character, the TCR chain to be analyzed. One of \code{A} or \code{B}. Default value is \code{A}.
#' @return a data.table of 9 columns: \code{sample_id} name, \code{V} V gene, \code{J} J gene,\code{VJ} V-J gene combinations, \code{aaCDR3} amino acid CDR3 sequence, \code{ntCDR3} nucleotide CDR3 sequence, \code{ntClone} full clone sequence in nucleotide, \code{aaClone} full clone sequence in amino acid, \code{count} clone count. Clones are eliminated if a STOP codon (*) is detected in aaCDR3 chain, if the ntCDR3 length is not divisible by 3 or if the ntCDR3 sequence is ambiguous (contains a "N" base).
#' @export
#' @keywords internal
#'
parseImmunoseq <- function(path, chain = c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL")) {
  if (path == "" | missing(path))  stop("Empty file name.")
  tab=V=ntCDR3=aaCDR3=J=sample_id=vMaxResolved=jMaxResolved <- NULL
  tab <- data.table::fread(path)
 
  namelist <- colnames(tab)
  aacdr3 <- "aminoAcid"
  ntcdr3 <- "nucleotide"
  keep.cols <- c(ntcdr3, aacdr3, "vMaxResolved", "jMaxResolved", "count (templates/reads)")
  ch <- match.arg(chain)
  tab <- switch(ch,
                TRA = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("TCRA", vMaxResolved)]
                },
                TRB = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("TCRB", vMaxResolved)]
                },
                TRD = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("TCRD", vMaxResolved)]
                },
                TRG = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("TCRG", vMaxResolved)]
                },
                IGH = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("IGH", vMaxResolved)]
                },
                IGK = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("IGK", vMaxResolved)]
                },
                IGL = {
                  tab[, ..keep.cols][!(get(aacdr3) == "")][grep("IGL", vMaxResolved)]
                }
  )
  tab[, vMaxResolved := gsub("\\*..", "", vMaxResolved)][, jMaxResolved := gsub("\\*..", "", jMaxResolved)]
  data.table::setnames(tab, c( "ntCDR3", "aaCDR3", "V", "J", "count"))
  sName <- gsub(".tsv|.txt|.gz|.zip|.tar|.csv", "", basename(path))
  tab$sample_id <- sName
  out <- tab[, c("VJ", "aaClone", "ntClone") := list(paste(V, J), paste(V, aaCDR3, J), paste(V, ntCDR3, J))][, c("sample_id", "ntCDR3", "aaCDR3", "V", "J", "VJ", "aaClone","ntClone" ,"count")]
  return(out)
}

#' @title parse MiAIRR format
#'
#' @description parse MiAIRR format
#'
#' @details function imports the MiAIRR Rearrangement format which is a tab-delimited file format (.tsv) that defines the required and optional annotations for rearranged adaptive immune receptor sequences
#'
#' @param path full path to the aligned file. Files can be loaded as gzipped.
#' @param chain character, the TCR chain to be analyzed. One of \code{A} or \code{B}. Default value is \code{A}.
#' @return a data.table of 9 columns: \code{sample_id} name, \code{V} V gene, \code{J} J gene,\code{VJ} V-J gene combinations, \code{aaCDR3} amino acid CDR3 sequence, \code{ntCDR3} nucleotide CDR3 sequence, \code{ntClone} full clone sequence in nucleotide, \code{aaClone} full clone sequence in amino acid, \code{count} clone count. Clones are eliminated if a STOP codon (*) is detected in aaCDR3 chain, if the ntCDR3 length is not divisible by 3 or if the ntCDR3 sequence is ambiguous (contains a "N" base).
#' @keywords internal
#' @export
#' 
parseAIRRC <- function(path, chain = c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL")) {
  if (path == "" | missing(path))  stop("Empty file name.")
  tab=V=ntCDR3=aaCDR3=J=sample_id <- NULL
  tab <- data.table::fread(path)

  namelist <- colnames(tab)
  n_count <- ifelse(length(grep("duplicate_count", namelist)) > 0, "duplicate_count", "consensus_count")
  keep.cols <- c("junction", "junction_aa", "v_call", "j_call", n_count)

  tab <- tab[, ..keep.cols]
  data.table::setnames(tab, c("ntCDR3", "aaCDR3", "V", "J", "count"))
  tab[, V := {t1 = gsub("\\*..", "", V); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
  tab[, J := {t1 = gsub("\\*..", "", J); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
  tab[, aaCDR3 := gsub("\\_", "*", aaCDR3)]

  ch <- match.arg(chain)
  tab <- switch(ch,
                TRA = {
                  tab[grepl("TRA", V) & grepl("TRA", J), ]
                },
                TRB = {
                  tab[grepl("TRB", V) & grepl("TRB", J), ]
                },
                TRG = {
                  tab[grepl("TRG", V) & grepl("TRG", J), ]
                },
                TRD = {
                  tab[grepl("TRD", V) & grepl("TRD", J), ]
                },
                IGH = {
                  tab[grepl("IGH", V) & grepl("IGH", J), ]
                },
                IGK = {
                  tab[grepl("IGK", V) & grepl("IGK", J), ]
                },
                IGL = {
                  tab[grepl("IGL", V) & grepl("IGL", J), ]
                })
  sName <- gsub(".tsv|.txt|.gz|.zip|.tar|.csv", "", basename(path))
  tab$sample_id <- sName
  out <- tab[, c("VJ", "aaClone", "ntClone") := list(paste(V, J), paste(V, aaCDR3, J), paste(V, ntCDR3, J))][, c("sample_id", "V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone","ntClone" ,"count")]
  return(out)
}


#' @title read in format
#'
#' @description alternative
#'
#' @param path full path to the aligned file. Files can be loaded as gzipped.
#' @param chain character, the TCR chain to be analyzed. One of \code{A} or \code{B}. Default value is \code{A}.
#' @param keep.ambiguous a boolean, whether ambiguous sequences containing an "N" in the CDR3 nucleotide should be left in. Default is \code{FALSE}.
#' @param keep.unproductive a boolean, whether unproductive sequences should be left in. Default is \code{FALSE}.
#' Unproductive sequences include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of nucleotides in the ntCDR3 column
#'
#' - sequences containing stop codons: aaCDR3 sequences with a "*" or "~" symbols.should be kept. Default is \code{FALSE}.
#'
#' @param aa.th an interger, indicates the maximum number of amino acids deviating from the mean length to tolerate. Default is 8. In this case, all amino acid CDR3s sequences with length falling outside the range of mean+- 8 are filtered out.
#' @param outFiltered write in the directory path a data frame containing the filtered out reads. Default path is getwd().
#' @return a data.table of 9 columns: \code{sample_id} name, \code{V} V gene, \code{J} J gene,\code{VJ} V-J gene combinations, \code{aaCDR3} amino acid CDR3 sequence, \code{ntCDR3} nucleotide CDR3 sequence, \code{ntClone} full clone sequence in nucleotide, \code{aaClone} full clone sequence in amino acid, \code{count} clone count. Clones are eliminated if a STOP codon (*) is detected in aaCDR3 chain, if the ntCDR3 length is not divisible by 3 or if the ntCDR3 sequence is ambiguous (contains a "N" base).
#' @export
#' @keywords internal
#'

readInFormat <- function(path,
                         chain = c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL"),
                         keep.ambiguous = FALSE, keep.unproductive = FALSE,
                         aa.th = NULL, outFiltered = FALSE) {
  if (path == "" | missing(path))  stop("Empty file name.")
  tab=V=ntCDR3=aaCDR3=J=sample_id <- NULL
  tab <- data.table::fread(path)

  coltab <- c("sample_id", "V", "J", "ntCDR3", "aaCDR3", "count")

  if (!all(grepl(paste(colnames(tab), collapse = "|"), coltab))) {
    stop("Column names of clonotype table must be sample_id, V, J, ntCDR3, aaCDR3, and count")
  }

  ch <- match.arg(chain)
  tab <- switch(ch,
                TRA = {
                  tab[grepl("TRA", V) & grepl("TRA", J), ]
                },
                TRB = {
                  tab[grepl("TRB", V) & grepl("TRB", J), ]
                },
                TRG = {
                  tab[grepl("TRG", V) & grepl("TRG", J), ]
                },
                TRD = {
                  tab[grepl("TRD", V) & grepl("TRD", J), ]
                },
                IGH = {
                  tab[grepl("IGH", V) & grepl("IGH", J), ]
                },
                IGK = {
                  tab[grepl("IGK", V) & grepl("IGK", J), ]
                },
                IGL = {
                  tab[grepl("IGL", V) & grepl("IGL", J), ]
                })
  tab <- tab[, c("VJ", "aaClone", "ntClone") := list(paste(V, J), paste(V, aaCDR3, J), paste(V, ntCDR3, J))][, c("sample_id",  "V", "J", "VJ", "ntCDR3", "aaCDR3","aaClone","ntClone" ,"count")]

  out <- filterClones(tab, keep.ambiguous = keep.ambiguous,
                          keep.unproductive = keep.unproductive, aa.th = aa.th,
                          outFiltered = outFiltered)

  return(out)
}


#' @title Adaptive immune repertoire filtering
#'
#' @description This function filters adaptive immune repertoires that were not aligned by a AnalyzAIRR-supported tool. The loaded clonotype table must however contain the required columns:
#'
#'  - sample_id: sample names
#'
#'  - V: Variable gene name
#'
#'  - J: Joining gene name
#'
#'  - aaCDR3: amino acid CDR3 sequence
#'
#'  - ntCDR3: nucleotide CDR3 sequence
#'
#'  - aaClone: Full clone sequence including the V gene, the amino acid CDR3 sequence and the J gene
#'
#'  - VJ: V-J gene combination using V and J gene names
#'
#'  - ntClone: Full clone sequence including the V gene, the nucleotide CDR3 sequence and the J gene
#'
#'  - count: the occurrence of the clone
#'
#'  and applies at least one of the following filters:
#'
#'  - filter out ambiguous clones
#'
#'  - filter out unproductive sequences including out-of-frame sequences with stop codons
#'
#'  - filter out singletons, i.e clones with an occurrence of 1
#'
#' @param raw a single clonotype table containing the previously cited required columns
#' @param keep.ambiguous a boolean, whether ambiguous sequences containing an "N" in the CDR3 nucleotide should be left in. Default is \code{FALSE}.
#' @param keep.unproductive a boolean, whether unproductive sequences should be left in. Default is \code{FALSE}.
#' Unproductive sequences include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of nucleotides in the ntCDR3 column
#'
#' - sequences containing stop codons: aaCDR3 sequences with a "*" or "~" symbols.should be kept. Default is \code{FALSE}.
#'
#' @param aa.th an interger, indicates the maximum number of amino acids deviating from the mean length to tolerate. Default is 8. In this case, all amino acid CDR3s sequences with length falling outside the range of mean+- 8 are filtered out.
#' @param outFiltered write in the directory path a data frame containing the filtered out reads. Default path is getwd().
#' @return a filtered clonotype table.
#' @export
#' @keywords internal
#' @examples
#'
#' l <- list.files(system.file(file.path('extdata/MiAIRR'),
#'                      package = 'AnalyzAIRR'),
#'                      full.names = TRUE)
#' path <- l[1]
#'
#' miairr_table <- parseAIRRC(path = path, chain = "TRA")
#'
#' out <- filterClones(miairr_table,
#'                         keep.ambiguous = FALSE,
#'                         keep.unproductive = TRUE,
#'                         aa.th = 8,
#'                         outFiltered = FALSE)
#'
#'
#'
filterClones <- function(raw,
                             keep.ambiguous = FALSE,
                             keep.unproductive = FALSE,
                             aa.th = NULL,
                             outFiltered = FALSE) {
  aaCDR3=ntCDR3 <- NULL
  if (missing(raw)) stop("a clonotype table is required.\n")
  if (!data.table::is.data.table(raw)) stop("a data.table is expected.\n")

  aa.size <- raw[, nchar(aaCDR3)]
  if (!is.null(aa.th)){
  aa.distr <- table(aa.size)
  aa.length <- as.numeric(names(aa.distr))
  aa.mean <- aa.length[which.max(aa.distr)]
  keep.length <- aa.size <= (aa.mean + aa.th) & aa.size >= (aa.mean - aa.th)
  indx <- keep.length
  } else {
    keep.length <- aa.size > 0
    indx <- keep.length
  }
  if (!keep.ambiguous)
    indx <- indx & !raw[, grepl("N", ntCDR3)]
  if (!keep.unproductive)
    indx <- indx & !raw[, nchar(ntCDR3) %% 3 > 0 | grepl("\\*", aaCDR3) | grepl("\\~", aaCDR3)]
  if (outFiltered) {
    #data.table::fwrite(raw[!indx, ], file=file.path(getwd(), paste("filtered_", raw$sample_id[1], ".csv", sep="")), row.names=TRUE)
  out <- list(raw[indx, ], raw[!indx, ])
  } else {
    out <- raw[indx, ]

}
  return(out)
}

#' @title Read a single alignment file
#'
#' @description parse and filter clonotype tables.
#'
#' @details this function is a wrapper of parse functions & filterClones
#'
#' @param path a path to an output file of the supported aligners.
#' @param fileFormat a character, the name of the aligner. Should be one of "immunoseq", "MiXCR" or "MiAIRR". Default is MiXCR.
#' @param chain a character, indicates which TCR chain to import: \code{A} for the alpha chain and \code{B} for the beta chain. Default is \code{A}.
#' @param keep.ambiguous a boolean, whether ambiguous sequences (containing a STOP codon) should be kept. Default is \code{FALSE}.
#' @param keep.unproductive a boolean, whether unproductive sequences (Euclidean dividion of aa length by 3 > 0) should be kept. Default is \code{FALSE}.
#' @param aa.th an interger, indicates the maximum number of amino acids deviating from the mean length to tolerate. Default is NULL.
#' @param outFiltered write in the directory path a data frame containing the columns of the input clonotype table and the filterd out reads as rows. Default path is getwd().
#' @return a filtered clonotype table.
#' @export
#' @keywords internal
#' @examples
#'
#' l <- list.files(system.file(file.path('extdata/MiAIRR'),
#'                      package = 'AnalyzAIRR'),
#'                      full.names = TRUE)
#'
#'
#' dataset <- readAIRR(path = l[1],
#'                           fileFormat = "MiAIRR",
#'                           chain = "TRA",
#'                           keep.ambiguous = FALSE,
#'                           keep.unproductive = FALSE,
#'                          outFiltered = FALSE)
#'
readAIRR <- function(path, fileFormat=c("MiXCR", "immunoseq", "MiAIRR"),
                           chain=c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL"),
                           keep.ambiguous=FALSE, keep.unproductive=FALSE, aa.th=NULL,
                           outFiltered=FALSE) {
  if (!file.exists(path)) stop("Full path to alignment file is required.")
  type <- match.arg(fileFormat)
  ch <- match.arg(chain)
  raw <- switch(type,

                MiXCR = {
                  parseMiXCR(path, chain = ch)
                },
                immunoseq = {
                  parseImmunoseq(path, chain = ch)
                },
                `MiAIRR` = {
                  parseAIRRC(path, chain = ch)
                })
  out <- filterClones(raw,
                          keep.ambiguous = keep.ambiguous,
                          keep.unproductive = keep.unproductive,
                          aa.th = aa.th,
                          outFiltered = outFiltered)

  return(out)
}


#' @title Building of a RepSeqExperiment object
#'
#' @description This function builds a RepSeqExperiment object from alignment files generated by MiXCR and immunoseq, as well as files formatted following the AIRR standards guidelines.
#'
#' A set of filters can be applied on the loaded files during the building process including:
#'
#'  - filtering out ambiguous sequences
#'
#'  - filtering out unproductive sequences including the ones with stop codons, namely out-of-frame sequences
#'
#'  - filtering out singletons, i.e clones with an occurrence of 1
#'
#'  - filtering out short or extensively long amino acid CDR3 sequences
#'
#' @param fileList a list of paths to the alignment files. File format can be one of the following: .tsv, txt.gz, .zip, or .tar
#' @param fileFormat a character vector specifying the format of the input files, i.e. the tool that was used to generate/align the files. Should be one of "MiXCR", "immunoseq", or "MiAIRR".
#' @param chain a character vector indicating a single TCR or Ig chain to analyze. The vector can be one of the following:
#'
#' - "TRA", "TRB","TRG" or "TRD" for the TCR repertoires
#'
#' - "IGH","IGK" or "IGL" for the BCR repertoires
#'
#' @param sampleinfo a data frame containing:
#'
#' - a column with the sample_ids. Ids should match the base names of the corresponding files as well as their order in the list. This column should be assigned as row.names when the metadata file is loaded (see the example below).
#'
#' - any additional columns with relevant information for the analyses. Group columns must be transformed into factors after loading (see the example below).
#'
#' No specific column names nor order are required
#' @param keep.ambiguous a boolean indicating whether or not ambiguous sequences containing an "N" in their nucleotide CDR3 should be left in. Default is \code{FALSE}.
#' @param keep.unproductive a boolean indicating whether or not unproductive sequences should be left in. Default is \code{FALSE}.
#' Unproductive sequences include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of nucleotides in their CDR3s
#'
#' - sequences containing stop codons: aaCDR3 sequences with a "*" or "~" should be kept. Default is \code{FALSE}
#'
#' @param filter.singletons a boolean indicating whether or not ntClones (V+ntCDR3+J) with an occurrence of 1 should be filtered out. Default is \code{FALSE}
#' @param aa.th an integer determining the CDR3 amino acid sequence length limits, i.e. the maximum number of amino acids deviating from the mean length that is accepted. The default value is 8, which keeps amino acid CDR3s with a length that falls inside the following range: mean length-8 =< aa.th <= mean length+ 8.
#' @param outFiltered a boolean indicating whether or not to write in the oData slot of the RepSeqexperiment object a data frame containing the filtered out reads.
#' @param cores an integer indicating the number of CPU cores to be used by the function. The process is paralleled and can thus be used on multiple cores. Default is 1.
#'
#' @return an object of class \code{\link{RepSeqExperiment}} that is  used in all the analytical metrics proposed by the AnalyzAIRR package. See \code{\link{RepSeqExperiment-class}} for more details.
#' @export
#' @examples
#' l <- list.files(system.file(file.path('extdata/MiAIRR'),
#'                      package = 'AnalyzAIRR'),
#'                      full.names = TRUE)
#'
#' metaData <- read.table(system.file(file.path('extdata/sampledata.txt'),
#'                          package='AnalyzAIRR'),
#'                          sep = "\t",
#'                          row.names = 1, header = TRUE)
#'
#' metaData$cell_subset <- factor(metaData$cell_subset)
#' metaData$sex <- factor(metaData$sex)
#'
#' dataset <- readAIRRSet(fileList = l[7:8],
#'                        fileFormat = "MiAIRR",
#'                        chain = "TRA",
#'                        sampleinfo = metaData[7:8,],
#'                        filter.singletons = FALSE,
#'                        aa.th=8,
#'                        outFiltered = FALSE,
#'                        cores=3)
#'
readAIRRSet <- function(fileList, fileFormat = c("MiXCR", "immunoseq",  "MiAIRR"),
                             chain = c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL"),
                             sampleinfo = NULL,
                             keep.ambiguous = FALSE,
                             keep.unproductive = FALSE,
                             filter.singletons = FALSE,
                             aa.th = NULL,
                             outFiltered = FALSE,
                             cores = 1L) {
  sample_id=clone=V=VJ=V1 <- NULL
  cores <- min(parallel::detectCores()-1, cores)
  cat("Running on", cores, "cores...\n")
  if (length(fileList) == 0) stop("Empty list of files, please check folder path.\n")
  if (!is.null(sampleinfo) && length(fileList) != nrow(sampleinfo)) stop("The number of samples in fileList and sampleinfo do not match.\n")
  
  snames <- gsub(".tsv|.txt|.gz|.zip|.tar|.csv", "", basename(fileList))
  parser <- match.arg(fileFormat)
  ch <- match.arg(chain)
  if (cores > 1) {
    Sys.sleep(0.1)
    cat("Loading and filtering sequences")
    cl <- parallel::makeCluster(cores, type = "SOCK", rscript_args = "--vanilla", useXDR = TRUE)
    parallel::clusterExport(cl = cl, varlist = c("filetype", "readAIRR",
                                                 "parseMiXCR", "parseImmunoseq", "parseAIRRC", "filterClones", "setnames"))
    repList <- pbapply::pblapply(cl = cl,
                                 fileList,
                                 readAIRR,
                                 fileFormat = parser,
                                 chain = ch,
                                 keep.ambiguous = keep.ambiguous,
                                 keep.unproductive = keep.unproductive,
                                 aa.th = aa.th,
                                 outFiltered=outFiltered)


    parallel::stopCluster(cl)
  } else {
    cat("Loading and filtering sequences...\n")
    repList <- lapply(fileList, readAIRR, fileFormat = parser,
                      chain = ch,
                      keep.ambiguous = keep.ambiguous,
                      keep.unproductive = keep.unproductive,
                      aa.th = aa.th, outFiltered=outFiltered)
  }
  cat("Assembling sequences...\n")
  if (outFiltered) {
   repdata=lapply(seq_along(repList), function(x) repList[[x]][[1]])
  } else {
    repdata=repList
  }

  countobj <- data.table::rbindlist(repdata)
  indx <- which(unlist(lapply(repdata, nrow)) == 0)
  if (length(indx) > 0) {
    message(length(indx), "clonotype table(s) have no records after filtering: ", paste0(snames[indx], collapse=", "), ".")
    message("These sample(s) will be excluded.")
  }
  sdata <- data.frame(sample_id = snames, row.names = snames, stringsAsFactors = TRUE)
  if (!is.null(sampleinfo)) sdata <- data.frame(base::merge(sdata, sampleinfo, by = 0, sort = FALSE), row.names = 1)
  if (nrow(sdata) != length(fileList)) stop("Number of files to import differ from the number of samples in sampleinfo file.")
  #sdata[vapply(sdata, FUN = is.character, FUN.VALUE = logical(1))] <- lapply(sdata[vapply(sdata, FUN = is.character, FUN.VALUE = logical(1))], as.factor)
  stats <- data.frame(countobj[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c( "V", "J", "VJ","ntCDR3", "aaCDR3","aaClone","ntClone"), by = "sample_id"], row.names = 1)
  pastek <- function(y) list( chao1 = .chao1(y)$chao.est, iChao = iChao(y))
  out <- data.frame(countobj[, .(count = sum(count)), by = c("sample_id", "ntClone")][, pastek(count), by = "sample_id"], row.names = 1)
  stats <- data.frame(base::merge(stats, out, by = 0, sort=FALSE), row.names = 1, stringsAsFactors = TRUE)

  sdata <- data.frame(base::merge(sdata, stats, by = 0, sort=FALSE), row.names = 1, stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)),]
  
  cat("Creating a RepSeqExperiment object...\n")
  x.hist <- data.frame(history = c(paste0("data directory=", dirname(fileList[1])),
                                   paste0("readAIRRSet; cores=", cores,
                                          "; fileFormat=", parser,
                                          "; chain=", ch,
                                          "; ambiguous ", keep.ambiguous,
                                          "; unprod ", keep.unproductive,
                                          "; filter.singletons ", filter.singletons,
                                          "; aa threshold=", aa.th)),
                       stringsAsFactors = FALSE)

  if (outFiltered == TRUE) {
  out <- methods::new("RepSeqExperiment",
                      assayData = countobj,
                      metaData = sdata,
                      otherData = list(filtered=data.table::rbindlist(lapply(seq_along(repList), function(x) repList[[x]][[2]]))),
                      History = x.hist)
  } else {

    out <- methods::new("RepSeqExperiment",
                        assayData = countobj,
                        metaData = sdata,
                        otherData = list(),
                        History = x.hist)
  }

  if (filter.singletons) {
    cat ("Removing singleton sequences...")
    out <- filterCount(out,level="ntClone", n = 1)
  }
  
  oData(out) <- c(oData(out), label_colors=list(plotColors(out)))
  
  cat("Done.\n")
  return(out)
}


#' @title An alternative method for the building of a RepSeqExperiment object
#' @description This function can be used to build a RepSeqExperiment object using aligned files that were not produced by a AnalyzAIRR-supported aligning tool. The loaded clonotype table must however contain the following required columns:
#'
#'  - sample_id: sample names
#'
#'  - V: Variable gene name
#'
#'  - J: Joining gene name
#'
#'  - aaCDR3: amino acid CDR3 sequence
#'
#'  - ntCDR3: nucleotide CDR3 sequence
#'
#'  - aaClone: Full clone sequence including the V gene, the amino acid CDR3 sequence and the J gene
#'  
#'  - VJ: V-J gene combination using V and J gene names
#'
#'  - ntClone: Full clone sequence including the V gene, the nucleotide CDR3 sequence and the J gene
#'
#'  - count: the occurrence of the clone
#'
#' Clonotype tables must only contain a single chain. No paired-chain analysis are provided by the AnalyzAIRR package.
#' Pre-filtered files obtained using the \code{\link{filterClones}} function can be used as input.
#'
#' @param clonotypetab a single clonotype table containing the previously cited columns
#' @param sampleinfo a data frame containing:
#'
#' - a column with the sample names. Names should match the base names of the corresponding files and their order. This column should be assigned as row.names when the metadata file is loaded. See the example below.
#'
#' - any additional columns with relevant information for the subsequent analyses. Group columns must be transformed into factors after loading. See the example below.
#'
#' No specific column names are required.
#' @return an object of class \code{RepSeqExperiment} that is  used in all the analytical metrics proposed by the AnalyzAIRR package. See \code{\link{RepSeqExperiment-class}} for more details.
#' @export
#' @keywords internal
#'
RepSeqExp <- function(clonotypetab, sampleinfo = NULL){
  coltab <- c("sample_id", "V", "J", "ntCDR3", "aaCDR3", "count")
  if (missing(clonotypetab)) stop("clonotype table is missing, a clonotype table is expected.")
  if (!is.data.table(clonotypetab)) setDT(clonotypetab)

  if (!all(grepl(paste(colnames(clonotypetab), collapse = "|"), coltab))) {
    stop("Column names of clonotype table must contain sample_id, V, J, ntCDR3, aaCDR3, VJ, count")
  }

  clonotypetab <- clonotypetab[, c("VJ", "aaClone", "ntClone") := list(paste(V, J), paste(V, aaCDR3, J), paste(V, ntCDR3, J))][, c("sample_id", "V", "J", "VJ","ntCDR3", "aaCDR3",  "aaClone", "ntClone" ,"count")]

  stats <- clonotypetab[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), by = "sample_id", .SDcols = c( "V", "J", "VJ", "ntCDR3", "aaCDR3","aaClone", "ntClone")]
  sNames <- unique(clonotypetab$sample_id)
  if (is.null(sampleinfo)) {
    sampleinfo <- data.frame(stats, row.names = sNames)
  } else {
    sampleinfo <- data.frame(cbind(sampleinfo, stats), row.names = sNames)
  }

  pastek <- function(y) list( chao1 = .chao1(y)$chao.est, iChao = iChao(y))
  out <- clonotypetab[, .(count = sum(count)), by = c("sample_id", "ntClone")][, pastek(count), by = "sample_id"]
  
  sdata <- data.frame(base::merge(sdata, out, by = "sample_id", sort=FALSE), stringsAsFactors = TRUE)
  rownames(sdata)<- sdata$sample_id
  
  
  x.hist <- data.frame(history = paste0("RepSeqExp; clononotypetab=",
                                        deparse(substitute(clonotypetab)), "; sampleinfo=",
                                        deparse(substitute(sampleinfo))),
                       stringsAsFactors = FALSE)

  out <- methods::new("RepSeqExperiment",
                      assayData = clonotypetab,
                      metaData = sampleinfo,
                      otherData = list(),
                      History = x.hist)
  oData(out) <- c(oData(out), label_colors=list(plotColors(out)))
  return(out)
}

#' @title An alternative method for the building of a RepSeqExperiment object
#'
#' @description This function allows the loading of alignment files in other formats than the ones supported by the \code{\link{readAIRRSet}} function.
#' These files should contain the following column names:
#' - sample_id: Sample names
#'
#' - ntCDR3: nucleotide CDR3 sequence
#'
#' - aaCDR3: amino acid CDR3 sequence
#'
#' - V: Variable gene name (IMGT nomenclature)
#'
#' - J: Joining gene name (IMGT nomenclature)
#'
#' - count: the occurrence of each sequence
#'
#' The order of the columns must be respected. Input files can contain additional columns that won't however be taken into account in the generated RepSeqExperiment object.
#'
#' Similarly to the \code{\link{readAIRRSet}} function, \code{\link{readFormatSet}} applies a series of filters during the building process.
#' The possible filters that can be applied include:
#'
#'  - filtering out ambiguous sequences
#'
#'  - filtering out unproductive sequences including the ones with stop codons, namely out-of-frame sequences
#'
#'  - filtering out singletons, i.e clonotypes with an occurrence of 1
#'
#'  - filtering out short or extensively long amino acid CDR3 sequences
#'
#' @param fileList a list of paths to the alignment files. File format can be one of the following: .tsv, txt.gz, .zip, or .tar
#' @param chain a character vector indicating a single TCR or Ig chain to analyze. The vector can be one of the following:
#'
#' - "TRA", "TRB","TRG" or "TRD" for the TCR repertoires
#'
#' - "IGH","IGK" or "IGL" for the BCR repertoires
#'
#' @param sampleinfo a data frame containing:
#'
#' - a column with the sample_ids. Ids should match the base names of the corresponding files as well as their order in the list. This column should be assigned as row.names when the metadata file is loaded (see the example below).
#'
#' - any additional columns with relevant information for the analyses. Group columns must be transformed into factors after loading (see the example below).
#'
#' No specific column names nor order are required
#' @param keep.ambiguous a boolean indicating whether or not ambiguous sequences containing an "N" in their nucleotide CDR3 should be left in. Default is \code{FALSE}.
#' @param keep.unproductive a boolean indicating whether or not unproductive sequences should be left in. Default is \code{FALSE}.
#' Unproductive sequences include:
#'
#' - out-of-frame sequences: sequences with frame shifts based on the number of nucleotides in their CDR3s
#'
#' - sequences containing stop codons: aaCDR3 sequences with a "*" or "~" should be kept. Default is \code{FALSE}
#'
#' @param filter.singletons a boolean indicating whether or not ntClones (V+ntCDR3+J) with an occurrence of 1 should be filtered out. Default is \code{FALSE}
#' @param aa.th an integer determining the CDR3 amino acid sequence length limits, i.e. the maximum number of amino acids deviating from the mean length that is accepted. The default value is 8, which keeps amino acid CDR3s with a length that falls inside the following range: mean length-8 =< aa.th <= mean length+ 8.
#' @param outFiltered a boolean indicating whether or not to write in the oData slot of the RepSeqexperiment object a data frame containing the filtered out reads.
#' @param cores an integer indicating the number of cores to be uses by the function. The process is paralleled and can thus be used on multiple cores. Default is 1.
#'
#' @return an object of class \code{\link{RepSeqExperiment}} that is  used in all the analytical metrics proposed by the AnalyzAIRR package. See \code{\link{RepSeqExperiment-class}} for more details.
#' @export
#' @examples
#' l <- list.files(system.file(file.path('extdata/informat'),
#'                      package = 'AnalyzAIRR'),
#'                      full.names = TRUE)
#'
#'
readFormatSet <- function(fileList,
                         chain = c("TRA", "TRB","TRG","TRD","IGH","IGK","IGL"),
                         sampleinfo = NULL, keep.ambiguous = FALSE,
                         keep.unproductive = FALSE, filter.singletons = FALSE,
                         aa.th = 8,  outFiltered = TRUE, 
                         cores = 1L) {

  sample_id = clone = V = VJ = V1 <- NULL
  if (length(fileList) == 0) stop("Empty list of files, please check folder path.\n")
  if (!is.null(sampleinfo) && length(fileList) != nrow(sampleinfo)) stop("The number of samples in fileList and sampleinfo do not match.\n")
  
  snames <- gsub(".tsv|.txt|.gz|.zip|.tar|.csv", "", basename(fileList))
  ch <- match.arg(chain)
  cores <- min(parallel::detectCores()-1, cores)
  cat("Running on", cores, "cores...\n")

  if (cores > 1) {
    Sys.sleep(0.1)
    cat("Loading and filtering sequences...\n")
    cl <- parallel::makeCluster(cores, type = "SOCK", rscript_args = "--vanilla", useXDR = TRUE)
    parallel::clusterExport(cl = cl, varlist = c("filetype", "readInFormat",
                                                 "filterClones",
                                                 "setnames"))
    repList <- pbapply::pblapply(cl = cl,
                                 fileList,
                                 readInFormat,
                                 chain = ch,
                                 keep.ambiguous = keep.ambiguous,
                                 keep.unproductive = keep.unproductive,
                                 aa.th = aa.th,
                                 outFiltered=outFiltered)

    parallel::stopCluster(cl)
  } else {
    cat("Loading and filtering sequences...\n")
    repList <- lapply(fileList,
                      readInFormat,
                      chain=ch,
                      keep.ambiguous = keep.ambiguous,
                      keep.unproductive = keep.unproductive,
                      aa.th = aa.th,
                      outFiltered=outFiltered)
  }

  if (outFiltered) {
    repdata=lapply(seq_along(repList), function(x) repList[[x]][[1]])
  } else {
    repdata=repList
  }
  cat("Assembling sequences...\n")
  tab <- data.table::rbindlist(repdata)
  indx <- which(unlist(lapply(repdata, nrow)) == 0)
  if (length(indx) > 0) {
    message(length(indx), "clonotype table(s) have no records after filtering: ",
            paste0(snames[indx], collapse = ", "), ".")
    message("These sample(s) will be excluded.")
  }

  sdata <- data.frame(sample_id = snames, row.names = snames,
                      stringsAsFactors = TRUE)
  if (!is.null(sampleinfo))
    sdata <- data.frame(base::merge(sdata, sampleinfo, by = 0,
                                    sort = FALSE), row.names = 1)
  sdata[vapply(sdata, FUN = is.character, FUN.VALUE = logical(1))] <- lapply(sdata[vapply(sdata, FUN = is.character, FUN.VALUE = logical(1))], as.factor)

  if (nrow(sdata) != length(fileList))
    stop("Number of files to import differ from the number of samples in sampleinfo file.")

  stats <- data.frame(tab[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ","ntCDR3", "aaCDR3", "aaClone","ntClone"), by = "sample_id"], row.names = 1)
  sdata <- data.frame(base::merge(sdata, stats, by = 0, sort = FALSE),
                      row.names = 1, stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)), ]

  pastek <- function(y) list( chao1 = .chao1(y)$chao.est, iChao = iChao(y))
  out <- tab[, .(count = sum(count)), by = c("sample_id", "ntClone")][, pastek(count), by = "sample_id"]
  
  sdata <- data.frame(base::merge(sdata, out, by = "sample_id", sort=FALSE), stringsAsFactors = TRUE)
  rownames(sdata)<- sdata$sample_id
  
  cat("Creating a RepSeqExperiment object...\n")

  x.hist <- data.frame(history = c(paste0("data directory=" , dirname(fileList[1])), paste0("readFormatSet;
                                          cores=", cores, "; chain=", ch, "; ambiguous ",
                                          keep.ambiguous, "; unprod ", keep.unproductive, "; filter.singletons ",
                                          filter.singletons, "; aa threshold=", aa.th)), stringsAsFactors = FALSE)

  if (outFiltered) {
    out <- methods::new("RepSeqExperiment",
                        assayData = tab,
                        metaData = sdata,
                        otherData = list(filtered=data.table::rbindlist(lapply(seq_along(repList), function(x) repList[[x]][[2]]))),
                        History = x.hist)
  } else {
    out <- methods::new("RepSeqExperiment",
                        assayData = tab,
                        metaData = sdata,
                        otherData = list(),
                        History = x.hist)
  }


  if (filter.singletons) {
    cat("Removing singleton sequences...")
    out <- filterCount(out, level = "ntClone", n = 1)
  }

  oData(out) <- c(oData(out), label_colors=list(plotColors(out)))
  return(out)
}


#' @title A method for the building of a RepSeqExperiment object using single cell data
#'
#' @description This function
#'
#' @param path a path to the acell ranger filtered contigs file
#'
#' @return an object of class \code{\link{RepSeqExperiment}}. See \code{\link{RepSeqExperiment-class}} for more details.
#' @export
#' @keywords internal
#' 
formatSingleCell<-function(path) 
{
  if (path == "" | missing(path)) 
    stop("Empty file name.")
  tab = V = ntCDR3 = aaCDR3 = J = sample_id = cell_barcode = clonotype_ID <- NULL
  tab <- data.table::fread(path)


  summ<- tab %>% 
  count(barcode) %>%
  filter(n>1)
  tab <- dplyr::inner_join(tab,summ)
  
  namelist <- colnames(tab)
  aacdr3 <- "cdr3"
  ntcdr3 <- "cdr3_nt" 
  vHits <- "v_gene"
  jHits <- "j_gene"
  cellbarcode= "barcode"
  clonotypeid = "raw_clonotype_id"

  keep.cols <- c(ntcdr3, aacdr3, vHits, jHits,cellbarcode,clonotypeid)
  tab <- tab[, ..keep.cols]
  tab <- tab[, `:=`(c("clonotype_id","VJ", "aaClone", "ntClone"), list(paste(raw_clonotype_id),
                                                        paste(v_gene,j_gene),
                                           paste(v_gene, cdr3, j_gene), 
                                           paste(v_gene, cdr3_nt, j_gene)))]
  tab_d<- tab %>%
    dplyr::group_by(clonotype_id) %>%
    dplyr::summarize("ntCDR3"=paste0(unique(cdr3_nt),collapse = '; '),
                     "aaCDR3"=paste0(unique(cdr3),collapse = '; '),
                     "V"=paste0(unique(v_gene),collapse = '; '),
                     "J"=paste0(unique(j_gene),collapse = '; '),
                     "cell_barcode"=paste0(unique(barcode),collapse = '; '),
                     "VJ"=paste0(unique(VJ),collapse = '; '),
                     "aaClone"=paste0(unique(clone),collapse = '; '),
                     "ntClone"=paste0(unique(ntClone),collapse = '; ')) %>%
    dplyr::group_by(clonotype_id) %>%
    dplyr::mutate(count=length(strsplit(cell_barcode, ';')[[1]]))
  
  sName <- stringr::str_split(basename(path),pattern="_filtered_contig_annotations.csv")[[1]][1]
  tab_d$sample_id <- sName  
  out <- tab_d[,   c("sample_id", "V", "J", "VJ","ntCDR3", "aaCDR3",  "aaClone", 
                 "ntClone", "count","clonotype_id","cell_barcode")]
   return(out)
}


#' @title A method for the building of a RepSeqExperiment object using single cell data
#'
#' @description This function
#'
#' @param fileList a path to the cell ranger filtered contigs file
#' @param sampleinfo a path to the cell ranger filtered contigs file
#' @param filter.singletons a path to the cell ranger filtered contigs file
#' @param cores a path to the acell ranger filtered contigs file
#'
#' @return an object of class \code{\link{RepSeqExperiment}}. See \code{\link{RepSeqExperiment-class}} for more details.
#' @export
#' @keywords internal
#' 
readSingleCellSet <- function(fileList,
                              sampleinfo = NULL, 
                              filter.singletons = FALSE,
                              cores = 1L) {
  snames <- unlist(stringr::str_split(basename(fileList),pattern="_filtered_contig_annotations.csv", simplify = TRUE)[,1])
  
  if (length(fileList) == 0) stop("Empty list of files, please check folder path.\n")
  cores <- min(parallel::detectCores()-1, cores)
  cat("Running on", cores, "cores...\n")
  
  if (cores > 1) {
    Sys.sleep(0.1)
    cat("Loading and filtering sequences...\n")
    cl <- parallel::makeCluster(cores, type = "SOCK", rscript_args = "--vanilla", useXDR = TRUE)
    parallel::clusterExport(cl = cl, varlist = c("filetype", "formatSingleCell",
                                                   "readInFormat", "filterClones", "setnames"))
    repList <- pbapply::pblapply(cl = cl,
                                 fileList,
                                 formatSingleCell)
    parallel::stopCluster(cl)
  } else {
    cat("Loading and filtering sequences...\n")
    repList <- lapply(fileList,
                      formatSingleCell )
  }
  repdata=repList
  cat("Assembling sequences...\n")
  tab <- data.table::rbindlist(repdata)
  indx <- which(unlist(lapply(repdata, nrow)) == 0)
  if (length(indx) > 0) {
    message(length(indx), "clonotype table(s) have no records after filtering: ",
            paste0(snames[indx], collapse = ", "), ".")
    message("These sample(s) will be excluded.")
  }
  
  sdata <- data.frame(sample_id = snames, row.names = snames,
                      stringsAsFactors = TRUE)
  if (!is.null(sampleinfo)){
    sdata <- data.frame(base::merge(sdata, sampleinfo, by = 0,
                                    sort = FALSE), row.names = 1)
  sdata[vapply(sdata, FUN = is.character, FUN.VALUE = logical(1))] <- lapply(sdata[vapply(sdata, FUN = is.character, FUN.VALUE = logical(1))], as.factor)
  }
  if (nrow(sdata) != length(fileList))  stop("Number of files to import differ from the number of samples in sampleinfo file.")
  
  stats <- data.frame(tab[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c( "V", "J", "VJ","ntCDR3", "aaCDR3","aaClone","ntClone"), by = "sample_id"], row.names = 1)
  sdata <- data.frame(base::merge(sdata, stats, by = 0, sort = FALSE),
                      row.names = 1, stringsAsFactors = TRUE)
  sdata <- sdata[match(rownames(sdata), rownames(stats)), ]
  
  cat("Creating a RepSeqExperiment object...\n")
  
  x.hist <- data.frame(history = c(paste0("data directory=" , dirname(fileList[1])), paste0("readSingleCellSet;
                                          cores=", cores,  "; filter.singletons ", filter.singletons)), stringsAsFactors = FALSE)

    out <- methods::new("RepSeqExperiment",
                        assayData = tab,
                        metaData = sdata,
                        otherData = list(),
                        History = x.hist)

  if (filter.singletons) {
    cat("Removing singleton sequences...")
    out <- filterCount(out, level = "ntClone", n = 1)
  }
  
    cat("Done.\n")
    return(out)
}


