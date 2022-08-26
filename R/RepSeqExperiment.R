# without this, there is a warning during check
# data.table variables
utils::globalVariables(c("J", ".", "clone", "sample_id", "V", "J", "VJ","clonotype_id","cell_barcode" ,"..cols", "..n.."))

#------------------------------------------------------------------
# Define the class
#------------------------------------------------------------------
#' @title The class RepSeqExperiment
#'
#' @description An S4 class object enclosing the adaptive immune repertoire data, and which is used in all the analytical metrics proposed by the AnalyzAIRR package. The RepSeqExperiment object is composed of 4 slots, each containing different information.
#'
#' @slot assayData a data.table binding all clonotype tables and containing the following columns:
#'
#'  - sample_id: sample names
#'
#'  - V: Variable gene name
#'
#'  - J: Joining gene name
#'
#'  - CDR3aa: CDR3 amino acid sequence
#'
#'  - CDR3nt: CDR3 nucleotide sequence
#'
#'  - clone: Full clonotype sequence including the V gene, the amino acid CDR3 sequence and the J gene
#'
#'  - VJ: V-J gene combination
#'
#'  - count: the occurrence of the clonotype, i.e the clone sequence
#'
#' @slot metaData a data frame containing sample information specified during the building of the \code{\linkS4class{RepSeqExperiment}} object. Each row represents a sample, and columns the possible information or group that can be attributed to the samples such as the cell population, the donor's age, sex, etc...
#' @slot otherData a list of meta data containing any other information than the user does not wish to include in the analysis.
#' @slot History a data frame registering all operations performed on \code{\linkS4class{RepSeqExperiment}} object, such as the filtering functions and the normalization.
#' @rdname RepSeqExperiment-class
#' @name RepSeqExperiment-class
#' @exportClass RepSeqExperiment
#'
RepSeqExperiment <- setClass("RepSeqExperiment",
    representation = representation(
        assayData = "data.table",
        metaData = "data.frame",
        otherData = "list",
        History = "data.frame"
        )
)

#' Method assay.
#'
#' @param object a RepSeqExperiment object
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod assay
setGeneric("assay", function(object) standardGeneric("assay"))

#' Method assay<-.
#'
#' @param value either numeric, character or data frame
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod assay<-
setGeneric("assay<-", function(object, i, j, value) standardGeneric("assay<-"))

#' Method oData.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod oData
setGeneric("oData", function(object) standardGeneric("oData"))

#' Method oData<-.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod oData<-
setGeneric("oData<-", function(object, value) standardGeneric("oData<-"))

#' Method mData.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod mData
setGeneric("mData", function(object) standardGeneric("mData"))

#' Method mData<-.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod mData<-
setGeneric("mData<-", function(object, value) standardGeneric("mData<-"))

#' Method History.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod History
setGeneric("History", function(object) standardGeneric("History"))

#' Method History<-.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod History<-
setGeneric("History<-", function(object, value) standardGeneric("History<-"))

#------------------------------------------------------------------
# get methods
#------------------------------------------------------------------

# \code{assay} get count data
# @title The method assay is defined in the class RepSeqExperiment
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]
# @return a data.table of assay of clone features (row) across samples (columns).
# @name RepSeqExperiment
# @rdname RepSeqExperiment
# @exportMethod assay

#' @rdname RepSeqExperiment-class
#' @aliases assay
#' @aliases assay,RepSeqExperiment-method
setMethod(f = "assay",
    signature = "RepSeqExperiment",
    definition = function(object) object@assayData
)

# \code{assay<-} set count data
# @title The method assay<- is defined in the class RepSeqExperiment
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]
# @exportMethod assay<-
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases assay<-
#' @aliases assay<-,RepSeqExperiment-method
setReplaceMethod(f = "assay",
    signature = "RepSeqExperiment",
    definition = function(object, i, j, value) {
        set(object@assayData, i, j, value)
        #object@assayData <- value
        object
        }
)

#------------------------------------------------------------------
# other data
#------------------------------------------------------------------
#' @rdname RepSeqExperiment-class
#' @aliases oData
#' @aliases oData,RepSeqExperiment-method
setMethod(f = "oData",
    signature = "RepSeqExperiment",
    definition = function(object) object@otherData
)


#' @rdname RepSeqExperiment-class
#' @aliases oData<-
#' @aliases oData<-,RepSeqExperiment-method
setReplaceMethod(f = "oData",
    signature = "RepSeqExperiment",
    definition = function(object, value) {
        object@otherData <- value
    object
    }
)

#------------------------------------------------------------------
# meta data
#------------------------------------------------------------------

# \code{mData} get metadata
# @title The method mData is defined in the class RepSeqExperiment
# @param object an RepSeqExperiment object [\code{\linkS4class{RepSeqExperiment}}]
# @return a data frame of sample information, samples are in rows and parameters are in columns.
# @exportMethod mData
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases mData
#' @aliases mData,RepSeqExperiment-method
setMethod(f = "mData",
    signature = "RepSeqExperiment",
    definition = function(object) object@metaData
)

# \code{mData<-} update sample data
# @title The method History is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]
# @exportMethod mData<-
# @rdname RepSeqExperiment-class


#' @rdname RepSeqExperiment-class
#' @aliases mData<-
#' @aliases mData<-,RepSeqExperiment-method
setReplaceMethod(
    f = "mData",
    signature = "RepSeqExperiment",
    definition = function(object, value) {
        object@metaData <- value
    object
    }
)
#------------------------------------------------------------------
# History
#------------------------------------------------------------------
# get history of the object
# @title The method History is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]
# @return a data frame of annotation of clonotype features.
# @exportMethod History
# @rdname RepSeqExperiment


#' @rdname RepSeqExperiment-class
#' @aliases History
#' @aliases History,RepSeqExperiment-method
setMethod(f = "History",
    signature = "RepSeqExperiment",
    definition = function(object) object@History
)

# set history
# @title The method History is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]
# @exportMethod History
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases History<-
#' @aliases History<-,RepSeqExperiment-method
setReplaceMethod(f = "History",
    signature = "RepSeqExperiment",
    definition = function(object, value) {
        object@History <- value
        object
        }
)

#------------------------------------------------------------------
# display object
#------------------------------------------------------------------
# display the object
# @title The method show is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]
# @exportMethod show
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases show,RepSeqExperiment-method
setMethod("show", "RepSeqExperiment",
function(object) {
    cts <- assay(object)
    sNames <- unique(cts, by = "sample_id")$sample_id
    V <- sort(unique(cts, by = "V")$V)
    J <- sort(unique(cts, by = "J")$J)
    m <- cts[, .(n = uniqueN(sample_id), V = uniqueN(V), J = uniqueN(J), s = uniqueN(clonotype), VJ = uniqueN(VJ))]
	cat("An object of class \"", class(object), "\"\n", sep="")
	if (m$n < 4) {
	    cat("Sample_ids              :", sNames[seq_len(m$n)], "\n")
	} else {
	    cat("Sample_ids              :", sNames[seq_len(3)], "...", sNames[m$n],"\n")
	}
	cat("Number of sequences :", m$s, "\n")
    cat("Number of V genes           :", m$V,  "\n")
	cat("Number of J genes           :", m$J, "\n")
})

#' @rdname RepSeqExperiment-class
#' @aliases names
#' @aliases names,RepSeqExperiment-method
setMethod(f = "names",
    signature(x="RepSeqExperiment"),
    definition = function(x) {
        rownames(mData(x))
    }
)
#modified by GPI
#' @rdname RepSeqExperiment-class
#' @aliases names<-
#' @aliases names<-,RepSeqExperiment-method
setReplaceMethod(f = "names",
    signature(x="RepSeqExperiment", value="ANY"),
    definition = function(x, value) {
        oldnames <- rownames(mData(x))
        rownames(mData(x)) <- value
        snames <- unique(assay(x)[["sample_id"]])
        for (l in seq_len(length(snames))) {
            set(assay(x), i=which(assay(x)[["sample_id"]] == snames[l]), j="sample_id", value=value[l])
        }
    History(x) <- data.frame(rbind(History(x),
                data.frame(history=paste(date(), "- updated sample names", paste0(snames, collapse=", "), "from", paste0(snames, collapse=", "), "using names()"))))
    x
    }
)
# modified by GPI
#' @rdname RepSeqExperiment-class
#' @name RepSeqExperiment-class
setValidity("RepSeqExperiment", function(object) {
    msg <- NULL
    valid <- TRUE
    if (!identical(unique(assay(object)$sample_id), rownames(mData(object)))) {
        valid <- FALSE
        msg <- c(msg, "column sample_id in assayData must contain metaData row names.")
    }
    if (!any(assay(object)$count %% 1 == 0)) {
        valid <- FALSE
        msg <- c(msg, "some count in assay are not integers.")
    }
    if (valid) TRUE else msg
})


#------------------------------------------------------------------
# Subsetting
#------------------------------------------------------------------

# @title The method [ is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
#' Wrapper functions
#'
#' \code{[]} Extract parts of the \code{\linkS4class{RepSeqExperiment}} object
#' @param i indice(s) of clonotype(s) to extract
#' @param j indice(s) of sample(s) to extract
#' @param drop If \code{TRUE} the result is coerced to the lowest possible dimension
#' @return an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @rdname RepSeqExperiment-class
#' @aliases [,RepSeqExperiment-method
#' @export
#
setMethod(
    f = "[",
    signature(x = "RepSeqExperiment", i = "ANY", j = "ANY"),
    definition = function(x, i, j, drop) {
    	if (missing(j)) {
    		j <- seq_len(nrow(mData(x)))
    	}
    	if (!is.character(j)) s <- rownames(mData(x))[j] else s <- j
        cts <- copy(assay(x))
        cts <- cts[sample_id %in% s, ]
    	out <- new("RepSeqExperiment",
    				assayData=cts,
    				metaData=droplevels(mData(x)[j, , drop=FALSE]),
    				otherData=oData(x),
    				History=data.frame(rbind(History(x),
    				    data.frame(history=paste0("subet by [ ", length(j)," samples were selected from orignal object RepSeqExperiment."), stringsAsFactors=FALSE)))
	   			)
	   out
})

#' \code{is.RepSeqExperiment} check whether an object is [\code{\linkS4class{RepSeqExperiment}}]
#' @param x an object
# @return TRUE if x is an object of class [\code{\linkS4class{RepSeqExperiment}}].
#' @name is.RepSeqExperiment
#' @rdname RepSeqExperiment-class
#' @export
is.RepSeqExperiment <- function(x) inherits(x, "RepSeqExperiment")
