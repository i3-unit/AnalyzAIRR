utils::globalVariables(c("ranks", "exp_shannon", "chao1"))

#' @title Calculation of diversity indices
#'
#' @description This function computes a set of diversity indices of a
#' repertoire level for each sample.
#'
#' The calculated indices are the following:
#'
#' - Shannon index: Calculates the proportional abundance of species in a
#' repertoire (Shannon, 1948).
#'
#' - Simpson index: Takes into account the number of species present as well as
#' their abundance. It gives relatively little weight to the rare species and
#' more weight to the frequent ones (Simpson, 1949).
#'
#' - Inverse Simpson index: Is the effective number of species that is obtained
#' when the weighted arithmetic mean is used to quantify average proportional
#' abundance of species.
#'
#' - Berger-Parker index: Expresses the proportional importance of the most
#' abundant species. This metric is highly biased by sample size and richness
#' (Berger and Parker 1970).
#'
#' - Gini coefficient: Measures the degree of inequality in a distribution of
#' abundances (Gini, 1921).
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the
#' diversity should be estimated. Should be one of "aaClone","ntClone", "V",
#' "J", "VJ", "ntCDR3" or "aaCDR3".
#' @return a data table with the diversity indices calculated for each sample.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' diversityIndices(RepSeqData, level="V")
#'
diversityIndices <- function(x,  level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3")) {
    if (missing(x)) stop("x is missing. An object of class RepSeqExperiment is expected.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    pastek <- function(y) list(shannon = .diversity(y, index = "shannon"),
                               simpson = .diversity(y, index = "simpson"),
                               invsimpson = .diversity(y, index = "invsimpson"),
                               bergerparker= .renyiCal(y, alpha=Inf, hill = FALSE),
                               gini = .gini(y)
                               #chao1 = .chao1(y)$chao.est,
                               #chao1.se = .chao1(y)$chao.se,
                               #iChao = iChao(y)
                               )
    out <- copy(assay(x))[, .(count = sum(count)), by = c("sample_id", levelChoice)][, pastek(count), by = "sample_id"]
    return(out)
}

.diversity <- function(x, index=c("shannon", "simpson", "invsimpson"), norm=FALSE, base=exp(1)) {
    if (missing(x)) stop("data set x is required.")
    x <- x/sum(x)
    id <- match.arg(index)
    if (id == "shannon") {
        x <- -x * log(x, base)
        } else {
            x <- x * x
        }
    H <- sum(x, na.rm = TRUE)
    if (norm) H <- H/log(sum(x>0), base)
    if (id == "simpson") {
        H <- 1 - H
    } else if(id == "invsimpson"){
        H <- 1 / H
    } 

    H2 <- round(H, digits = 3)
    return(H2)
}

#' @title Diversity index
#'
#' @description calculate a sublist of diversity indices on gene segments.
#'
#' @details function computes the diversity indices for level \code{VJ},
#' \code{V} or \code{J} only.
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param index character, the diversity index to calculate. Should be one of
#' "chaowor","shannon","invsimpson","simpson" or "gini".
#' @param level character, the level of the repertoire to estimate. 
#' Should be one of "V", "J" or "VJ".
#' @param norm boolean, whether to compute the normalized diversity index,
#' which divides the computed value by the log(S), S being the number of
#' observed species in the repetoire. Default is TRUE.
#' @param base a positive or complex number, the base with respect to which
#' logarithms are computed. Default is e=exp(1).
#' @return a data.table
#' @export
#' @keywords internal
#' @examples
#'
#' data(RepSeqData)
#'
#' divLevel(x = RepSeqData, index = "shannon", level = "V", norm = TRUE)
#'
#' divLevel(x = RepSeqData, index = "simpson", level = "J", norm = FALSE)
#'
divLevel <- function(x, index = c("shannon", "simpson", "invsimpson"), level = c("VJ", "V", "J"), norm = TRUE, base = exp(1)) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected.")
    divIndex <- match.arg(index)
    levelChoice <- match.arg(level)
    res <- assay(x)[, .diversity(count, index=index, base=base, norm=norm), by=c("sample_id", levelChoice)]
    out <- data.table::dcast(res, as.formula(paste0(levelChoice, "~sample_id")), value.var="V1", fill=0)
    return(out)
}

# compute Renyi or hill number
#
# function computes Renyi (Hill's numbers) indices according to \alpha.
# @param x a vector of counts.
# @param alpha a number between 0 and infinity, alpha is Renyi's parameter.
# @param hill a boolean if TRUE the Hill's index is computed.
# @return a number
# @export
.renyiCal <- function(x, alpha=2, hill = FALSE) {
    if (missing(x)) stop("x is missing.")
    if (!is.numeric(x)) stop("x musts be a numerical vector.")
    x <- x/sum(x)
    if (alpha != 0 && alpha != 1 && alpha != Inf) {
        res <- log(sum(x^alpha))/(1 - alpha)
        } else {
            if (alpha == 0) {
                res <- log(sum(x > 0))
            } else if (alpha == Inf) {
                    res <- -log(max(x))
                    } else {
                        res <- sum(-x * log(x, exp(1)), na.rm=TRUE)
                        }
            }
    if (hill) {
      res <- exp(res)
      return(res)
      } else {
        res2<- round(res, digits = 2)
        return(res2)
      }

}

#' @title Computation of the Renyi or the Hill index
#'
#' @description This function computes the Renyi values at any repertoire level
#' for all the samples in the RepSeqExperiment object.
#'
#' The alpha values for which the Renyi is to be estimated can be personalized,
#' thus allowing to focus on certain indices such as the Shannon index for
#' alpha=1 or the Simpson index for alpha=2.
#'
#' @param x an object of class RepSeqExperiment
#' @param alpha a numerical vector specifying the alpha values to compute. If
#' not specified, the following values are estimated: c(0, 0.25, 0.5, 1, 2, 4,
#' 8, 16, 32, 64, Inf).
#' @param level a character specifying the level of the repertoire to be taken
#' into account when calculate VJ usages. Should be one of "aaClone","ntClone",
#' "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param Hill, boolean indicating whether to plot the Renyi or the Hill Index
#' @return a table with Renyi values calculated for all alpha in each sample.
#' @details Hill diversity uses a parameter q to adjust how species are count. 
#' If q=0, it just counts the number of species (species richness). If q=1, 
#' it gives equal weight to all species, balancing rare and common ones. If q=2, 
#' it focuses more on the most common species.
#' Rényi diversity is a generalization of the Shannon index. It uses a parameter 
#' alpha that works similarly. When alpha is small, it gives more weight to rare 
#' species. When alpha is large, it focuses more on the dominant species. Alpha =1 
#' is an approximation of the Shannon index; alpha = 2 corresponds to the Simpson
#' index and alpha=Inf corresponds to the Berger-Parker index.
#' The main difference is that Hill diversity is designed to give “effective
#' species numbers” i.e. how many equally abundant species would give the same
#' diversity, while Rényi diversity is based on entropy and doesn't directly
#' translate to species counts. Both are useful for exploring how diversity
#' changes depending on whether you care more about rare or common species.
#'
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' generalizedDiversity(RepSeqData, level = "V", alpha = 1)
#'
generalizedDiversity <- function(x, alpha = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), 
                       level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                       Hill=FALSE) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is exprected.")
    levelChoice <- match.arg(level)
    out <- copy(assay(x))[, .(count=sum(count)), by=c("sample_id", levelChoice)][, lapply(alpha, function(y) .renyiCal(count, y, hill = Hill)), by="sample_id"]
    
    data.table::setnames(out, c("sample_id", alpha))
    out<-  data.table::melt(out, id.vars = "sample_id")
    out <- data.table::dcast(out, variable ~ sample_id)
    return(out)
}

# compute Gini's coefficient
#
#
#
.gini <- function(x) {
    x <- sort(x)
    n <- length(x)
    out <- 1/n * (n + 1 - 2 * sum((n + 1 - seq_len(n)) * x)/sum(x))
    out2 <- round(out, digits = 2)
    out2
}

# compute Chao1 indices
#
# function computes Chao indices
#
# @param x a vector of counts
# @return a list containing estimated Chao1's diversity and its standard error.
.chao1 <- function(x) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    S <- sum(x > 0)
    #if (f2 > 0) est <- S + f1^2/(2*f2)
    #if (f2 == 0) est <- S + f1*(f1 - 1)/2
    est <- S + f1*(f1 - 1)/(2 * (f2 + 1))
    r <- f1/f2
    chao.var <- f2 * ((r^4)/4 + r^3 + (r^2)/2)
    chao.se <- sqrt(chao.var)
    
    return(list(chao.est=round(est, 2), chao.se=chao.se))
}

#' @title Improved Chao1
#'
#' @description calculate the improved Chao1 index
#'
#' @details function computes the improved version of Chao1 (Chao and Lin, 2012)
#'
#' @param x a vector of count. Can be a list clone counts within a given sample.
#' @return the computed iChao1 value
#' @export
#' @keywords internal
#' @examples
#' x <- rbinom(20, 5, 0.5)
#' iChao(x)
iChao <- function(x) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    f3 <- sum(x == 3)
    f4 <- sum(x == 4)
    if (f4 == 0) f4 <- 1
    n <- sum(x)
    p1 <- (n-3)/n
    p2 <- (n-3)/(n-1)
    est <- .chao1(x)$chao.est + p1*f3/(4*f4) * max(f1 - p2*f2*f3/(2*f4), 0)
    return(round(est, 2))
}

#' @title Adjusted Chao1
#'
#' @description calculate the adjusted Chao1 index
#'
#' @details function computes the adjusted Chao1 index for sampling without
#' replacement
#'
#' @param x a vector of count. Can be a list clone counts within a given sample.
#' @return the computed Chao1 value
#' @export
#' @keywords internal
#' @examples
#' x <- rbinom(20, 5, 0.5)
#' chaoWor(x)
chaoWor <- function(x) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    S <- sum(x > 0)
    n <- sum(x)
    q <- n/S
    dn1 <- n*2*f2/(n-1)
    dn2 <- q*f1/(q-1)
    est <- S + (f1^2)/(dn1 + dn2)
    return(est)
}

# @title Rarefaction curve
# compute rarefaction data
#
# function compute
# @param x a numeric vector
# @return a rarefied vector
# @export
rarefyDT <- function(x) {
    if (missing(x)) stop("x is missing, a vector of integers is expected.")
    lib.size <- sum(x)
    breakp <- lib.size * 0.1
    if (breakp < 10) {
        steps <- pretty(seq_len(lib.size))
        steps[length(steps)] <- lib.size
        } else {
            step1 <- pretty(seq_len(breakp), n = 10)
            step2 <- pretty(breakp:sum(x), n = 10)
            step2[length(step2)] <- lib.size
            steps <- c(step1, step2[step2>0])
            }
    xx <- vegan::rarefy(x, sample=steps)
    output <- data.frame(x = attr(xx, "Subsample"), y=as.double(xx))
    output$x <- as.double(output$x)
    return(output)
}

# @title Calculation of rarefaction values
#
# @description This function computes rarefaction values for each sample.
#
# Rarefaction is a measure of species richness. It assesses the number of
# observed clones for sets of N of sequences in a sample by randomly
# re-sampling a number of sequences multiple times and calculating the mean
# number of observed clones.
#
# @param x an object of class \code{\linkS4class{RepSeqExperiment}}
# @return a data.table of values that can be represented graphically as
# rarefaction curves. The function \code{\link{plotRarefaction}} can be used
# for that purpose.
# @export
# @keywords internal
# @examples
#
# data(RepSeqData)
#
# rarefactionTab(x = RepSeqData)
#
# rarefactionTab <- function(x) {
#     if (missing(x)) stop("x is required. An object of class RepSeqExperiment is expected.")
#     if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected")
#     cts <- data.table::copy(assay(x))
#     raretab <- cts[, round(rarefyDT(count),2), by = sample_id]
#     return(raretab)
# }

#' @title Down-sampling of repertoires
#'
#' @description This function downsamples all repertoires in the dataset.
#'
#' Users can choose the value to which all the samples are down-sampled.
#' If not specified, the lowest number of sequences across all samples within
#' the dataset will be used.
#'
#' This strategy can be applied when studying different cell subsets with
#' significant differences in their repertoire sizes.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sample.size an integer indicating the desired down-sampled size.
#' The default is the smallest repertoire size among all samples of the dataset.
#' @param rngseed a integer used as seed for a reproducible result.
#' @param replace a boolean indicating if the resampling should be performed
#' with or without replacement. Default is TRUE.
#' @param verbose a boolean indicating whether or not to show the details of
#' every computation step within the function. Default is TRUE.
#' @return a new \code{\linkS4class{RepSeqExperiment}}  with the downsized data.
#'
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' RepSeqData_ds<- sampleRepSeqExp(x = RepSeqData,
#'                                 rngseed = 1234,
#'                                 replace = TRUE)
#'
#' RepSeqData_ds<- sampleRepSeqExp(x = RepSeqData,
#'                                 rngseed = FALSE,
#'                                 replace = FALSE)
#'
sampleRepSeqExp <- function(x, sample.size = min(mData(x)$nSequences), rngseed = FALSE, replace = TRUE, verbose = TRUE) {
    if (missing(x)) stop("x is missing, an object of class RepSeqExperiment is expected")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected")
    sampleinfo <- mData(x) %>% dplyr::select(-chao1 ,-iChao)

    if (as(rngseed, "logical")) {
        set.seed(rngseed)
        if (verbose) {
            cat("`set.seed(", rngseed, ")` was applied for reproducibility of the random subsampling. Please save this seed for data reproduction \n")
            cat("You can try `set.seed(", rngseed, "); .Random.seed` for the full vector.\n")
        }
    } else if (verbose) {
        cat("You set `rngseed` to FALSE. Make sure you've set & saved\n",
            " the random seed of your session for reproducibility.\n",
            "See `?set.seed`\n")
    }
    if (length(sample.size) > 1) {
        warning("`sample.size` had more than one value. ", "Using only the first value. \n ... \n")
        sample.size <- sample.size[1]
    }
    if (sample.size <= 0) {
        stop("sample.size = ", sample.size, " is less than or equal to zero. ", "Need positive sample size to work.")
    }
    if (sample.size > max(sampleinfo$nSequences)) {
        stop("sample.size = ", sample.size, " is larger than the largest library size.")
    }
    if (verbose) {
        message("Down-sampling to ", sample.size, " sequences...")
    }
    if (min(sampleinfo$nSequences) < sample.size) {
        lib.drop <- rownames(sampleinfo)[which(sampleinfo$nSequences < sample.size)]
        if (verbose) {
            message(length(lib.drop), " samples were removed because they contained fewer sequences than the fixed sample.size.")
            message("The first 5 removed samples are: \n")
            message(lib.drop[seq_len(min(5, length(lib.drop)))])
            message("...")
        }
        x <- dropSamples(x, lib.drop)
    }
    cts <- data.table::copy(assay(x))
    if (replace) {
        cts[, count := {rareout = numeric(length(count));
                        tmp1 = sample(seq_len(length(count)), sample.size, replace = TRUE, prob = count/sum(count));
                        tmp2 = table(tmp1);
                        rareout[as(names(tmp2), "integer")] <- tmp2;
                        list(rarevec=rareout)}, by = sample_id]
    } else {
        cts[,count := {rareout = numeric(length(count));
                   tmp1 = rep(seq_len(.N), count);
                   tmp2 = sample(tmp1, 100, replace = FALSE);
                   tmp3 = table(tmp2);
                   rareout[as(names(tmp3), "integer")] <- tmp3;
                   list(rarevec=rareout)}, by = sample_id]
    }
    cts <- cts[count>0]
    stats <- data.frame(cts[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone" , "ntClone"), by = "sample_id"], row.names = 1)
    sampleinfo <- data.frame(merge(sampleinfo[, setdiff(colnames(sampleinfo), colnames(stats))], stats, by = 0, sort=FALSE), row.names = 1)

    
	x.hist <- paste("A down-sampling to", sample.size , "with replacement set to", replace, "was performed")
		if (verbose) message("Creating a RepSeqExperiment object...\n")
	out <- methods::new("RepSeqExperiment",
	                      assayData = cts,
	                      metaData = sampleinfo,
	                      otherData = oData(x),
	                      History = rbind(History(x), x.hist))
    cat("Done.\n")
    return(out)
}


#' @title Shannon normalization of repertoires
#'
#' @description This function allows to normalize repertoires using the Shannon
#' entropy as a threshold, thus eliminating rare sequences.
#'
#' It is described in (Chaara et al., 2018) and is used to eliminate
#' “uninformative” sequences resulting from experimental noise.
#'
#' It is particularly efficient when applied on small samples as it corrects
#' altered count distributions caused by a high-sequencing depth.
#'
#' Sequences that are eliminated with this method are stored in the \code{oData}
#' slot.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @return a new \code{\linkS4class{RepSeqExperiment}} object with the
#' normalized data.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' RepSeqData_sh<- ShannonNorm(x = RepSeqData)

ShannonNorm <- function(x) {
  if (missing(x)) stop("x is missing, an object of class RepSeqExperiment is expected")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected")
  sampleinfo <- mData(x) %>% dplyr::select(-chao1 ,-iChao)
  exp_shannon <- NULL
  row_number <- NULL
  out <- copy(assay(x))[, .(count = sum(count)), by=c("sample_id", "ntClone")]
  out[, exp_shannon := .renyiCal(count, 1, hill = TRUE), by = "sample_id"]

  keep <- out[order(-count), .SD[data.table::frank(-count, ties.method = "first") <= exp_shannon], by = sample_id]

  res <- copy(assay(x))[keep, on = c("ntClone", "sample_id", "count")][, c("exp_shannon") := NULL]
  # res[, sample_id := factor(sample_id, levels = rownames(sampleinfo))]
  res <- res[order(sample_id)]
  
  filtered <- out[!keep, on = .(sample_id, ntClone)]
  
  stats <- data.frame(res[, c(.(nSequences = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("V", "J", "VJ", "ntCDR3", "aaCDR3", "aaClone" , "ntClone"), by = "sample_id"], row.names = 1)
  sampleinfo <- data.frame(merge(sampleinfo[, setdiff(colnames(sampleinfo), colnames(stats))], stats, by = 0, sort=FALSE), row.names = 1)
  
  x.hist <- "A Shannon-based normalization was performed"

  message("Creating a RepSeqExperiment object...\n")
  out <- methods::new("RepSeqExperiment",
                      assayData = res,
                      metaData = sampleinfo,
                      otherData = c(oData(x), normalisation=list(filtered[,c(1:3)])),
                      History = rbind(History(x), x.hist))
  cat("Done.\n")
  return(out)
}


