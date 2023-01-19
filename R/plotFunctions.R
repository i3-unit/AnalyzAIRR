utils::globalVariables(c("..adj.rr.label..", "Group", "group", "method", "ste", "value", "y","grp","variable","stats","adj.rr.label","p.value.label", "where", "to", "CELL_META", "padj", "log2FoldChange", "BHpvalue", "rn"))

#' @title  Visualization of the VJ combination usage
#'
#' @description  This function calculates the count or frequency of each combination by taking into account the weight of the chosen repertoire level: either "clone" or "clonotype".
#'
#' It offers two types of visualization of the calculated VJ usage in a given sample.
#'
#' - For prop=1, a heatmap with all the VJ combinations is plotted, with the V genes as columns and the J genes as rows.
#'
#' - For prop!=1, a circos plot is drawn showing the top proportions of VJ combinations specified in the prop parameter.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param scale a character specifying whether to plot the VJ usage in "count" or "frequency".
#' @param level a character specifying the level of the repertoire to be taken into account when calculate VJ usages. Should be one of "clone" or "clonotype".
#' @param prop a numeric indicating the proportions of top VJ combinations to be plotted. It ranges from 0 to 1.
#' @examples
#'
#' data(RepSeqData)
#' snames <- rownames(mData(RepSeqData))
#'
#' plotVJusage(x = RepSeqData, sampleName = NULL, scale = "count", level="clone", prop=1)
#'
#'
#' @export
#'
plotVJusage <- function(x, sampleName = NULL, scale = c("count", "frequency"),
                        level = c("clone","clonotype"),
                        prop=1) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  sNames <- rownames(mData(x))
  levelChoice <- match.arg(level)
  if (is.null(sampleName) || sampleName == "") {
    index <- sNames[1]
    cat("Plot the first sample in the dataset:", index)
  } else {
    index <- sampleName
    if (length(sampleName) > 1) {
      index <- sampleName[1]
    }
    if (is.na(match(index, sNames))) stop("Sample ", index, " not found in x.")
  }
  cts <- data.table::copy(assay(x))
  cts_b <- cts[sample_id == index, ][, .(count = sum(count)), by=c("sample_id", levelChoice, "V","J")]
  tmp <- data.table::dcast(cts_b, J~V, value.var = "count", fun.aggregate = sum)
  data2plot <- data.frame(tmp, row.names = 1)
  if(scale == "frequency"){
    data2plot <- prop.table(data2plot)
  }
 if (0!=prop & prop<1){
  data2plot$to<- rownames(data2plot)
  data2plot_m<- data2plot %>%
                reshape2::melt() %>%
                dplyr::filter(value!=0)%>%
                  dplyr::arrange(desc(value) ) %>%
                  dplyr::slice(seq_len(floor(prop*nrow(.)))) %>%
                dplyr::relocate(to, .after = variable)

  p<- circlize::chordDiagram(data2plot_m, annotationTrack =  "grid", annotationTrackHeight = 0.03,
               preAllocateTracks = list(track.height=0.2), big.gap = 20)
  p1<- circlize::circos.track(track.index=1, panel.fun=function(x,y){
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], circlize::CELL_META$sector.index,
                facing="clockwise", niceFacing = TRUE, adj=c(0,0.5), cex=0.45)
  }, bg.border=NA)
 } else {

  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    p1<-ComplexHeatmap::pheatmap(data2plot,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                          cluster_rows = FALSE, cluster_cols = FALSE, name = " ",
                          silent =TRUE, angle_col="90",fontsize =4, scale="column")
    return(p1)
    
  }
   }

}


#' @title Visualization of CDR3 spectratyping
#'
#' @description CDR3 spectratyping allows large-scale analysis of the repertoire diversity based on the CDR3 amino acid length usage.
#'
#' With this technique, a polyclonal repertoire is represented by an eight-peak bell-shaped profile, that is disturbed following an immune response.
#'
#' This function plots a histogram of the CDR3 length distribution in a single sample.
#' It is also possible to add the information on the V gene usage per CDR3 length.
#' To do so, the proportion of the most used V genes within the repertoire can be specified with the prop parameter, thus showing their relative proportions within each CDR3 length.
#'
#' Moreover, this function integrates a statistical method that can be applied on the spectratypes in order to compare peak distributions across different groups of samples called ISEApeaks (Colette and Six, 2002).
#' The strategy can be applied using the \code{\link{perturbationScore}} function.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param scale a character specifying whether to plot the length usage in "count" or "frequency".
#' @param prop a numeric indicating the proportions of top V genes to be plotted. It ranges from 0 to 1.
#' @export
#' @seealso \code{\link{perturbationScore}}
#' @examples
#'
#' data(RepSeqData)
#' snames <- rownames(mData(RepSeqData))
#' 
#' plotSpectratyping(x = RepSeqData, sampleName = NULL, scale = "count", prop=1)
#'
#' plotSpectratyping(x = RepSeqData, sampleName = NULL, scale = "frequency", prop=0)
#'
plotSpectratyping <- function(x, sampleName = NULL,
                              scale = c("count", "frequency"),
                              prop=0) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  CDR3length=pep=sample_id=CDR3aa=N=frequency <- NULL
  sNames <- rownames(mData(x))
  if (is.null(sampleName)) {
    index <- sNames[1]
    cat("Plot the first sample in the dataset:", index)
  } else {
    index <- sampleName
    if (length(sampleName) > 1) {
      index <- sampleName[1]
    }
    if (is.na(match(index, sNames))) {
      stop("Sample ", index, " not found in x.")
    }
  }
  scl <- match.arg(scale)

  cts <- data.table::copy(assay(x))
  legend_v<- cts[sample_id == index,][, .(.N), by = .(V) ]

  if (prop == 0){
    legend_v<- NULL

  } else {
    legend_v<- legend_v %>%
      dplyr::arrange(desc(N) ) %>%
      dplyr::slice(seq_len(floor(prop*length(unique(cts$V)))))
  }

  if (scl == "count") {
    data2plot <- cts[sample_id == index,][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)]
    data2plot$V<- ifelse(data2plot$V %in% legend_v$V,  data2plot$V , "Other")
     p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = N, fill = V), color="black")

  }else if (scl == "frequency") {
    data2plot <- cts[sample_id == index,][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,frequency := prop.table(N)]
    data2plot$V<- ifelse(data2plot$V %in% legend_v$V,  data2plot$V , "Other")
     p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = frequency, fill = V), color="black")

  }
  p1 <- p + ggplot2::geom_bar(stat = "identity", position="stack", color="black", size=0.1) +
    ggplot2::xlab("CDR3 length (aa)")+ggplot2::ylab(paste(scl))+
    ggplot2::scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(data2plot$V))))+
    theme_RepSeq()+
    ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size=8))

  return(p1)
}

#' @title Visualization of the Renyi index
#'
#' @description This function plots the Renyi values at any repertoire level for all the samples in the dataset.
#' The alpha values for which the Renyi is to be estimated can be personalized, thus allowing to focus on certain indices such as the Shannon index for alpha=1 or the Simpson index for alpha=2.
#' For instance, the most frequent sequence is attributed a rank of 1, and its relative abundance is plotted on the y-axis.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param alpha a numerical vector specifying the alpha values to compute. If not specified, the following values are estimated: c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf).
#' @param level a character specifying the level of the repertoire to be taken into account when calculating VJ usages. Should be one of clone","clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param grouped a boolean indicating whether or not the mean and se of samples belonging to the same experimental group specified in the ColorBy parameter should be calculated. Default is FALSE. The area under the curve (AUC) is calculated for each sample and compared between chosen groups using a Wilcoxon test with Holm's correction.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#'
#' @details The Renyi index is a generalization of the Shannon index. It represents the distribution of clonal expansions
#' as a function of the parameter alpha. At alpha=0, it equally considers all species including the rare ones,
#' whereas it up-weighs the abundant species with an increasing value of alpha. Alpha =1 is an approximation of the Shannon index; 
#' alpha = 2 corresponds to the Simpson index and alpha=Inf corresponds to the Berger-Parker index. The latter highlights the highest clonal expansion in a repertoire.
#'
#' @export
#' @examples
#'
#' data(RepSeqData)
#' plotRenyiIndex(x = RepSeqData, level = "V", colorBy = "sex")
#'
#' plotRenyiIndex(x = RepSeqData, level = "J", colorBy = "sex", grouped=TRUE)
#'
#'
plotRenyiIndex <- function(x, alpha = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf),
                              level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"),
                              colorBy=NULL, grouped=FALSE, label_colors=NULL) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (length(alpha) < 2) stop("At least 2 alpha values are needed.")
  if(is.null(colorBy)) stop("need to specify a group column from mData")

  if (is.null(label_colors)) {
    label_colors = plotColors(x) }

  variable=value <- NULL
  sdata <- mData(x)
  sNames <- rownames(sdata)
  levelChoice <- match.arg(level)

  tmp <- renyiIndex(x, alpha = alpha, level=levelChoice)

  data2plot <- data.table::melt(data = tmp, id.vars = "variable", measure.vars = sNames, variable.name = "sample_id")
  data2plot[, Group := lapply(.SD, function(x) sdata[x, colorBy]), .SDcols = "sample_id"]
  colnames(data2plot)[1]<-"variable2"
  
  
  if (grouped) {
        aucs <- data2plot %>%
          dplyr::group_by(sample_id) %>%
          dplyr::mutate(AUC=MESS::auc(x=as.numeric(variable2), y=value))
        
        se<- function(x) sqrt(var(x)/length(x))
        data2plot <-  data2plot %>% dplyr::group_by(Group, variable2) %>% dplyr::mutate(ste=se(value))
        data2plot <- data2plot %>% dplyr::group_by(Group, variable2, ste) %>% dplyr::summarize(mean=mean(value))
    
        data2plot<- setDT(data2plot)
        data2plot[, `:=`(alpha, as.numeric(as.character(variable2)))]
        
        auc_test <- data.frame(aucs) %>%
          rstatix::wilcox_test(formula = AUC ~ Group) %>%
          rstatix::adjust_pvalue(method="holm")
        
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = as.numeric(as.character(variable2)), y = mean)) +
            ggplot2::geom_line(ggplot2::aes(group = Group, color=Group), size = .8) +
            ggplot2::geom_point(ggplot2::aes( color=Group),  shape=21, size=2)+
            ggplot2::labs(subtitle=gsub('p', 'p.adj', get_test_label(auc_test, detailed=FALSE, p.col = "p.adj", type="text")))+
            ggplot2::geom_ribbon(
            ggplot2::aes(ymin=mean-ste, ymax=mean+ste,  fill=Group),
            alpha = 0.3,colour=NA )+
            ggplot2::xlab("alpha")  + ggplot2::ylab("Renyi's Entropy") +
            ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
            ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
            theme_RepSeq()+ggplot2::theme(legend.position="right", plot.subtitle = element_text(hjust=0.90, vjust=-10))
    } else {
      p<- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = variable2, y = value, color=Group)) +
        ggplot2::geom_line(ggplot2::aes(group = sample_id), size = .8) +
        ggplot2::geom_point( shape=21, size=2)+
        ggplot2::xlab("alpha")+
        ggplot2::ylab("Renyi's Entropy") +
        ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
        ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
        theme_RepSeq()+ggplot2::theme(legend.position="right")
    }
  return(p)
}

#' @title Visualization of the V or J gene usage
#'
#' @description This function plots the V or J gene usages for a given sample.
#'
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the gene usages. Should be one of clone" or "clonotype".
#' @param scale a character specifying whether to plot the gene usage in "count" or "frequency".
#' @export
#' @examples
#'
#' data(RepSeqData)
#' snames <- rownames(mData(RepSeqData))
#'
#' plotIndGeneUsage(x = RepSeqData, level = "V", sampleName = snames[1], scale = "count")
#'
#' plotIndGeneUsage(x = RepSeqData, level = "V", sampleName = snames[2], scale = "frequency")

plotIndGeneUsage <- function(x,  sampleName = NULL, level = c("V", "J"), scale = c("count", "frequency")) {
  frequency <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  if (is.null(sampleName)) {
    sName <- rownames(mData(x))[1]
    cat("Plot the first sample in the dataset:", sName)
  } else sName <- sampleName

  levelChoice <- match.arg(level)
  cts <- data.table::copy(assay(x))
  data2plot <- cts[sample_id == sName, lapply(.SD, sum), by = levelChoice, .SDcols = "count"]

  if (scale == "count") {
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x=levelChoice, y = "count", fill = levelChoice))
  }
  if (scale == "frequency") {
    data2plot[, frequency := prop.table(count)]
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x=levelChoice, y = "frequency", fill = levelChoice))
  }
  cols<-RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category =="qual" & RColorBrewer::brewer.pal.info$colorblind == "TRUE", ]
  p1 <- p + ggplot2::geom_bar(stat = "identity", show.legend=FALSE) +
    ggplot2::scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nrow(unique(data2plot[, 1]))))+
    theme_RepSeq()+
    ggplot2::xlab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,size=8),
                   axis.text.y = ggplot2::element_text(size=8))


  return(p1)
}

#' @title Compare V or J gene distributions
#'
#' @description This function compares the V or J gene usages between given groups.
#'
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the gene usages. Should be one of clone" or "clonotype".
#' @param scale a character specifying whether to plot the gene usage in "count" or "frequency".
#' @param groupBy character indicating the groups to be compared. If not specified, no comparative statistical tests will be performed, and calculated values for each sample_id will be represented.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotGeneUsage(x = RepSeqData, level = "J", scale = "count", groupBy = "cell_subset")
#'
#'
#'
plotGeneUsage <- function(x, level = c("V", "J"), scale = c("count", "frequency"), 
                             groupBy = NULL, label_colors = NULL) {
  frequency <- NULL
  sdata<-mData(x)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(groupBy)) stop("a group column from the metadata is expected.")
  if (is.null(label_colors)) {
    label_colors = plotColors(x)
  }
  
  levelChoice <- match.arg(level)
  scaleChoice <- match.arg(scale)
  
  cts <- data.table::copy(assay(x))
  data2plot <- cts[, lapply(.SD, sum), by = c(levelChoice, "sample_id"), .SDcols = "count"]
  data2plot <- data2plot %>% dplyr::group_by(sample_id) %>% dplyr::mutate(frequency=count/sum(count)) %>% setDT()

  data_exp <- data2plot %>% tidyr::expand(sample_id,get(levelChoice))
  colnames(data_exp)[2] <- paste(levelChoice)
  data2plot <- dplyr::left_join(data_exp,data2plot, by=c("sample_id",paste(levelChoice)))
  data2plot <- data2plot %>% 
    replace(is.na(.),0) %>% setDT()
  
  if(length(groupBy)==1){
    data2plot[, Group := lapply(.SD, function(x) sdata[x, groupBy]), .SDcols = "sample_id"]
  } else if(length(groupBy)==2){
    data2plot[, Group := lapply(.SD, function(x) sdata[x, groupBy[[1]]]), .SDcols = "sample_id"]
    data2plot[, Group2 := lapply(.SD, function(x) sdata[x, groupBy[[2]]]), .SDcols = "sample_id"]
  } else if(length(groupBy)==3){
    data2plot[, Group := lapply(.SD, function(x) sdata[x, groupBy[[1]]]), .SDcols = "sample_id"]
    data2plot[, Group2 := lapply(.SD, function(x) sdata[x, groupBy[[2]]]), .SDcols = "sample_id"]
    data2plot[, Group3 := lapply(.SD, function(x) sdata[x, groupBy[[3]]]), .SDcols = "sample_id"]
  }
  
    stat.test <- data2plot %>%
      dplyr::select(all_of(levelChoice), sample_id, all_of(scaleChoice), tidyr::starts_with("Gr")) %>%
      dplyr::group_by(get(levelChoice)) %>%
      rstatix::wilcox_test(formula=as.formula(paste(paste(scaleChoice),"Group",sep="~"))) %>%
      rstatix::adjust_pvalue() %>%
      rstatix::add_significance() %>%
      rstatix::add_y_position() 
  
  data2plot.summary<- data2plot %>%
         dplyr::group_by_at(vars(tidyr::starts_with("Gr"), paste(levelChoice))) %>%
          dplyr::summarise(
            sd = sd(get(scaleChoice), na.rm = TRUE),
            mean = mean(get(scaleChoice)),
            n = n(),
            se = sd / sqrt(n)
          ) %>%
    dplyr::rename(levelChoice=paste(levelChoice))
  
  p <-ggplot2::ggplot(data2plot.summary, 
                      aes(x = levelChoice, y = mean ))+
       ggplot2::geom_bar(aes(fill=Group), stat="identity", position = position_dodge()) +
      ggplot2::geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=Group), width=.2,position=position_dodge(.9))+
    {if(length(groupBy)==2)list(ggplot2::facet_grid(Group2~.))} +
    {if(length(groupBy)==3)list(ggplot2::facet_grid(Group2~Group3))} +
    ggplot2::scale_fill_manual(values=label_colors[[groupBy[[1]]]])+
    theme_RepSeq()+
    ggplot2::xlab("") +
    ggplot2::ylab(paste(scaleChoice)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,size=8),
                   axis.text.y = ggplot2::element_text(size=8),
                   legend.position = "top")+
    ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                               tip.length = 0,x="get(levelChoice)", hide.ns = TRUE)
  
  return(p)
}



#' @title Visualization of the clonal distribution per count interval in a single ample
#'
#' @description This function plots two histograms of the clonal distribution per a set of fixed count intervals in a given sample. Intervals are: 1, ]1, 10], ]10, 100], ]100, 1000], ]1000, 10000] and ]10000, Inf].
#'
#' The upper plot calculates the proportion of each interval in the whole repertoire, whereas the lower plot shows the cumulative frequency of the sequences within each interval.
#'
#' This allows a global view of the repertoire fraction contributing the most to the repertoire. For instance, top sequences belonging to the highest interval often constitute a low fraction of the whole repertoire but contribute more significantly in terms of cumulative frequency in view of their high occurrence.
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the clonal distribution. Should be one of clone","clonotype", "CDR3nt" or "CDR3aa".
#' @export
#' @examples
#'
#' data(RepSeqData)
#' snames <- rownames(mData(RepSeqData))
#'
#' plotIndCountIntervals(x = RepSeqData, level="CDR3aa",sampleName = snames[1])
#'
#' plotIndCountIntervals(x = RepSeqData, level="clone", sampleName = NULL)
#'
plotIndCountIntervals <- function(x, sampleName = NULL, level = c("clone","clonotype", "CDR3nt","CDR3aa")){
  interval=percent <- NULL
  levelChoice<- match.arg(level)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(sampleName)) {
    sName <- rownames(mData(x))[1]
    cat("Plot the first sample in the dataset:", sName)
  } else sName <- sampleName
  data2plot <- data.table::copy(assay(x))
  f <- function(x){
    if(x == 1) "1"
    else if(x <= 10) "]1, 10]"
    else if(x <= 100) "]10, 100]"
    else if(x <= 1000) "]100, 1000]"
    else if(x <= 10000) "]1000, 10000]"
    else "]10000, Inf]"
  }
  ff<- data.frame(interval=c("1","]1, 10]","]10, 100]", "]100, 1000]","]1000, 10000]","]10000, Inf]"))
  data2plot <- data2plot[sample_id == sName, lapply(.SD, sum), by = levelChoice, .SDcols = "count"][,interval := unlist(lapply(count, f))]

  data2plot_b <- data2plot %>% dplyr::group_by(interval) %>% dplyr::summarize(sum=dplyr::n()) %>% dplyr::mutate(freq=sum/sum(sum))
  data2plot<- data2plot[, lapply(.SD, sum), by = interval, .SDcols = "count"][,percent := prop.table(count)]
  data2plot<- merge(data2plot[,c(1,3)], data2plot_b[,c(1,3)], by="interval")
  if(nrow(data2plot)<6){
    rown= which(!ff$interval %in% data2plot$interval)
    add<- data.frame(interval=ff[rown,], percent=0, freq=0)
    data2plot<- rbind(data2plot, add)
  }
  breaks <- unique(data2plot[, interval])
  plotBreaks <- breaks[order(nchar(breaks), breaks)]
  data2plot<- reshape2::melt(data2plot, id.vars="interval")

  p1 <- ggplot2::ggplot(data = data2plot[data2plot$variable == "percent",], ggplot2::aes(x = interval, y =value) , fill = "gray", alpha=.7) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), alpha=.8) +
    ggplot2::scale_x_discrete(limits=plotBreaks) +
    ggplot2::ylim(0, 1) +
    ggplot2::geom_text(ggplot2::aes(y = value, label = paste(100*round(value,3), "%")),  size=3, position = ggplot2::position_dodge(.8), vjust=-1) +
    ggplot2::xlab("")+
    theme_RepSeq()+
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                    axis.text.x=ggplot2::element_text( size=8))+
    ggplot2::ggtitle("Cumulative frequency")

  p2 <- ggplot2::ggplot(data = data2plot[data2plot$variable == "freq",], ggplot2::aes(x = interval, y =value ) , fill = "gray", alpha=.7) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), alpha=.8) +
    ggplot2::scale_x_discrete(limits=plotBreaks) +
    ggplot2::ylim(0, 1) +
    ggplot2::geom_text(ggplot2::aes(y = value, label = paste(100*round(value,3), "%")),  size=3, position = ggplot2::position_dodge(.8), vjust=-1) +
    ggplot2::xlab("")+
    theme_RepSeq()+
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text( size=8))+
    ggplot2::ggtitle("Distribution")

  g <- gridExtra::arrangeGrob(p2, p1 ,nrow=2)
  grid::grid.draw(g)
}


#' @title Visualization of the CDR3 spectratyping per V gene
#'
#' @description CDR3 spectratyping allows large-scale analysis of the repertoire diversity based on the CDR3 amino acid length usage.
#' With this technique, a polyclonal repertoire is represented by an eight-peak bell-shaped profile that is disturbed following an immune response.
#'
#' This function plots a histogram of the CDR3 length distribution per V gene in a single sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param scale a character specifying whether to plot the length usage in "count" or "frequency".
#' @param prop a numeric indicating the proportions of top V genes to be plotted. It ranges from 0 to 1.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' snames <- rownames(mData(RepSeqData))
#'
#' plotSpectratypingV(x = RepSeqData, sampleName = snames[1], scale = "count",prop=0.05)
#'
#' plotSpectratypingV(x = RepSeqData, sampleName = snames[1], scale = "frequency")
#'
#'
plotSpectratypingV <- function(x, sampleName = NULL, scale = c("count", "frequency"),
                               prop=0) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  frequency=N=CDR3length=pep=sample_id=CDR3aa <- NULL
  sNames <- rownames(mData(x))
  if (is.null(sampleName)) {
    index <- sNames[1]
    cat("Plot the first sample in the dataset:", index)
  } else {
    index <- sampleName
    if (length(sampleName) > 1) {
      index <- sampleName[1]
    }
    if (is.na(match(index, sNames))) {
      stop("Sample ", index, " not found in x.")
    }
  }
  scl <- match.arg(scale)
  cts <- data.table::copy(assay(x))

  legend_v<- cts[sample_id == index,][, .(.N), by = .(V) ]
  if (prop == 0){
    legend_v<- legend_v

  } else {
    legend_v<- legend_v %>%
      dplyr::arrange(desc(N) ) %>%
      dplyr::slice(seq_len(ceiling(prop*length(unique(cts$V)))))
  }

  if (scl == "count"){
    tmp <- cts[sample_id == index,][, CDR3length := nchar(CDR3aa)]
    data2plot <- tmp[, .(.N), by = .(V, CDR3length)]
    data2plot<- data2plot[data2plot$V %in% legend_v$V,]
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = N))
  }
  if (scl == "frequency"){
    data2plot <- cts[sample_id == index,][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,percent := prop.table(N), by = V]
    data2plot<- data2plot[data2plot$V %in% legend_v$V,]

    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = percent))
  }

  p <- p + ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_continuous(breaks = data2plot[, unique(CDR3length)]) +
    # ggplot2::labs(subtitle = paste("CDR3 length distribution per V gene in",paste(index))) +
    ggplot2::xlab("CDR3 length (aa)")+
    ggplot2::ylab(paste(scl))+
    ggplot2::facet_wrap(~ factor(V, levels = naturalsort::naturalsort(unique(V))), ncol = 6) +
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = ggplot2::element_text(size=9),
                   axis.title = ggplot2::element_text(size=10),
                   strip.text.x = ggplot2::element_text(size = 6),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size=6))
  return(p)
}

#' @title Visualization of repertoire dissimilarities
#'
#' @description This function assesses the repertoire dissimilarities by plotting a heatmap of a squared distance matrix computed between samples using a specific dissimilarity method.
#'
#'
#' It calculates a list of dissimilarity indices, each taking into account different parameters. The proposed methods include:
#' 
#'  The Jaccard similarity: a measure of similarity between sample sets defined as the size of the intersection divided by the size of the union of the sample sets.
#'
#'  The Morisita-Horn similarity: a measure of similarity that tends to be over-sensitive to abundant species.
#'
#' The function also performs a hierarchical clustering on the calculated distance scores.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the indices are computed. Should be one of "clone","clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param method a character specifying the distance method to be computed. Should be one of the following: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis."
#' @param clustering a character specifying the clustering method to be used.
#' @param binary a boolean indicating whether or not to transform the data into a presence/absence data. Default is FALSE
#' @param label_colors a list of colors for each factor column in metaData. See \code{\link{plotColors}}. If NULL, default colors are used.
#'
#' @details Details on the calculated indices as well as the clustering methods can be found in the vegan package: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotDissimilarityMatrix(x = RepSeqData, level = "V", method = "euclidean")
#'
plotDissimilarityMatrix <- function(x, level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"),
                                    method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
                                               "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
                                               "binomial", "chao", "cao", "mahalanobis"),
                                    clustering=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty" ,
                                                "median", "centroid" ),
                                    binary = FALSE,
                                    label_colors=NULL) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  variable <- NULL
  levelChoice <- match.arg(level)
  methodChoice <- match.arg(method)
  clust<- match.arg(clustering)
  cols <- c("sample_id", levelChoice, "count")
  tmp <- data.table::copy(assay(x))[, ..cols]
  sdata <- mData(x)
  sNames <- rownames(sdata)
  groups <- sdata[, unlist(lapply(sdata, is.factor)), drop = FALSE]

  if (is.null(label_colors)) {
    label_colors= plotColors(x)
  }

  dat <- data.table::dcast(data = tmp, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)

  simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag = TRUE, upper = TRUE, binary = binary), .SDcols=sNames]

  p <- ComplexHeatmap::pheatmap(as.matrix(simmat),
                     cluster_rows = TRUE, cluster_cols = TRUE, name = " ",  col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     treeheight_row = 0L, clustering_distance_rows = simmat, clustering_distance_cols = simmat,
                     annotation_col=groups, show_colnames=FALSE, labels_col = sNames, annotation_colors = label_colors,
                     show_rownames=FALSE, clustering_method = clust, silent = FALSE, fontsize =7)
  return(p)
}



#' @title Visualization of the clonal distribution as a function of the rank
#'
#' @description This function plots the clonal distribution, at any repertoire level, as a function of the occurrence rank within each sample.
#'
#' For instance, the most frequent sequence is attributed a rank of 1, and its relative abundance is plotted on the y-axis.
#'  
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account to plot the clonal distribution. Should be one of "clone", "clonotype", "CDR3nt" or "CDR3aa".
#' @param scale a character specifying whether to plot the clonal abundance in "count" or "frequency".
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param grouped a boolean indicating whether or not the mean and se of samples belonging to the same experimental group specified in the ColorBy parameter should be computed. Default is FALSE. The area under the curve (AUC) is calculated for each sample and compared between chosen groups using a Wilcoxon test with Holm's correction.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @export
#' @examples
#'
#' data("RepSeqData")
#'
#' plotRankDistrib(x = RepSeqData, level="clone",colorBy = "sex", scale = "count")
#'
#' plotRankDistrib(x = RepSeqData,level="clonotype", colorBy = "cell_subset", scale = "frequency")
#'
#'
plotRankDistrib <- function(x, level = c("clone","clonotype", "CDR3nt","CDR3aa"),
                             scale=c("count","frequency"),
                             colorBy = NULL, grouped=FALSE,
                             label_colors=NULL){
  group <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if(is.null(colorBy)) stop("need to specify a column from mData. Can be the sample_id column")

  if (is.null(label_colors)) {
    label_colors = plotColors(x)
  }

  scl <- match.arg(scale)
  levelChoice <- match.arg(level)
  sName <- rownames(mData(x))

  sdata <- mData(x)
  counts <- data.table::copy(assay(x))
  counts <- counts[, .(count = sum(count)), by=c("sample_id", levelChoice)][, rank := lapply(.SD, frankv, ties.method = "first", order = -1L), by = "sample_id", .SDcols = "count"]
  counts[, group := lapply(.SD, function(x) sdata[x, colorBy]), .SDcols = "sample_id"]
  
  if (grouped) {
    aucs <- counts %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(AUC=MESS::auc(x=rank, y=count))
    
    se<- function(x) sqrt(var(x)/length(x))

    if (scl == "count"){
      counts[, ste := lapply(.SD, se), by = .( group, rank), .SDcols = "count"]
      counts <- counts[, lapply(.SD, mean), by = .( group,rank, ste), .SDcols = "count"]

    } else if (scl == "frequency"){
      counts <- counts %>% dplyr::group_by(sample_id) %>% dplyr::mutate(count=count/sum(count))
      counts <- data.table::setDT(counts)
      counts[, ste := lapply(.SD, se), by = .( group,rank), .SDcols = "count"]
      counts <- counts[, lapply(.SD, mean), by = .( group,rank, ste), .SDcols = "count"]

    }
    
    auc_test <- data.frame(aucs) %>%
      rstatix::wilcox_test(formula = AUC ~ group) %>%
      rstatix::adjust_pvalue(method="holm")
    
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = rank, y = count, colour = group)) +
    ggplot2::geom_point(shape=21) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = count-ste, ymax = count+ste,  fill=group),
        alpha = 0.3, colour = NA)+
      ggplot2::geom_line()+
      ggplot2::labs(subtitle=gsub('p', 'p.adj', get_test_label(auc_test, detailed=FALSE, p.col = "p.adj", type="text")))+
      ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
      ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
      ggplot2::scale_x_log10()+
      ggplot2::ylab(paste0("mean ", scl))+
      theme_RepSeq()+
      ggplot2::theme(legend.position = "right", plot.subtitle = element_text(hjust=0.90, vjust=-10))

  } else {

    if (scl == "frequency"){
      counts <- counts %>% dplyr::group_by(sample_id) %>% dplyr::mutate(count = count/sum(count))
      counts <- data.table::setDT(counts)
    } else {
      counts <- counts
    }

    p <-  ggplot2::ggplot(counts, ggplot2::aes(x = rank, y = count, colour = group)) +
          ggplot2::geom_point(shape=21) +
          ggplot2::geom_line(ggplot2::aes(group=sample_id))+
          ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
          ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
          ggplot2::scale_x_log10()+
          theme_RepSeq()+
          ggplot2::theme( legend.position = "right")+
          ggplot2::ylab(paste(scl))

  }

  return(p)
}


#' @title Visualization of the inter-repertoire sharing on an Venn diagram
#'
#' @description The repertoire sharing at any level evaluates the degree of convergence between repertoires and experimental conditions.
#'
#' This function plots the number of shared sequences, at any repertoire level, between samples belonging for instance, to a same experimental group.
#' No limitations in the number of sample_ids is imposed, however for visualization purposes, a maximum of 5 sample_ids should be ideally specified.
#'
##' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "clone","clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param sampleNames a vector of character indicating the sample_ids of the repertoires to be analyzed. If not specified, the first three samples in the dataset are analyzed.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' snames <- rownames(mData(RepSeqData))[1:4]
#'
#' plotVenn(x = RepSeqData, level = "V", sampleNames = snames)
#'
#' plotVenn(x = RepSeqData, level = "J", sampleNames = NULL)
#'
plotVenn <- function(x, level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"), sampleNames = NULL) {
  ..sampleNames <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  levelChoice <- match.arg(level)
  sdata <- mData(x)
  label_colors <- plotColors(x)

  if (is.null(sampleNames)) sampleNames <- rownames(sdata)[seq_len(3)]
  if (length(sampleNames) == 1) stop("At least 2 sampleNames must be provided.")

  counts <- data.table::copy(assay(x))[sample_id %in% sampleNames]
  list<-list()
  for ( i in unique(counts$sample_id)){
    counts_i<- counts %>% filter(sample_id==i)
    list[[i]]<- unique(counts_i[[levelChoice]])
  }

 plot <- ggVennDiagram::ggVennDiagram(
   list,  label_alpha = 0,
   label_percent_digit = 1 )+
   ggplot2::scale_fill_distiller(palette = "RdBu")+
   ggplot2::scale_color_manual(values=label_colors$sample_id)
 
 

 
 plot <- ggVennDiagram::plot_venn(list, label_alpha = 0,edge_lty="solid",
                                  show_intersect=F,label="both",label_percent_digit = 1,
                                  set_size=3,edge_size=.5 ,set_color="black",
                                  label_geom="label", label_color="black",
                                  label_size=2.5)+
         ggplot2::scale_color_manual(values=label_colors$sample_id[sampleNames])+
         ggplot2::scale_fill_distiller(palette = "RdBu")+
         ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2))
 
  return(plot)
}


#' @title Visualization of the inter-repertoire sharing on a scatter plot
#'
#' @description This function computes a correlation between a pair of repertoires at any repertoire level.
#' It allows a comparative analysis of the occurrence of a given repertoire level between a pair of samples, which goes beyond the simple calculation of the sharing degree.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be analyzes. Should be one of "clone", "clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param sampleNames a vector of character indicating the two sample_ids of the repertoires to be analyzed. If not specified, the first two samples in the dataset are analyzed.
#' @param scale an integer indicating whether to plot a regular or a logarithmic scale.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' plotScatter(x = RepSeqData,
#'             level = "V",
#'             sampleNames = c("tripod-30-813", "tripod-30-815"),
#'             scale = "log")
#'
plotScatter <- function(x, sampleNames = NULL,
                        level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"),
                        scale = c("frequency", "log")) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (length(sampleNames) != 2) stop("Two sampleNames are required.")
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
  add_freq<- apply(data2plot[,-1], 2, freqs)
  data2plot<- cbind(data2plot[,1], as.data.frame(add_freq))


  sampleNames <- vapply(sampleNames, function(s) {
    if (!grepl("^`", s)) {
      s <- paste("`", s, sep="", collapse="")
    }
    if (!grepl("`$", s)) {
      s <- paste(s, "`", sep="", collapse="")
    }
  }
  , FUN.VALUE = character(1))
  data2plot<- data.frame(data2plot, check.names = FALSE)
  formulas <- y ~ x
  p <- ggplot2::ggplot(data2plot, ggplot2::aes_string(x = sampleNames[1], y = sampleNames[2])) +
    ggplot2::geom_count(size = 1.5, shape=21, alpha=.5) +
    ggplot2::geom_smooth(method="lm", se=FALSE, linetype="dashed", color="red")+
    theme_RepSeq()+
    ggpmisc::stat_poly_eq(ggplot2::aes(label=paste(ggplot2::after_stat(adj.rr.label),ggplot2::after_stat(p.value.label),sep = "~~~~")),
                          formula=formulas,
                          parse = TRUE)+
    ggplot2::theme(legend.position = "right")

    p2 <- p + ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  if (scale == "log") return(p2) else return(p)
}


#' @title Visualization of the diversity indices
#'
#' @description This function plots and compares a chosen diversity index calculated on a selected repertoire level between groups of samples.
#'
#' The calculated indices can be one of the following:
#'
#' - Shannon index: Calculates the proportional abundance of species in a repertoire.
#'
#' - Simpson index: Takes into account the number of species present as well as their abundance. It gives relatively little weight to the rare species and more weight to the frequent ones
#'
#' - Inverse Simpson index: Is the effective number of species that is obtained when the weighted arithmetic mean is used to quantify average proportional abundance of species.
#'
#' - Berger-Parker index: Expresses the proportional importance of the most abundant species. This metric is highly biased by sample size and richness (Berger and Parker 1970).
#'
#' - Gini coefficient: Measures the degree of inequality in a distribution of abundances.
#'
#' - Chao1: Estimates undetected species using the information on the rarest species (the numbers of singletons and doubletons).
#'
#' - Improved Chao1: an extension of Chao1 which uses additional information, namely, the numbers of tripletons and quadrupletons.
#'
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param index a character specifying the diversity index to be estimated. Should be one of "chao1", shannon","invsimpson","simpson" or "gini".
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "clone","clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param groupBy a character indicating one or multiple groups to be compared. A Wilcoxon test is thus performed and adjusted p-values using the Holm method are shown. Colors are attributed to the different groups within the first column, and a facet is applied on the second column. If not specified, no statistical tests will be performed, and calculated values for each sample_id will be represented. 
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotDiversity(x = RepSeqData, level = "V", groupBy = "sex", index="shannon")
#' 
#' plotDiversity(x = RepSeqData, level = "clone", groupBy = c("cell_subset", "sex"), index="shannon")
#'
plotDiversity <- function(x, index=c("chao1","shannon","simpson", "invsimpson","bergerparker", "gini","iChao"),
                          level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"),
                          groupBy = NULL, label_colors=NULL){
  group <- NULL
  group2 <- NULL
  group3 <- NULL
  sdata<-mData(x)
  levelChoice <- match.arg(level)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  if (is.null(label_colors)) {
    label_colors <- plotColors(x)
  }

  diversity <- diversityIndices(x, level=levelChoice)
  diversity_m <- diversity %>% dplyr::select(sample_id, paste(index)) %>% dplyr::rename(method=paste(index))
  if (!is.null(groupBy)){
    if(length(groupBy)==1){
      diversity_m[, group := lapply(.SD, function(x) sdata[x, groupBy] ), .SDcols = "sample_id"]
    } else if(length(groupBy)==2){
      diversity_m[, group := lapply(.SD, function(x) sdata[x, groupBy[[1]]] ), .SDcols = "sample_id"]
      diversity_m[, groupB := lapply(.SD, function(x) sdata[x, groupBy[[2]]] ), .SDcols = "sample_id"]
    } else if(length(groupBy)==3){
      diversity_m[, group := lapply(.SD, function(x) sdata[x, groupBy[[1]]] ), .SDcols = "sample_id"]
      diversity_m[, groupB := lapply(.SD, function(x) sdata[x, groupBy[[2]]] ), .SDcols = "sample_id"]
      diversity_m[, groupC := lapply(.SD, function(x) sdata[x, groupBy[[3]]] ), .SDcols = "sample_id"]
    }

    if(length(groupBy)==1){
      stat.test <- ggpubr::compare_means(data=diversity_m, formula=method ~ group) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(data=diversity_m, formula=method ~ group) 
    } else if(length(groupBy)==2){
      stat.test <- diversity_m %>%
        group_by(groupB) %>%
        ggpubr::compare_means(formula=method ~ group) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(data=diversity_m, formula=method ~ group)
    } else if(length(groupBy)==3){
      stat.test <- diversity_m %>%
        group_by(groupB, groupC) %>%
        ggpubr::compare_means(formula=method ~ group) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(data=diversity_m, formula=method ~ group)
    }
    
    p<- ggplot2::ggplot(diversity_m,ggplot2::aes(x = group, y = method, fill = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_point(shape=21, position = ggplot2::position_dodge(width=.8))+
      {if(length(groupBy)==2)list(ggplot2::facet_grid(~groupB))} +
      {if(length(groupBy)==3)list(ggplot2::facet_grid(groupB~groupC))} +
      ggplot2::xlab("")+
      ggplot2::ylab(paste(index))+
      ggplot2::scale_color_manual(values=label_colors[[groupBy[[1]]]])+
      ggplot2::scale_fill_manual(values=label_colors[[groupBy[[1]]]])+
      theme_RepSeq()+
      ggplot2::theme( legend.position = "none")+
      ggprism::add_pvalue(stat.test,  label = "p.adj.signif", bracket.nudge.y = .01,tip.length = 0, step.increase = 0.1)
    
}  else {

  p<- ggplot2::ggplot(diversity_m,ggplot2::aes(x = sample_id, y = method, fill = sample_id)) +
    ggplot2::geom_bar(stat="identity", size=.5, color="black") +
    ggplot2::xlab("")+
    ggplot2::ylab(paste(index))+
    ggplot2::scale_fill_manual(values=label_colors[["sample_id"]])+
    theme_RepSeq()+
    ggplot2::theme( legend.position = "none",
                    axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

}
  return(p)
}


#' @title Visualization of the rarefaction curves
#'
#' @description This function plots the rarefaction curve for each sample within the dataset.
#'
#' Rarefaction is a measure of species richness. The curves plot the number of clonotypes against the number of sequences in a sample, each being obtained by randomly re-sampling a number of sequences multiple times and representing the mean number of found clonotypes.
#'
#' Calculation can be obtained using the \code{\link{rarefactionTab}} function which computes rarefaction values for each sample.
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @return none
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotRarefaction(x = RepSeqData, colorBy = "sex")
#'
#' plotRarefaction(x = RepSeqData, colorBy = "cell_subset")
#'
#'
plotRarefaction <- function(x, colorBy=NULL, label_colors=NULL){
  if(is.null(colorBy)) stop("need to specify a group column from mData")
  raretab<- rarefactionTab(x)

  if (is.null(label_colors)) {
    label_colors= plotColors(x) 
    }

  sdata <- mData(x)
  raretab[, group := lapply(.SD, function(x) sdata[x, colorBy] ), .SDcols = "sample_id"]

  p <- ggplot2::ggplot(data = raretab, ggplot2::aes(x = x, y = y, fill = group, color = group)) +
    ggplot2::geom_line(ggplot2::aes(group=sample_id)) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      x = "Number of sequences",
      y = "Number of clonotypes") +
    ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
    ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
    ggrepel::geom_text_repel(data=raretab %>% group_by(sample_id) %>% dplyr::slice_max(x),
      nudge_x = -0.2, direction = "y", hjust = "left", ggplot2::aes(label = sample_id)
    ) +
    theme_RepSeq()+
    ggplot2::theme(legend.position = "none")+
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .1))
  

  return(p)
}


#' @title Visualization of the clonal distribution per count interval
#'
#' @description This function plots two histograms of the clonal distribution per a set of count intervals in all the samples within the dataset. Intervals are: 1, ]1, 10], ]10, 100], ]100, 1000], ]1000, 10000] and ]10000, Inf].
#'
#' The plot titled "Distribution" calculates the proportion of each interval in the whole repertoire, whereas the one titled "Cumulative frequency" shows the cumulative frequency of the sequences within each interval.
#'
#' This could allow a global view of the repertoire fraction contributing the most to the repertoire. For instance, top sequences belonging to the highest interval often constitute a low fraction of the whole repertoire but contribute more significantly in terms of cumulative frequency in view of their high occurrence.
#'
#' Samples can be statistically compared in each interval using the \code{groupBy} parameter.
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the clonal distribution. Should be one of clone","clonotype", "CDR3nt" or "CDR3aa".
#' @param groupBy a character indicating one or multiple column names in mData. Colors are thus attributed to the different groups within the first column, and a facet is applied on the second column. Statistical tests are performed between the chosen groups. The chosen column must be of class factor.
#' @param label_colors a list of colors for each variable in groupBy See \code{\link{plotColors}}. If NULL, default colors are used.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotCountIntervals(x = RepSeqData, level="CDR3aa")
#'
#' plotCountIntervals(x = RepSeqData, level="CDR3nt",  groupBy="cell_subset")
#' 
#' plotCountIntervals(x = RepSeqData, level="CDR3nt",  groupBy=c("cell_subset", "sex"))
#'
plotCountIntervals <- function(x, level = c("clone","clonotype", "CDR3nt","CDR3aa"),
                               groupBy=NULL, label_colors=NULL){
  interval=percent <- NULL
  levelChoice<- match.arg(level)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  data2plot <- data.table::copy(assay(x))
  f <- function(x){
    if(x == 1) "1"
    else if(x <= 10) "]1, 10]"
    else if(x <= 100) "]10, 100]"
    else if(x <= 1000) "]100, 1000]"
    else if(x <= 10000) "]1000, 10000]"
    else "]10000, Inf]"
  }

  if (is.null(label_colors)) {
    label_colors = plotColors(x)
  }

  ff<- data.frame(interval=c("1","]1, 10]","]10, 100]", "]100, 1000]","]1000, 10000]","]10000, Inf]"))
  data2plot <- data2plot[, lapply(.SD, sum), by = c(levelChoice, "sample_id"), .SDcols = "count"][,interval := unlist(lapply(count, f)), by="sample_id"]

  data2plot_b <- data2plot %>% dplyr::group_by(interval, sample_id) %>% dplyr::summarize(sum=dplyr::n())
  data2plot_b <- data2plot_b %>% dplyr::group_by( sample_id) %>% dplyr::mutate(freq=sum/sum(sum))
  data2plot <- data2plot[, lapply(.SD, sum), by = c("interval","sample_id"), .SDcols = "count"][,percent := prop.table(count), by="sample_id"]
  data2plot <- merge(data2plot[,c(1,2,4)], data2plot_b[,c(1,2,4)], by=c("interval", "sample_id"))

  data2plot <- setDT(data2plot)
  sdata <- mData(x)
  breaks <- unique(as.data.frame(data2plot)[, "interval"])

  if(!is.null(groupBy)){
    plotBreaks <- c("1","]1, 10]","]10, 100]", "]100, 1000]" ,"]1000, 10000]","]10000, Inf]")
    
    if(length(groupBy)==1){
      data2plot[, grp := lapply(.SD, function(x) sdata[x, groupBy]), .SDcols = "sample_id"]
      data2plot<- reshape2::melt(data2plot, id.vars=c("interval","sample_id","grp"))
      
    } else if(length(groupBy)==2){
      data2plot[, grp := lapply(.SD, function(x) sdata[x, groupBy[[1]]]), .SDcols = "sample_id"]
      data2plot[, grp2 := lapply(.SD, function(x) sdata[x, groupBy[[2]]]), .SDcols = "sample_id"]
      data2plot<- reshape2::melt(data2plot, id.vars=c("interval","sample_id","grp", "grp2"))
      
    } else if(length(groupBy)==3){
      data2plot[, grp := lapply(.SD, function(x) sdata[x, groupBy]), .SDcols = "sample_id"]
      data2plot[, grp2 := lapply(.SD, function(x) sdata[x, groupBy[[2]]]), .SDcols = "sample_id"]
      data2plot[, grp3 := lapply(.SD, function(x) sdata[x, groupBy[[3]]]), .SDcols = "sample_id"]
      data2plot<- reshape2::melt(data2plot, id.vars=c("interval","sample_id","grp", "grp2", "grp3"))
    }
    
    df<- vector("list")
    for(i in unique(data2plot$sample_id)){
        dfi<- data2plot %>% dplyr::filter(sample_id==i)

        if(length(unique(dfi$interval)) != length(unique(data2plot$interval))){
          rown= which(!unique(data2plot$interval) %in% dfi$interval)
          if(length(groupBy)==1){
            add<- data.frame(interval=unique(data2plot$interval)[rown], sample_id=i,grp=unique(dfi$grp),variable="percent", value=0)
            addb<- data.frame(interval=unique(data2plot$interval)[rown], sample_id=i,grp=unique(dfi$grp),variable="freq", value=0)  
          } else if(length(groupBy)==2){
            add<- data.frame(interval=unique(data2plot$interval)[rown], sample_id=i,grp=unique(dfi$grp), grp2=unique(dfi$grp2), variable="percent", value=0)
            addb<- data.frame(interval=unique(data2plot$interval)[rown], sample_id=i,grp=unique(dfi$grp), grp2=unique(dfi$grp2), variable="freq", value=0)
          } else if(length(groupBy)==3){
            add<- data.frame(interval=unique(data2plot$interval)[rown], sample_id=i,grp=unique(dfi$grp), grp2=unique(dfi$grp2), grp3=unique(dfi$grp3),variable="percent", value=0)
            addb<- data.frame(interval=unique(data2plot$interval)[rown], sample_id=i,grp=unique(dfi$grp), grp2=unique(dfi$grp2), grp3=unique(dfi$grp3),variable="freq", value=0)
          }
          
          
          df[[i]]<- rbind(dfi, add,addb)
        } else {
          df[[i]]<- dfi
        }
      }
    
    data2plot<- plyr::ldply(df, data.frame, .id = NULL)
    data2plot<- setDT(data2plot)
    
    if(length(groupBy)==1){
      stat.test1 <- data2plot %>%
        dplyr::filter(variable == "percent") %>%
        dplyr::group_by(interval) %>%
        rstatix::wilcox_test(value  ~ grp) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position()
      
      stat.test2 <- data2plot %>%
        dplyr::filter(variable== "freq") %>%
        dplyr::group_by(interval) %>%
        rstatix::wilcox_test(value  ~ grp) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position()
      
    } else if(length(groupBy)==2){
      stat.test1 <- data2plot %>%
        dplyr::filter(variable == "percent") %>%
        dplyr::group_by(interval, grp2) %>%
        rstatix::wilcox_test(value  ~ grp) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position()
      
      stat.test2 <- data2plot %>%
        dplyr::filter(variable== "freq") %>%
        dplyr::group_by(interval, grp2) %>%
        rstatix::wilcox_test(value  ~ grp) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position()
      
    } else if(length(groupBy)==3){
      stat.test1 <- data2plot %>%
        dplyr::filter(variable == "percent") %>%
        dplyr::group_by(interval, grp2, grp3) %>%
        rstatix::wilcox_test(value  ~ grp) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position()
      
      stat.test2 <- data2plot %>%
        dplyr::filter(variable== "freq") %>%
        dplyr::group_by(interval, grp2, grp3) %>%
        rstatix::wilcox_test(value  ~ grp) %>%
        rstatix::adjust_pvalue() %>%
        rstatix::add_significance() %>%
        rstatix::add_y_position()
    }
    

    p1 <- ggplot2::ggplot(data = data2plot[data2plot$variable == "percent",],
                         ggplot2::aes(x = factor(interval, levels=plotBreaks),
                                      y = .data[["value"]]), alpha=.7) +
      ggplot2::geom_boxplot(ggplot2::aes(fill=.data[["grp"]]),outlier.shape = NA) +
      ggplot2::geom_point(ggplot2::aes(fill=.data[["grp"]]),shape = 21, position=ggplot2::position_dodge(width=.8)) +
      {if(length(groupBy)==2)list(ggplot2::facet_grid(~grp2))} +
      {if(length(groupBy)==3)list(ggplot2::facet_grid(grp2~grp3))} +
      ggplot2::labs(subtitle = "Cumulative frequency")+
      ggplot2::xlab("")+ggplot2::ylab("")+
      ggplot2::scale_color_manual(values = label_colors[[groupBy[[1]]]]) +
      ggplot2::scale_fill_manual(values = label_colors[[groupBy[[1]]]]) +
      theme_RepSeq() +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size=8),
                     axis.text.y = ggplot2::element_text(size=8))+
        ggpubr::stat_pvalue_manual(stat.test1, label = "p.adj.signif",
                            tip.length = 0,x="interval")

    p2 <- ggplot2::ggplot(data = data2plot[data2plot$variable == "freq",],
                          ggplot2::aes(x = factor(interval, levels=plotBreaks),
                                       y = .data[["value"]]), alpha=.7) +
      ggplot2::geom_boxplot(ggplot2::aes(fill=.data[["grp"]]),outlier.shape = NA) +
      ggplot2::geom_point(ggplot2::aes(fill=.data[["grp"]]),shape = 21, position=ggplot2::position_dodge(width=.8)) +
      {if(length(groupBy)==2)list(ggplot2::facet_grid(~grp2))} +
      {if(length(groupBy)==3)list(ggplot2::facet_grid(grp2~grp3))} +
      ggplot2::labs(subtitle = "Distribution")+
      ggplot2::xlab("")+ggplot2::ylab("")+
      ggplot2::scale_color_manual(values = label_colors[[groupBy[[1]]]]) +
      ggplot2::scale_fill_manual(values = label_colors[[groupBy[[1]]]]) +
      theme_RepSeq() +
      ggplot2::theme(legend.position = "right",
                     legend.direction = "vertical",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size=8),
                     axis.text.y = ggplot2::element_text(size=8),
                     legend.background = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size=8),
                     legend.justification = "center" )+
      ggpubr::stat_pvalue_manual(stat.test2, label = "p.adj.signif",
                                 tip.length = 0,x="interval")

    legend<-lemon::g_legend(p2)
    g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p2 + ggplot2::theme(legend.position="none"),
                                                        p1 + ggplot2::theme(legend.position="none"),
                                                        nrow=1),
                                 legend, nrow=1,widths=c(8, 1.3))

     } else{
   data2plot<- reshape2::melt(data2plot, id.vars=c("interval","sample_id"))
   colorBreaks <-  c("1"="#FFD92F","]1, 10]"="#A6D854","]10, 100]"= "#E78AC3" , "]100, 1000]"="#8DA0CB" ,"]1000, 10000]"="#FC8D62","]10000, Inf]"="#66C2A5")
   plotBreaks <- c("1","]1, 10]","]10, 100]", "]100, 1000]" ,"]1000, 10000]","]10000, Inf]")


    p1 <- ggplot2::ggplot(data = data2plot[data2plot$variable == "percent",], ggplot2::aes(x = sample_id, y =value, fill=factor(interval, levels=rev(plotBreaks))) ,  alpha=.7) +
      ggplot2::geom_bar(stat = "identity",  alpha=.8) +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_fill_manual(values=colorBreaks)+
      theme_RepSeq()+
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size=8),
                      axis.text.y = ggplot2::element_text(size=8),
                      legend.position = "none")+
      ggplot2::labs(subtitle = "Cumulative frequency")+
      ggplot2::xlab("")+ggplot2::ylab("")

    p2 <- ggplot2::ggplot(data = data2plot[data2plot$variable == "freq",], ggplot2::aes(x = sample_id, y =value, fill=factor(interval, levels=rev(plotBreaks))) ,  alpha=.7) +
      ggplot2::geom_bar(stat = "identity",  alpha=.8) +
      ggplot2::ylim(0, 1) +
      ggplot2::scale_fill_manual(values=colorBreaks)+
      theme_RepSeq()+
      ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,size=8),
                      axis.text.y = ggplot2::element_text(size=8),
                      legend.position = "right",
                      legend.direction = "vertical",
                      legend.background = ggplot2::element_blank(),
                      legend.text = ggplot2::element_text(size=8),
                      legend.justification = "center")+
      ggplot2::labs(subtitle = "Distribution")+
      ggplot2::xlab("")+ggplot2::ylab("")

    legend<-lemon::g_legend(p2)
    g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p2 + ggplot2::theme(legend.position="none"),
                                   p1 + ggplot2::theme(legend.position="none"),
                                   nrow=1),
                       legend, nrow=1,widths=c(8, 1.3))
     }

}


#' @title Visualization of basic statistics
#'
#' @description This function plots the statistics in the mData slot, calculated for each sample during the building of the RepSeqExperiment object.
#'
#' These statistics include:
#'
#' - nSequences: the total number of sequences in a sample
#'
#' - CDR3nt: the number of unique nucleotide CDR3s
#'
#' - CDR3aa: the number of unique amino acid CDR3s
#'
#' - V: the total number of V genes expressed in each sample
#'
#' - J: the total number of J genes
#'
#' - VJ: the total number of V-J gene combinations
#'
#' - clone: the number of unique clones
#'
#' - clonotype: the number of unique clonotypes
#'
#' They can be compared between groups of samples or simply plotted for each sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param groupBy a character indicating one or multiple groups to be compared. A Wilcoxon test is thus performed and adjusted p-values using the Holm method are shown. Colors are attributed to the different groups within the first column, and a facet is applied on the second column. If not specified, no statistical tests will be performed, and calculated values for each sample_id will be represented. 
#' @param label_colors a list of colors for each variable in groupBy See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param stat a character specifying the statistic to plot. Should one of the statistics in the metaData slot.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotStatistics(x = RepSeqData, groupBy = "sex", stat = "V")
#' 
#' plotStatistics(x = RepSeqData, groupBy = c("cell_subset", "sex"), stat = "V")
#'
plotStatistics <- function(x, stat = c("nSequences", "clone", "clonotype","V", "J","VJ", "CDR3aa", "CDR3nt"),
                           groupBy = NULL, label_colors = NULL){
  group <- NULL
  groupB <- NULL
  groupC <- NULL
  sdata<-mData(x)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  if (is.null(label_colors)) {
    label_colors= plotColors(x)
  }

  if (!is.null(groupBy)){
    
    if(length(groupBy)==1){
      sdata_m<- sdata %>% dplyr::select(sample_id, paste(groupBy),paste(stat))  %>% dplyr::rename(stats=paste(stat))
      colnames(sdata_m)[2]<- "group"
    } else if(length(groupBy)==2){
      sdata_m<- sdata %>% dplyr::select(sample_id, paste(groupBy),paste(stat))  %>% dplyr::rename(stats=paste(stat))
      colnames(sdata_m)[c(2,3)]<- c("group", "groupB")
    } else if(length(groupBy)==3){
      sdata_m<- sdata %>% dplyr::select(sample_id, paste(groupBy),paste(stat))  %>% dplyr::rename(stats=paste(stat))
      colnames(sdata_m)[c(2, 3, 4)]<- c("group", "groupB", "groupC")
    }
    
    if(length(unique(sdata_m[,3]))>1){
      
      if(length(groupBy)==1){
        stat.test <- sdata_m %>%
          rstatix::wilcox_test( stats ~ group) %>%
          rstatix::adjust_pvalue() %>%
          rstatix::add_significance() %>%
          rstatix::add_xy_position()
      } else if(length(groupBy)==2){
        stat.test <- sdata_m %>%
          dplyr::group_by(groupB) %>%
          rstatix::wilcox_test(stats ~ group) %>%
          rstatix::adjust_pvalue() %>%
          rstatix::add_significance() %>%
          rstatix::add_xy_position()
      } else if(length(groupBy)==3){
        stat.test <- sdata_m %>%
          dplyr::group_by(groupB, groupC)  %>%
          rstatix::wilcox_test( stats ~ group) %>%
          rstatix::adjust_pvalue() %>%
          rstatix::add_significance() %>%
          rstatix::add_xy_position()
      }
      
      p <- ggplot2::ggplot(sdata_m,ggplot2::aes(x = group, y = stats, fill = group)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_point(shape=21,position=ggplot2::position_dodge(width=.8) )+
        {if(length(groupBy)==2)list(ggplot2::facet_grid(~groupB))} +
        {if(length(groupBy)==3)list(ggplot2::facet_grid(groupB~groupC))} +
        ggplot2::xlab("")+
        ggplot2::ylab(paste(stat))+
        ggplot2::scale_color_manual(values=label_colors[[groupBy[[1]]]])+
        ggplot2::scale_fill_manual(values=label_colors[[groupBy[[1]]]])+
        theme_RepSeq()+
        ggplot2::theme(legend.position = "none")+
        ggprism::add_pvalue(stat.test,  label = "p.adj.signif",bracket.nudge.y = .01,tip.length = 0, step.increase = 0.1)
      
    } else {
      p <- ggplot2::ggplot(sdata_m,ggplot2::aes(x = group, y = stats, fill = group)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_point(shape=21,position=ggplot2::position_dodge(width=.8) )+
        ggplot2::xlab("")+
        ggplot2::ylab(paste(stat))+
        ggplot2::scale_color_manual(values=label_colors[[groupBy]])+
        ggplot2::scale_fill_manual(values=label_colors[[groupBy]])+
        theme_RepSeq()+
        ggplot2::theme( legend.position = "none")

    }
  }  else {
    sdata_m<- sdata %>% dplyr::select(sample_id, paste(stat))  %>% dplyr::rename(stats=paste(stat))

    p<- ggplot2::ggplot(sdata_m,ggplot2::aes(x = sample_id, y = stats, fill = sample_id)) +
      ggplot2::geom_bar(stat="identity", size=.5, color="black") +
      ggplot2::xlab("")+
      ggplot2::ylab(paste(stat))+
      ggplot2::scale_fill_manual(values=label_colors[["sample_id"]])+
      theme_RepSeq()+
      ggplot2::theme( legend.position = "none",
                      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

  }
  return(p)
}



#' @title Visualization of the perturbation scores
#'
#' @description plots the calculated perturbation scores using the \code{\link{perturbationScore}} function onto a heatmap.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param ctrl.names a vector of characters indicating the sample_ids to be used as controls.
#' @param distance a character specifying the distance method to be used in the calculation of the perturbation scores. Should be one of the following: "manhattan", "euclidean", "canberra", "minkowski" or "maximum".
#' @param p an integer indicating the power of the Minkowski distance. Default is 2.
#' @param order a character specifying the column name in mData to be used as reference to order the sample_ids in the heatmap, as no hierarchical clustering is performed in this analysis.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' plotPerturbationScore(x = RepSeqData,
#'                       ctrl.names = c("tripod-30-813" ,"tripod-31-846"  ,"tripod-35-970"),
#'                       distance = "manhattan",
#'                       order="cell_subset")
#'                       
plotPerturbationScore <- function(x, ctrl.names=NULL,
                                  distance = c("manhattan", "euclidean", "canberra",
                                               "minkowski", "maximum"),
                                    p=2, order="cell_subset",
                                    label_colors=NULL) {
  if (missing(x))
    stop("x is missing, an object of class RepSeqExperiment is expected.")
  if (!is.RepSeqExperiment(x))
    stop("an object of class RepSeqExperiment is expected.")
  snames <- rownames(mData(x))
  if (missing(ctrl.names))
    stop("ctrl.names is missing. A vector of characters containing the names of control samples is expected.")

  sdata <- mData(x)
  groups <- sdata[, unlist(lapply(sdata, is.factor)), drop = FALSE]

  if (is.null(label_colors)) {
    label_colors = plotColors(x)
  }

  per<- perturbationScore(x, ctrl.names, distance =distance)
  sdata<- sdata[order(sdata[[order]]),]
  per <- per[,sdata[["sample_id"]]]
  sNames <- rownames(sdata)

  p <- ComplexHeatmap::pheatmap(t(per),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     annotation_colors = label_colors,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     show_colnames = TRUE, annotation_row = groups, name = " ",
                     show_rownames = FALSE, clustering_method = "ward.D", 
                     silent = FALSE, fontsize = 4.5,angle_col="90")
  return(p)
}


#' @title Visualization of differential expression in a volcano plot
#'
#' @description This function plots differentially expressed repertoire levels calculated using the \code{\link{diffExpGroup}} in a volcano plot.
#'
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "clone","clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param group a vector of character indicating the column name in the mData slot, as well as the two groups to be compared.
#' @param top an integer indicating the top n significant labels to be shown on the plot. Default is 10.
#' @param FC.TH an integer indicating the log2FoldChange threshold. Default is 2.
#' @param PV.TH an integer indicating the adjusted pvalue threshold. Default is 0.05.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' plotVolcano(x = RepSeqData,
#'             level = "V",
#'             group = c("cell_subset", "amTreg", "nTreg"),
#'             top = 10,
#'             FC.TH = 1,
#'             PV.TH = 0.05)
#'
plotVolcano <- function(x,
                        level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"),
                        group = c("cell_subset", "amTreg", "nTreg"),
                        FC.TH=2, PV.TH=0.05,  top=10 ){

  if (missing(x))
    stop("x is missing, an object of class RepSeqExperiment is expected.")
  if (!is.RepSeqExperiment(x))
    stop("an object of class RepSeqExperiment is expected.")

  dds <- .toDESeq2(x, colGrp = group[1], level = level)
  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq2::DESeq(dds, fitType = 'local')
  res <- DESeq2::results(dds, contrast = group)
  res <- as.data.frame(res[order(res$padj),])

  degTab <- res[!is.na(res$padj), ]
  setDT(degTab, keep.rownames = TRUE)
  degTab[, group := "Other"]
  degTab[padj < PV.TH & log2FoldChange >= FC.TH, group := "Over-expression"]
  degTab[padj < PV.TH & log2FoldChange <= FC.TH, group := "Down-expression"]
  degTab[padj < PV.TH & abs(log2FoldChange) < FC.TH, group := "Other"]
  degTab[padj > PV.TH & abs(log2FoldChange) >= FC.TH, group := "Other"]
  degTab[, BHpvalue := -log10(padj)]
  degTab[group == "Over-expression" | group == "Down-expression", labels := rn]
  degTab[is.na(degTab$labels)]$labels <- ""

  x_legend <- list(title = paste0("log2(",paste(group[2]),"/",paste(group[3]), ")"))
  y_legend <- list(title = "-log10(padj)")

  fc_limits <- round(max(abs(degTab$log2FoldChange)), 0)

  p <- ggplot2::ggplot()+
            ggplot2::geom_point(data = degTab, ggplot2::aes(x = log2FoldChange, y = BHpvalue, fill = group, group=rn), shape = 21, size = 2)+
            ggplot2::labs(x = x_legend,
                          y = y_legend)+
            ggrepel::geom_text_repel(data = degTab[seq_len(top)], ggplot2::aes(x = log2FoldChange,
                                                               y = BHpvalue,
                                                               label = labels), size = 3,
                                     size=5, max.overlaps = 9999999999)+
            ggplot2::geom_hline(yintercept = -log10(PV.TH), linetype = "dashed")+
            ggplot2::geom_vline(xintercept = c(-FC.TH, FC.TH), linetype = "dashed")+
            ggplot2::scale_x_continuous(limits = c(-fc_limits, fc_limits),
                                        breaks = seq(-fc_limits, fc_limits))+
            ggplot2::scale_fill_manual(values = c("Other"="gray",
                                                  "Over-expression"="red",
                                                  "Down-expression"="#6495ED"))+
            theme_RepSeq()

  return(p)
}


#' @title Multidimensional Scaling analysis and Principal Component Analysis
#'
#' @description This function can be used to visualize:
#'
#' - repertoire dissimilarities by performing a multidimensional scaling (MDS) on a pairwise distance matrix calculated between samples on a selected repertoire level, using a specific dissimilarity method.
#'  The proposed methods are detailed in the \code{\link{plotDissimilarityMatrix}} function.
#'
#' - differentially expressed repertoire levels calculated using the \code{\link{diffExpGroup}} function by performing a principal component analysis (PCA).
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "clone","clonotype", "V", "J", "VJ", "CDR3nt" or "CDR3aa".
#' @param method a character specifying the distance method to be computed. Should be specified only if the dim_method parameter is set to "MDS", and should be one of the following: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis."
#' @param colorBy a character specifying the column name in mData to be used to attribute group colors. If not specified, a color is attributed per sample_id.
#' @param label_colors a list of colors for each factor column in metaData. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param dim_method a character indicating the dimensional reduction method to be performed. Should be one of "PCA" or "MDS".
#' @details Details on the proposed dissimilarity indices can be found in the vegan package:
#'
#' https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist
#' 
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotDimReduction(x = RepSeqData,
#'                  level = "V",
#'                  method = "euclidean",
#'                  colorBy = "cell_subset",
#'                  dim_method="MDS")
#'
#' plotDimReduction(x = RepSeqData,
#'                 level = "J",
#'                 colorBy = "cell_subset",
#'                 dim_method="PCA")
#'
plotDimReduction <- function(x, level = c("clone","clonotype", "V", "J", "VJ", "CDR3nt","CDR3aa"),
                    method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
                               "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
                               "binomial", "chao", "cao", "mahalanobis"), colorBy = NULL, label_colors=NULL, dim_method=c("MDS", "PCA")) {

  if (missing(x)) stop("x is missing.")
  if(missing(dim_method)) stop("dim_method is missing.")
  if (missing(method) && dim_method=="MDS") stop("method is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  if (is.null(label_colors)) {
    label_colors = plotColors(x)
  }

  sdata <- mData(x)
  sNames <- rownames(sdata)

  if (is.null(colorBy)) {
    colorBy <- "sample_id"
  } else {
    if (length(grep(colorBy, colnames(sdata))) == 0) colorBy <- "sample_id" else colorBy <- colorBy
  }

  levelChoice <- match.arg(level)
  if(dim_method=="MDS"){

    variable <- NULL
    methodChoice <- match.arg(method)
    cols <- c("sample_id", levelChoice, "count")
    tmp <- data.table::copy(assay(x))[, ..cols]
    groups <- sdata[,unlist(lapply(sdata, is.factor)), drop = FALSE]
    dat <- data.table::dcast(data = tmp, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)
    simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag=TRUE, upper=TRUE), .SDcols = sNames]
    fit <- simmat %>% stats::cmdscale( k = 2) %>%  tidyr::as_tibble()
    colnames(fit) <- c("D1" ,"D2")

    fit <- fit %>%
      dplyr::mutate(groups=x@metaData[[colorBy]])

    mnmx<- list()
    for(i in unique(fit$groups)){
      outi <- car::dataEllipse(fit$D1[fit$groups==i],fit$D2[fit$groups==i], levels=c(0.95,0.95), draw = FALSE)
      rng <- do.call(rbind, lapply( outi, function(mtx) apply(mtx,2,range)))
      mnmx[[i]] <- apply(rng, 2, range)
    }
    mnmx <- plyr::ldply(mnmx, data.frame)

    ref= mnmx %>% 
      select(x,y) %>%
      abs(.) %>% max()

    p<-ggpubr::ggscatter(fit, x = "D1", y = "D2",
                      color = "black",
                      fill = "groups",
                      palette = unique(label_colors[[colorBy]]),
                      shape = 21,
                      conf.int = TRUE,
                      #size = 1.5,
                      conf.int.level=0.95,
                      ellipse = TRUE, 
                      repel = TRUE,
                      ggtheme = theme_RepSeq())+
        ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
        ggplot2::geom_vline(xintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
        ggplot2::scale_x_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))+
        ggplot2::scale_y_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))

  } else if(dim_method == "PCA"){

    dds <- .toDESeq2(x, colGrp = colorBy, level = levelChoice)
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::DESeq(dds, fitType = 'local')
    rsd <- DESeq2::rlog(dds)
    datapca <- DESeq2::plotPCA(rsd, intgroup = colorBy, returnData = TRUE)
    percentVar <- round(100 * attr(datapca, "percentVar"))
    
    mnmx<- list()
    for(i in unique(datapca$group)){
      outi <- car::dataEllipse(datapca$PC1[datapca$group==i],datapca$PC2[datapca$group==i], levels=c(0.95,0.95), draw = FALSE)
      rng <- do.call(rbind, lapply( outi, function(mtx) apply(mtx,2,range)))
      mnmx[[i]] <- apply(rng, 2, range)
    }
    mnmx <- plyr::ldply(mnmx, data.frame)
    
    ref= mnmx %>% 
      select(x,y) %>%
      abs(.) %>% max()
    

    p <- ggpubr::ggscatter(datapca, x = "PC1", y = "PC2",
                           color = "black",
                           fill = colorBy,
                           palette = unique(label_colors[[colorBy]]),
                           shape = 21,
                           conf.int = TRUE,
                           size = 1.5,
                           ellipse = TRUE,
                           repel = TRUE,
                           ggtheme = theme_RepSeq())+
      ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
      ggplot2::geom_vline(xintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
      theme_RepSeq()+
      ggplot2::scale_x_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))+
      ggplot2::scale_y_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))

  }

  return(p)

}
