utils::globalVariables(c("..adj.rr.label..", "Group", "group", "method", "ste", "value", "y","grp","variable","stats","adj.rr.label","p.value.label", "where", "se", "to", "CELL_META", "padj", "log2FoldChange", "BHpvalue", "rn",
                         "p.signif","xmin", "xmax", "shapes", "comb", "sample_id.1","ypos","id"))

#' @title  Visualization of the VJ combination usage
#'
#' @description  This function calculates the count or frequency of each combination by taking into account the weight of the chosen repertoire level: either "aaClone" or "ntClone".
#'
#' It offers two types of visualization of the calculated VJ usage in a given sample, either a circos plot or a heatmap.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param scale a character specifying whether to plot the VJ usage in "count" or "frequency".
#' @param level a character specifying the level of the repertoire to be taken into account when calculate VJ usages. Should be one of "aaClone" or "ntClone".
#' @param prop a numeric indicating the proportions of top VJ combinations to be plotted. It ranges from 0 to 1.
#' @param plot a character indicating the type of visualization in which the results will be represented, either a heatmap or a circos plot.
#' @keywords internal
#' @examples
#' 
#' plotVJusage(x = RepSeqData, sampleName = NULL, scale = "count", 
#'             level="aaClone", prop=0.1, plot="Circos")
#' 
#' plotVJusage(x = RepSeqData, sampleName = NULL, scale = "frequency", 
#'             level="ntClone", prop=1, plot="Heatmap")
#'
#' @export
#'
plotVJusage <- function(x, sampleName = NULL, scale = c("count", "frequency"),
                        level = c("aaClone","ntClone"),
                        prop=1,
                        plot=c("Circos","Heatmap")) {
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
  tmp <- cts_b[, .(value = sum(count)), by = .(V, J)]
  
  if(scale == "frequency"){
    tmp <- tmp[,value := prop.table(value)]
  }
  
  data2plot_m<- tmp %>%
                dplyr::filter(value!=0) %>%
                dplyr::arrange(desc(value) ) %>%
                dplyr::slice(seq_len(floor(prop*nrow(.)))) 
 
  if (plot=="Circos"){
  p<- circlize::chordDiagram(data2plot_m, annotationTrack =  "grid", 
                             annotationTrackHeight = 0.03,
                             preAllocateTracks = list(track.height=0.2), 
                             big.gap = 20)
  p1<- circlize::circos.track(track.index=1, panel.fun=function(x,y){
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], circlize::CELL_META$sector.index,
                facing="clockwise", niceFacing = TRUE, adj=c(0,0.5), cex=0.45)
  }, bg.border=NA)
  
 } else {
   
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    data2plot_m<- data2plot_m %>%
                  data.table::dcast(J~V) %>%
                  replace(is.na(.),0) %>%
                  tibble::column_to_rownames("J") 

    values<- c(unname(colSums(data2plot_m)), unname(rowSums(data2plot_m)))
    colors<-colorRampPalette(c( "#01665E", "#F6E8C3", "#B35806" ))(length(values))[as.numeric( as.factor(values))]
    names(values)<- colors
    
    column_ha = ComplexHeatmap::HeatmapAnnotation(V=ComplexHeatmap::anno_barplot(unname(colSums(data2plot_m)),
                                                                                 axis_param=list(gp=grid::gpar(fontsize = 7)),
                                                                                 border = FALSE, 
                                                    gp = grid::gpar(color = "black", fill= names(values)[ unname(values) %in% unname(colSums(data2plot_m))]) ),
                                                  annotation_name_gp=grid::gpar(fontsize = 7), annotation_name_side = "left",  annotation_name_rot = 0)
    row_ha = ComplexHeatmap::rowAnnotation( J = ComplexHeatmap::anno_barplot(unname(rowSums(data2plot_m)),
                                                                             axis_param=list(gp=grid::gpar(fontsize = 7)),
                                                                             border = FALSE,
                                             gp = grid::gpar( color = "black",  fill= names(values)[ unname(values) %in% unname(rowSums(data2plot_m))])),
                                            annotation_name_gp=grid::gpar(fontsize = 7))
    
    p1<-ComplexHeatmap::pheatmap(as.matrix(data2plot_m),col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                          cluster_rows = FALSE, cluster_cols = FALSE, name = " ",
                          row_names_side = "left",top_annotation = column_ha, 
                          right_annotation = row_ha, heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 7), direction = "horizontal"),
                         angle_col="90", fontsize =4.5, scale="column")
    return(ComplexHeatmap::draw(p1, heatmap_legend_side = "bottom" ))
    
    # return(p1)
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
#' @keywords internal
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
  CDR3length=pep=sample_id=aaCDR3=N=frequency=percent <- NULL
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
    data2plot <- cts[sample_id == index,][, CDR3length := nchar(aaCDR3)][,.(.N), by = .(V, CDR3length)]
    data2plot$V<- ifelse(data2plot$V %in% legend_v$V,  data2plot$V , "Other")
     p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = N, fill = V), color="black")

  }else if (scl == "frequency") {
    data2plot <- cts[sample_id == index,][, CDR3length := nchar(aaCDR3)][,.(.N), by = .(V, CDR3length)][,frequency := prop.table(N)]
    data2plot$V<- ifelse(data2plot$V %in% legend_v$V,  data2plot$V , "Other")
     p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = frequency, fill = V), color="black")

  }
  p1 <- p + ggplot2::geom_bar(stat = "identity", position="stack", color="black", linewidth=0.1) +
    ggplot2::xlab("CDR3 length (aa)")+ggplot2::ylab(paste(scl))+
    ggplot2::scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(data2plot$V))))+
    theme_RepSeq()+
    ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size=8), axis.title = ggplot2::element_text(size = 9))

  return(p1)
}

#' @title Visualization of the Hill or Renyi profiles
#'
#' @description This function plots the Renyi Entropy or Hill Diversity at any repertoire level for all the samples in the dataset.
#' The alpha parameters can be personalized, thus allowing to focus on certain indices such as the Shannon index for alpha=1 or the Simpson index for alpha=2.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param alpha a numerical vector specifying the alpha values to compute. If not specified, the following values are estimated: c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf).
#' @param Hill a boolean if TRUE the Hill's index is computed.
#' @param level a character specifying the level of the repertoire to be taken into account when calculating VJ usages. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param grouped a boolean indicating whether or not the mean and se of samples belonging to the same experimental group specified in the ColorBy parameter should be calculated. Grouping will be performed on the group chosen in the colorBy parameter. Default is FALSE.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' 
#' @details Hill diversity uses a parameter q to adjust how species are count. 
#' If q=0, it just counts the number of species (species richness). If q=1, 
#' it gives equal weight to all species, balancing rare and common ones. If q=2, 
#' it focuses more on the most common species.
#' Rényi diversity is a generalization of the Shannon index. It also uses a parameter 
#' alpha that works similarly. When alpha is small, it gives more weight to rare 
#' species. When alpha is large, it focuses more on the dominant species. Alpha =1 
#' is an approximation of the Shannon index; alpha = 2 corresponds to the Simpson
#' index and alpha=Inf corresponds to the Berger-Parker index.
#' The main difference is that Hill diversity is designed to give “effective
#' species numbers” i.e. how many equally abundant species would give the same
#' diversity, while Rényi diversity is based on entropy and doesn't directly
#' translate to species counts. Both are useful for exploring how diversity
#' changes depending on whether you care more about rare or common species.
#' @export
#' @examples
#'
#' data(RepSeqData)
#' 
#' plotGenDiversity(x = RepSeqData, 
#'                alpha = c(0,  1, 2, 8, 16, 32, 64), 
#'                level = "V", 
#'                colorBy = "sex")
#'                
#' plotGenDiversity(x = RepSeqData, 
#'                alpha = c(0,  1, 2, 8, 16, 32, 64), 
#'                level = "V", 
#'                colorBy = "sample_id")
#' 
#'
#' plotGenDiversity(x = RepSeqData, 
#'                 level = "J",
#'                  colorBy = "cell_subset", 
#'                  grouped=TRUE)
#' 
#' plotGenDiversity(x = RepSeqData,
#'               level = "J",
#'                colorBy = "sex", 
#'                facetBy= "cell_subset", 
#'                grouped=TRUE, 
#'                Hill=TRUE)

plotGenDiversity <- function(x, alpha = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf),
                              level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                              Hill = FALSE,
                              colorBy=NULL, 
                              facetBy=NULL,
                              grouped=FALSE,
                              label_colors=NULL) {
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (length(alpha) < 2) stop("At least 2 alpha values are needed.")
  if(is.null(colorBy)) stop("need to specify a group column from mData to assign colors.")
  if(colorBy=='sample_id' & grouped==TRUE) stop("grouped=TRUE is not compatible with colorBy='sample_id'")
  
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }
    
  sdata <- mData(x)
  sNames <- rownames(sdata)
  levelChoice <- match.arg(level)

  gNames <- names(sdata[, c(colorBy, facetBy), drop=FALSE])
  
  tmp <- generalizedDiversity(x, alpha = alpha, level=levelChoice, Hill=Hill )
  data2plot <- data.table::melt(data = tmp, id.vars = "variable", measure.vars = sNames, variable.name = "sample_id")
  data2plot <- data.table::dcast(data2plot,sample_id~ variable )
  data2plot[, (gNames) := lapply(gNames, function(x) sdata[, x] )][, "sample_id"]
  if(colorBy=='sample_id'){
    data2plot <- data.table::melt(data = data2plot, id.vars = c(gNames))
  } else{
    data2plot <- data.table::melt(data = data2plot, id.vars = c(gNames,'sample_id'))
    
  }
  data2plot<- data2plot %>% dplyr::rename(alpha=variable)

  
  if (grouped) {
    se<- function(x) sqrt(var(x)/length(x))
    
    data2plot <-  data2plot %>% dplyr::group_by(!!!rlang::syms(gNames), alpha) %>% dplyr::mutate(ste=se(value))
    data2plot <-  data2plot %>% dplyr::group_by(!!!rlang::syms(gNames), alpha, ste) %>% dplyr::summarize(mean=mean(value)) 

  
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = alpha, y = mean)) +
        ggplot2::geom_line(ggplot2::aes(group = .data[[colorBy]],color=.data[[colorBy]]), linewidth = .8)  +
        ggplot2::geom_point(ggplot2::aes( color=.data[[colorBy]]),  size=1.2) +
        ggplot2::geom_ribbon( ggplot2::aes(ymin=mean-ste, ymax=mean+ste,
                                           fill = .data[[colorBy]],
                                           group = .data[[colorBy]]), alpha = 0.3,colour=NA)+
        ggplot2::xlab("alpha")  +
        ggplot2::ylab("Renyi's Entropy")+
        {if(Hill==TRUE )list(ggplot2::ylab("Hill's Diversity"))} +
        {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]]))} +
         {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
        ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
        ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
        theme_RepSeq()
  
  } else {
    
    if(colorBy=="sample_id"){
      
      p<- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = alpha, y = value, color=.data[[colorBy]])) +
          ggplot2::geom_line(ggplot2::aes( group=.data[[colorBy]], color=.data[[colorBy]]), linewidth = .8)  +
          ggplot2::geom_point(  size=1.2) +
        ggplot2::xlab("alpha")+
        ggplot2::ylab("Renyi's Entropy")+
        {if(Hill==TRUE )list(ggplot2::ylab("Hill's Diversity"))} +
        {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
        {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                        cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
        ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
        theme_RepSeq()+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))+
        {if(length(label_colors[[colorBy]])>20)list(ggplot2::theme( legend.position = "none"))} 

    } else {

    p<- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = alpha, y = value,  
                                                       group =.data[["sample_id"]],
                                                       color=.data[[colorBy]])) +
        ggplot2::geom_line( linewidth = .8) +
        ggplot2::geom_point( size=1.5)+
        ggplot2::xlab("alpha")+
        ggplot2::ylab("Renyi's Entropy")+
        {if(Hill==TRUE )list(ggplot2::ylab("Hill's Diversity"))} +
        {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
        {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                        cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
        ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
        theme_RepSeq()+
        ggplot2::theme(  axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))+
        ggplot2::theme( legend.position = "right")
    }
  }
  return(p)
  }

#' @title Visualization of the V or J gene usage
#'
#' @description This function plots a heatmap of V and J combination usage within the sample of interest. 
#' Frequencies, represented by the color scale, are scaled column-wise. Barplots at the top and right side of the heatmap show the usage of each gene across the row or the column for the V and J genes, respectively. 
#'
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the gene usages. Should be one of "aaClone" or "ntClone".
#' @export
#' @examples
#'
#' data(RepSeqData)
#' snames <- rownames(mData(RepSeqData))
#'
#' plotIndGeneUsage(x = RepSeqData, level = "aaClone", sampleName = snames[1])
#'

# plotIndGeneUsage <- function(x,  sampleName = NULL, level = c("V", "J"), scale = c("count", "frequency")) {
#   frequency <- NULL
#   if (missing(x)) stop("x is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
# 
#   if (is.null(sampleName)) {
#     sName <- rownames(mData(x))[1]
#     cat("Plot the first sample in the dataset:", sName)
#   } else sName <- sampleName
# 
#   levelChoice <- match.arg(level)
#   cts <- data.table::copy(assay(x))
#   data2plot <- cts[sample_id == sName, lapply(.SD, sum), by = levelChoice, .SDcols = "count"]
# 
#   if (scale == "count") {
#     p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x=levelChoice, y = "count", fill = levelChoice))
#   }
#   if (scale == "frequency") {
#     data2plot[, frequency := prop.table(count)]
#     p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x=levelChoice, y = "frequency", fill = levelChoice))
#   }
#   cols<-RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category =="qual" & RColorBrewer::brewer.pal.info$colorblind == "TRUE", ]
#   p1 <- p + ggplot2::geom_bar(stat = "identity", show.legend=FALSE) +
#     ggplot2::scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nrow(unique(data2plot[, 1]))))+
#     theme_RepSeq()+
#     ggplot2::xlab("") +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,size=8),
#                    axis.text.y = ggplot2::element_text(size=8))
# 
# 
#   return(p1)
# }

plotIndGeneUsage <- function(x,  
                             sampleName = NULL ,
                             level = c("aaClone","ntClone")) {
  frequency <- NULL
  levelChoice<- match.arg(level)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")

  if (is.null(sampleName)) {
    sName <- rownames(mData(x))[1]
    cat("Plot the first sample in the dataset:", sName)
  } else sName <- sampleName

  
  cts <- data.table::copy(assay(x))
  cts_b <- cts[sample_id == sName, ][, .(count = sum(count)), by=c("sample_id", levelChoice, "V","J")]
  tmp <- data.table::dcast(cts_b, J~V, value.var = "count", fun.aggregate = sum)
  data2plot <- data.frame(tmp, row.names = 1)
  data2plot_m <- prop.table(data2plot)
  
  
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  
    values<- c(unname(colSums(data2plot_m)), unname(rowSums(data2plot_m)))
    colors<-colorRampPalette(c( "#01665E", "#F6E8C3", "#B35806" ))(length(values))[as.numeric( as.factor(values))]
    names(values)<- colors
    
    column_ha = ComplexHeatmap::HeatmapAnnotation(V=ComplexHeatmap::anno_barplot(unname(colSums(data2plot_m)),
                                                                                 axis_param=list(gp=grid::gpar(fontsize = 7)),
                                                                                 border = FALSE,
                                                                                 gp = grid::gpar(color = "black", fill= names(values)[ unname(values) %in% unname(colSums(data2plot_m))]) ),
                                                  annotation_name_gp=grid::gpar(fontsize = 7), annotation_name_side = "left",  annotation_name_rot = 0)
    row_ha = ComplexHeatmap::rowAnnotation( J = ComplexHeatmap::anno_barplot(unname(rowSums(data2plot_m)),
                                                                             axis_param=list(gp=grid::gpar(fontsize = 7)),
                                                                             border = FALSE,
                                                                             gp = grid::gpar( color = "black",  fill= names(values)[ unname(values) %in% unname(rowSums(data2plot_m))])),
                                            annotation_name_gp=grid::gpar(fontsize = 7))

    p1<-ComplexHeatmap::pheatmap(as.matrix(data2plot_m),col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                 cluster_rows = FALSE, cluster_cols = FALSE, name = " ",
                                 row_names_side = "left",top_annotation = column_ha,
                                 right_annotation = row_ha, heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 5), direction = "horizontal"),
                                 angle_col="90", fontsize =4.5, scale="column")

    return(ComplexHeatmap::draw(p1, heatmap_legend_side = "bottom"))
  }
}

#' @title Compare V or J gene distributions
#'
#' @description This function compares the V or J gene usages between given groups.
#'
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken 
#' into account when calculating the gene usages. Should be one of "aaClone" or 
#' "ntClone".
#' @param scale a character specifying whether to plot the gene usage in "count" or "frequency".
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param show_stats whether to statistically compare groups
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotGeneUsage(x = RepSeqData, 
#'               level = "J", 
#'               scale = "count", 
#'               colorBy = "cell_subset", 
#'               show_stats=TRUE )
#' 
#' plotGeneUsage(x = RepSeqData,
#'              level = "V", 
#'              scale = "frequency",
#'              colorBy = "cell_subset", 
#'              facetBy="sex")
#'

plotGeneUsage <- function(x, level = c("V", "J"), 
                          scale = c("count", "frequency"), 
                          colorBy=NULL,
                          facetBy=NULL, label_colors = NULL,
                          show_stats=FALSE) {
  frequency <- NULL
  sdata<-mData(x)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(colorBy)) stop("a group column from the metadata is expected.")
  if (colorBy=="sample_id") stop("colorBy='sample_id' is not compatible with this function.")
  
   if (is.null(label_colors)) {
     label_colors= oData(x)$label_colors 
   }
  
  scaleChoice <- match.arg(scale)
  levelChoice <- match.arg(level)

  cts <- data.table::copy(assay(x))
  data2plot <- cts[, lapply(.SD, sum), by = c(levelChoice, "sample_id"), .SDcols = "count"][, frequency := count / sum(count), by = sample_id]

  data_exp <- data.table::CJ(sample_id = unique(data2plot$sample_id),
                             levelChoice = unique(data2plot[[levelChoice]]))
  setnames(data_exp, "levelChoice", levelChoice)
  
  data2plot <- dplyr::left_join(data_exp,data2plot, by=c("sample_id",paste(levelChoice)))
  data2plot[is.na(data2plot)] <- 0
  
  data2plot<- merge(data2plot, sdata[,c('sample_id',colorBy,facetBy), drop=FALSE], by="sample_id")
  
  data2plot.summary <-  data2plot %>% 
                    dplyr::group_by(!!!rlang::syms(c(colorBy,facetBy,levelChoice))) %>%
                    dplyr::summarise(
                      sd = sd(get(scaleChoice), na.rm = TRUE),
                      mean = mean(get(scaleChoice)),
                      n = dplyr::n(),
                      se = sd / sqrt(n))
  
  if(show_stats==TRUE){
    stat.test<- .safe_kruskal_pairwise(
      df = data2plot,
      group_var = levelChoice,
      facet_var = facetBy,
      color_var = colorBy,
      value_var = scaleChoice
    )
  }
  
  p <- ggplot2::ggplot(data2plot.summary, 
                      ggplot2::aes(x = !!rlang::sym(levelChoice), y = mean )) +
       ggplot2::geom_bar(ggplot2::aes(fill=.data[[colorBy]]), stat="identity", position = ggplot2::position_dodge()) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-se, ymax=mean+se, group=.data[[colorBy]]), width=.2,position=ggplot2::position_dodge(.9))+
    {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
    {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                    cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
    ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
    theme_RepSeq()+
    ggplot2::xlab("") +
    ggplot2::ylab(paste(scaleChoice)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,size=8),
                   axis.text.y = ggplot2::element_text(size=8),
                   legend.position = "top")+
    {if (show_stats==TRUE)  ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif",
                                  size=4, hide.ns = TRUE)}
      
  return(p)
}


#' @title Visualization of the clonal distribution per interval in a single sample
#'
#' @description This function plots two histograms of the clonal distribution per a set of fixed intervals in a given sample.
#'
#' The upper plot calculates the proportion of each interval in the whole repertoire, whereas the lower plot shows the cumulative frequency of the sequences within each interval.
#'
#' This allows a global view of the repertoire fraction contributing the most to the repertoire. For instance, top sequences belonging to the highest interval often constitute a low fraction of the whole repertoire but contribute more significantly in terms of cumulative frequency in view of their high occurrence.
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the clonal distribution. Should be one of "aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param interval_scale whether intervals should be determined in count or frequency
#' @export
#' @examples
#'
#' data(RepSeqData)
#' 
#' snames <- rownames(mData(RepSeqData))
#'
#' plotIndIntervals(x = RepSeqData, level="aaCDR3", sampleName = snames[1],  interval_scale="count")
#'
#' plotIndIntervals(x = RepSeqData, level="aaClone", sampleName = NULL,  interval_scale="frequency")
#'
plotIndIntervals <- function(x, sampleName = NULL, 
                                  level = c("aaClone","ntClone", "ntCDR3","aaCDR3"),
                                 interval_scale=c("count", "frequency" )){
  interval=cumfreq=distribution=ct <- NULL
  levelChoice<- match.arg(level)
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(sampleName)) {
    sName <- rownames(mData(x))[1]
    cat("Plot the first sample in the dataset:", sName)
  } else sName <- sampleName
  
  data2plot <- data.table::copy(assay(x))
  data2plot<- data2plot[sample_id == sName, lapply(.SD, sum), by = levelChoice, .SDcols = "count"]
  
  data2plot <- if (interval_scale == "count") {
    intervals = c("1", "]1, 10[", "[10, 100[", "[100, 1000[", "[1000, 10000[", "[10000, Inf]")
    colorBreaks = c("1" = "#1F77B4B2", "]1, 10[" = "#FF7F0EB2", "[10, 100[" = "#2CA02CB2", 
                    "[100, 1000[" = "#D62728B2", "[1000, 10000[" = "#9467BDB2", "[10000, Inf]" = "#8C564BB2")
    f = function(x) {
      if (x == 1) "1"
      else if (x > 1 & x < 10) "]1, 10["
      else if (x < 100) "[10, 100["
      else if (x < 1000) "[100, 1000["
      else if (x < 10000) "[1000, 10000["
      else "[10000, Inf]"
    }
    data2plot[, `:=`(interval, unlist(lapply(get(interval_scale), f)))]
    
  } else {
    intervals = c("[0, 0.00001[", "[0.00001, 0.0001[", "[0.0001, 0.001[",
                  "[0.001, 0.01[","[0.01, 1]")
    colorBreaks = c("[0, 0.00001[" = "#1F77B4B2", "[0.00001, 0.0001[" = "#FF7F0EB2",
                    "[0.0001, 0.001[" = "#2CA02CB2", "[0.001, 0.01[" = "#D62728B2", 
                    "[0.01, 1]" = "#9467BDB2")
    f = function(x) {
      if (x < 0.00001) "[0, 0.00001["
      else if (x < 0.0001) "[0.00001, 0.0001["
      else if (x < 0.001) "[0.0001, 0.001["
      else if (x < 0.01) "[0.001, 0.01["
      else "[0.01, 1]"
    }
    data2plot[,frequency := prop.table(count)][, `:=`(interval, unlist(lapply(get(interval_scale), f)))]
  }

  tmp1<- data2plot[, .(ct = .N), by = .(interval)][, distribution := ct / sum(ct)][,ct := NULL]
  tmp2<- data2plot[, .(ct = sum(count)), by = .(interval)][, cumfreq := ct / sum(ct)][, ct := NULL]
  data2plot <- merge(tmp1, tmp2, by = "interval")

  plotBreaks<- intervals
  data2plot<- data.table::melt(data2plot, id.vars="interval")
  data2plot <- data2plot %>% 
                dplyr::group_by(variable) %>%
                dplyr::mutate(ypos = cumsum(value)- 0.5*value) %>%
                dplyr::mutate(variable=ifelse(variable=="cumfreq", "Cumulative frequency", "Distribution")) %>%
                dplyr::mutate(variable=factor(variable, levels=c("Distribution", "Cumulative frequency")))

  
  p <- ggplot2::ggplot(data = data2plot, 
                  ggplot2::aes(x = "", y =value, fill=factor(interval, levels=rev(as.character(plotBreaks))))) +
      ggplot2::geom_bar(stat="identity", width=1) +
      ggplot2::coord_radial("y", start=0, expand=FALSE)+
      ggrepel::geom_label_repel(data = data2plot[ data2plot$value!=0,],
                                ggplot2::aes(label = round(value,3), y = ypos), size=4,
                                show.legend = FALSE) +
      ggplot2::scale_fill_manual(
        values = colorBreaks,
        name = paste(interval_scale, "intervals") 
      )+
      theme_RepSeq()+
      ggplot2::theme(axis.text.theta = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
      ggplot2::facet_wrap(~variable, nrow=1)+
      ggplot2::theme(legend.position = "right")+
      ggplot2::xlab("")+
      ggplot2::ylab("")+
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))

  return(p)
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
#' @keywords internal
#' @examples
#'
#' snames <- rownames(mData(RepSeqData))
#'
#' plotSpectratypingV(x = RepSeqData, sampleName = snames[1], scale = "count", prop=0.05)
#'
#' plotSpectratypingV(x = RepSeqData, sampleName = snames[1], scale = "frequency")
#'
#'
plotSpectratypingV <- function(x, sampleName = NULL, 
                               scale = c("count", "frequency"),
                               prop=0) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  frequency=N=CDR3length=pep=sample_id=aaCDR3=percent <- NULL
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
    legend_v <- legend_v %>%
      dplyr::arrange(desc(N) ) %>%
      dplyr::slice(seq_len(floor(prop*length(unique(cts$V)))))
  }

  if (scl == "count"){
    tmp <- cts[sample_id == index,][, CDR3length := nchar(aaCDR3)]
    data2plot <- tmp[, .(.N), by = .(V, CDR3length)]
    data2plot<- data2plot[data2plot$V %in% legend_v$V,]
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = N))
  }
  if (scl == "frequency"){
    data2plot <- cts[sample_id == index,][, CDR3length := nchar(aaCDR3)][,.(.N), by = .(V, CDR3length)][,percent := prop.table(N), by = V]
    data2plot<- data2plot[data2plot$V %in% legend_v$V,]

    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = percent))
  }

  p <- p + ggplot2::geom_bar(stat = "identity", fill="lightgray", color="black") +
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

#' @title Visualization of repertoire dissimilarities in a heatmap
#'
#' @description This function assesses pairwise repertoire dissimilarities using a specific dissimilarity method.
#'
#' It calculates a list of dissimilarity indices, each taking into account different parameters. The proposed methods include:
#' 
#'  The Jaccard similarity: a measure of similarity between sample sets defined as the size of the intersection divided by the size of the union of the sample sets.
#'
#'  The Morisita-Horn similarity: a measure of similarity that tends to be over-sensitive to abundant species.
#'
#' The function also performs a hierarchical clustering on the calculated distance scores in case the results are represented on a heatmap.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the indices are computed. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param method a character specifying the distance method to be computed. Should be one of the following: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis."
#' @param clustering a character specifying the clustering method to be used in case a heatmap is plotted. If not, the parameter can be set to NULL.
#' @param binary a boolean indicating whether or not to transform the data into a presence/absence data. Default is FALSE
#' @param annotation_groups a vector indicating at least one column name in mData. Colors are thus attributed to the different groups within this column in the MDS. The chosen column must be of class factor.
#' @param label_colors a list of colors for each factor column in metaData. See \code{\link{plotColors}}. If NULL, default colors are used.
#'
#' @details Details on the calculated indices as well as the clustering methods can be found in the vegan package: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotDissHeatmap(x = RepSeqData, level = "V", method = "morisita", annotation_groups="sex",
#'                  clustering="ward.D")
#'
#' plotDissHeatmap(x = RepSeqData, level = "aaCDR3", method = "jaccard", 
#'                   annotation_groups=c("cell_subset", "sex"),  clustering="ward.D")

plotDissHeatmap <- function(x, level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                                    method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
                                               "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
                                               "binomial", "chao", "cao", "mahalanobis"),
                                    clustering=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty" ,
                                                "median", "centroid" ),
                                    binary = FALSE,
                                    annotation_groups=NULL,
                                    label_colors=NULL) {
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(annotation_groups)) stop("at least one group column for color is expected.")
  
  variable <- NULL
  levelChoice <- match.arg(level)
  methodChoice <- match.arg(method)
  clust<- match.arg(clustering)
  cols <- c("sample_id", levelChoice, "count")
  tmp <- data.table::copy(assay(x))[, ..cols]
  sdata <- mData(x)
  sNames <- rownames(sdata)

  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }

  dat <- data.table::dcast(data = tmp, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)
  simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag = TRUE, upper = TRUE, binary = binary), .SDcols=sNames]

  groups <- sdata[, annotation_groups, drop = FALSE]
    
  hide_sample_legend <- "sample_id" %in% annotation_groups
  
  # Create annotation object with conditional legend visibility
  ha <- HeatmapAnnotation(
    df = groups,
    col = label_colors,
    show_legend = if (hide_sample_legend & length(sNames)>20) {
      # Hide legend for "sample_id", show others
      stats::setNames(rep(TRUE, length(annotation_groups)), annotation_groups) %>% 
        replace(names(.) == "sample_id", FALSE)
    } else {
      TRUE  # Show all legends if "sample_id" is not in colorBy
    }
  )
  
  p <- ComplexHeatmap::pheatmap(as.matrix(simmat),
                     cluster_rows = TRUE, cluster_cols = TRUE, name = " ",  col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     treeheight_row = 0L, clustering_distance_rows = simmat, clustering_distance_cols = simmat,
                     top_annotation =ha, show_colnames=FALSE, annotation_colors = label_colors,
                     show_rownames=FALSE, clustering_method = clust, silent = FALSE, fontsize =4)

   return(p)
}

#' @title Visualization of repertoire dissimilarities in a multidimensional scaling (MDS) plot
#'
#' @description This function assesses pairwise repertoire dissimilarities using a specific dissimilarity method.
#'
#' It calculates a list of dissimilarity indices, each taking into account different parameters. The proposed methods include:
#' 
#'  The Jaccard similarity: a measure of similarity between sample sets defined as the size of the intersection divided by the size of the union of the sample sets.
#'
#'  The Morisita-Horn similarity: a measure of similarity that tends to be over-sensitive to abundant species.
#'
#' The function performs multidimensional scaling (MDS) on the calculated dissimilarity scores. This enables visualization of repertoire relationships in a reduced-dimensional space, highlighting similarities and differences among samples.
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the indices are computed. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param method a character specifying the distance method to be computed. Should be one of the following: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis."
#' @param binary a boolean indicating whether or not to transform the data into a presence/absence data. Default is FALSE
#' @param colorBy a vector indicating one column name in mData. Colors are thus attributed to the different groups within this column in the MDS. The chosen column must be of class factor.
#' @param shapeBy a vector indicating one column name in mData. Shapes are thus attributed to the different groups within this column
#' @param label_colors a list of colors for each factor column in metaData. See \code{\link{plotColors}}. If NULL, default colors are used.
#'
#' @details Details on the calculated indices as well as the clustering methods can be found in the vegan package: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotDissMDS(x = RepSeqData, level = "V", method = "euclidean", colorBy="sex")
#'
#' plotDissMDS(x = RepSeqData, level = "aaCDR3", method = "jaccard", 
#'                   colorBy="cell_subset", shapeBy="sex")

plotDissMDS <- function(x, level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                              method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
                                         "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
                                         "binomial", "chao", "cao", "mahalanobis"),
                              binary = FALSE,
                              colorBy=NULL,
                              shapeBy= NULL,
                              label_colors=NULL) {
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(colorBy)) stop("at least one group column for color is expected.")
  if( colorBy=="sample_id") stop("colorBy='sample_id' is not compatible with this function.")

  
  variable <- NULL
  levelChoice <- match.arg(level)
  methodChoice <- match.arg(method)
  cols <- c("sample_id", levelChoice, "count")
  tmp <- data.table::copy(assay(x))[, ..cols]
  sdata <- mData(x)
  sNames <- rownames(sdata)
  
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }
  
  dat <- data.table::dcast(data = tmp, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)
  simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag = TRUE, upper = TRUE, binary = binary), .SDcols=sNames]

  fit <- simmat %>% stats::cmdscale( k = 2) %>%  tibble::as_tibble(.name_repair="unique" )
  colnames(fit) <- c("D1" ,"D2")
  
  if(is.null(shapeBy)){
    fit <- fit %>%
            dplyr::mutate(!!colorBy := x@metaData[[colorBy]])
    
    p<-ggpubr::ggscatter(fit, x = "D1", y = "D2",
                         color = colorBy,
                         palette = label_colors[[colorBy]] ,
                         shape = 21,
                         conf.int = TRUE,
                         conf.int.level=0.95,
                         ellipse = TRUE, 
                         repel = TRUE,
                         ggtheme = theme_RepSeq())+
      ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
      ggplot2::geom_vline(xintercept = 0, color = "gray", size = 0.1, linetype = "dashed")
    
  } else if (!is.null(shapeBy)){
    fit <- fit %>%
      dplyr::mutate(!!colorBy := x@metaData[[colorBy]],
                    !!shapeBy := x@metaData[[shapeBy]])
    
    p<-ggpubr::ggscatter(fit, x = "D1", y = "D2",
                         color = colorBy,
                         palette = label_colors[[colorBy]] ,
                         shape = shapeBy,
                         conf.int = TRUE,
                         conf.int.level=0.95,
                         ellipse = TRUE, 
                         repel = TRUE,
                         ggtheme = theme_RepSeq())+
      ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
      ggplot2::geom_vline(xintercept = 0, color = "gray", size = 0.1, linetype = "dashed")
  } 

  return(p)
}


#' @title Visualization of the clonal distribution as a function of the rank
#'
#' @description This function plots the clonal distribution, at any repertoire level, as a function of the occurrence rank within each sample.
#'
#' For instance, the most frequent sequence is attributed a rank of 1, and its relative abundance is plotted on the y-axis.
#'  
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account to plot the clonal distribution. Should be one of "aaClone", "ntClone", "ntCDR3" or "aaCDR3".
#' @param scale a character specifying whether to plot the clonal abundance in "count" or "frequency" or "log" scale.
#' @param colorBy a character indicating one column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' @param ranks an integer specifying the number of top ranks to be plotted. Default is NULL, which plots all ranks.
#' @param grouped a boolean indicating whether or not the mean and se of samples belonging to the same experimental group specified in the ColorBy parameter should be computed. Grouping is performed on the chosen groupe in colorBy. Default is FALSE.
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @export
#' @examples
#' 
#' data(RepSeqData)
#' 
#' plotRankDistrib(x=RepSeqData, level = "aaCDR3",
#'                  scale="log",
#'                  ranks=10,
#'                  colorBy="sample_id" )
#'
#' plotRankDistrib(x=RepSeqData, level = "aaCDR3",
#'                  scale="log",
#'                  ranks=10,
#'                  colorBy="sex" )
#'                  
#' plotRankDistrib(x=RepSeqData, level = "aaCDR3",
#'                  scale="count",
#'                  ranks=10,
#'                  colorBy="cell_subset",
#'                  grouped=TRUE )
#'                                  
plotRankDistrib <- function(x, level = c("aaClone","ntClone", "ntCDR3","aaCDR3"),
                            scale=c("count","frequency", "log"),
                            ranks=NULL,
                            grouped=FALSE,
                            colorBy=NULL, 
                            facetBy=NULL,
                            label_colors=NULL){
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if(is.null(colorBy)) stop("need to specify a column from mData. Can be the sample_id column")
  if(colorBy=="sample_id" & grouped==TRUE) stop("grouped cannot be used when colorBy is set to sample_id")

  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }
  
  scl <- match.arg(scale)
  levelChoice <- match.arg(level)
    
  sName <- rownames(mData(x))
  sdata <- mData(x)

  counts <- data.table::copy(assay(x))
  counts <- counts[, .(count = sum(count)), by=c("sample_id", levelChoice)][, rank := lapply(.SD, frankv, ties.method = "first", order = -1L), by = "sample_id", .SDcols = "count"]
  
  counts<- merge(counts, sdata[,c(colorBy,facetBy, "sample_id","nSequences")],by="sample_id")

  if (grouped) {
    
    se<- function(x) sqrt(var(x)/length(x))
    
     if (scl == "frequency"){
      counts <- counts %>% dplyr::group_by(sample_id) %>% 
                dplyr::mutate(count=count/sum(count))
      counts <- data.table::setDT(counts)

      counts <-  counts %>% dplyr::group_by(!!!rlang::syms(c(colorBy,facetBy)), rank) %>% dplyr::mutate(ste=se(count))
      counts <-  counts %>% dplyr::group_by(!!!rlang::syms(c(colorBy,facetBy)), rank, ste) %>% dplyr::summarize(mean=mean(count))
      
    } else {
      counts <-  counts %>% dplyr::group_by(!!!rlang::syms(c(colorBy,facetBy)), rank) %>% dplyr::mutate(ste=se(count))
      counts <-  counts %>% dplyr::group_by(!!!rlang::syms(c(colorBy,facetBy)), rank, ste) %>% dplyr::summarize(mean=mean(count))
    }
    
    if(!is.null(ranks)) counts<- counts %>% dplyr::filter(rank<=ranks)
    
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = rank, y = mean, fill = .data[[colorBy]])) +
      ggplot2::geom_point(data=counts[counts$rank==1,],ggplot2::aes(), shape=21) +
      ggplot2::geom_ribbon(
                ggplot2::aes(ymin = mean-ste, ymax = mean+ste),
                alpha = 0.3, colour = NA)+
      ggplot2::geom_line(ggplot2::aes(colour = .data[[colorBy]]))+
      {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
      {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
      ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
      ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
      ggplot2::scale_x_log10()+
      {if(scl == "log")list(ggplot2::scale_y_log10())} +
      ggplot2::ylab(paste0("mean ", scl))+
      theme_RepSeq()+
      ggplot2::theme(legend.position = "right", plot.subtitle = ggplot2::element_text(hjust=0.90, vjust=-10))
    
  } else {
    
    if (scl == "frequency"){
      counts <- counts %>% dplyr::group_by(sample_id) %>% dplyr::mutate(count = count/sum(count))
      counts <- data.table::setDT(counts)
    } else {
      counts <- counts 
    }
    
    if(!is.null(ranks)) counts<- counts %>% dplyr::filter(rank<=ranks)
    
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = rank, y = count, color=.data[[colorBy]])) +
      ggplot2::geom_point(data=counts[counts$rank==1,], shape=21) +
      ggplot2::geom_line(ggplot2::aes(group = .data[['sample_id']], color=.data[[colorBy]]), linewidth = .8) +
      ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
     {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
      {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
      ggplot2::scale_x_log10()+
      {if(scl == "log")list(ggplot2::scale_y_log10())} +
     theme_RepSeq()+
     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))+
     {if(length(label_colors[[colorBy]])>20) ggplot2::theme( legend.position = "none") else ggplot2::theme( legend.position = "right")}+
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
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param sampleNames a vector of character indicating the sample_ids of the repertoires to be analyzed. If not specified, the first three samples in the dataset are analyzed.
#' @export
#' @examples
#' 
#'
#' data(RepSeqData)
#'
#' snames <- rownames(mData(RepSeqData))[1:4]
#'
#' plotVenn(x = RepSeqData, level = "V", sampleNames = snames)
#'
#' plotVenn(x = RepSeqData, level = "aaCDR3", sampleNames = NULL)
#'
plotVenn <- function(x, level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                     sampleNames = NULL) {
  ..sampleNames <- NULL
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  levelChoice <- match.arg(level)
  sdata <- mData(x)
  label_colors= oData(x)$label_colors 

  if (is.null(sampleNames)) sampleNames <- rownames(sdata)[seq_len(3)]
  if (length(sampleNames) == 1) stop("At least 2 sampleNames must be provided.")

  counts <- data.table::copy(assay(x))[sample_id %in% sampleNames]
  list<-list()
  for ( i in unique(counts$sample_id)){
    counts_i<- counts %>% dplyr::filter(sample_id==i)
    list[[i]]<- unique(counts_i[[levelChoice]])
  }

 
 plot <- ggVennDiagram::ggVennDiagram(list, label_alpha = 0, edge_lty="solid", 
                                  show_intersect=FALSE,label="both",label_percent_digit = 1,
                                  set_size=3,edge_size=.5 ,set_color="black",
                                  label_geom="label", label_color="black",
                                  label_size=2.5)+
         ggplot2::scale_color_manual(values=rep("gray",length(sampleNames)))+
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
#' @param level a character specifying the level of the repertoire to be analyzes. Should be one of "aaClone", "ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param sampleNames a vector of character indicating the two sample_ids of the repertoires to be analyzed.
#' @param scale a character specifying whether to plot in "count" , "frequency". or "log'
#' @param shiny default is FALSE. whether to plot the shiny compatible version
#' @export
#' @examples
#'
#' data(RepSeqData)
#' plotScatter(x = RepSeqData,
#'             level = "V",
#'             sampleNames = c("tripod-30-813", "tripod-30-815", "tripod-31-846"),
#'             scale = "log")
#'             
#' plotScatter(x = RepSeqData,
#'             level = "aaClone",
#'             sampleNames = c("tripod-30-813", "tripod-30-815"),
#'             scale = "frequency")
#'
#'
plotScatter <- function(x, sampleNames = NULL,
                        level = c("aaClone", "ntClone", "V", "J", "VJ", "ntCDR3", "aaCDR3"),
                        scale = c("frequency", "count", "log"),
                        shiny = FALSE) {
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected.")
  if (length(sampleNames) < 2) stop("At least two sampleNames are required.")
  
  levelChoice <- match.arg(level)
  scaleChoice <- match.arg(scale)
  cols <- c("sample_id", levelChoice, "count")
  
  cts <- data.table::copy(assay(x))
  cts$count <- as.numeric(cts$count)
  counts <- cts[sample_id %in% sampleNames, ..cols]
  
  data2plot <- data.table::dcast(counts, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)
  
  if (scaleChoice == "frequency") 
    data2plot <- cbind(data2plot[, 1], as.data.frame(apply(data2plot[, -1], 2, function(x) x / sum(x))))
 
  sNames <- vapply(sampleNames, function(s) {
    if (!grepl("^`", s)) s <- paste0("`", s)
    if (!grepl("`$", s)) s <- paste0(s, "`")
    # if (grepl('tripod', s)) {
    #   s <- gsub("-", " ", s)
    #   return(stringr::str_wrap(s, width = 6))
    # } else {
    #   return(s)
    # }
    
    # if (nchar(s) > 6) {
    #   return(stringr::str_wrap(s, width = 6))
    # } else {
    #   return(s)
    # }
  }, FUN.VALUE = character(1))
  
  
  data2plot <- data.frame(data2plot, check.names = FALSE)
  
  if (length(sampleNames) == 2) {
    if (shiny) {
      correlation <-  stats::cor.test(data2plot[, 2], data2plot[, 3], method = "pearson")
      p <- ggplot2::ggplot(data2plot, ggplot2::aes_string(x = sNames[1], y = sNames[2])) +
        ggplot2::geom_point(ggplot2::aes_string(label = levelChoice), size = 2, shape = 21, alpha = 0.5) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5, color = "red") +
        ggplot2::annotate("text", label = paste0("R = ", round(correlation$estimate, 2), "\n", "p = ", round(correlation$p.value, 3)),
                          hjust = 0, x = min(data2plot[[sNames[1]]]), y = max(data2plot[[sNames[1]]])) +
        theme_RepSeq()
    } else {
      p <- ggplot2::ggplot(data2plot, ggplot2::aes_string(x = sNames[1], y = sNames[2])) +
        ggplot2::geom_count(size = 2, shape = 21, alpha = 0.5) +
        smplot2::sm_statCorr(linetype = "dashed", linewidth = 0.5, color = "red", label.fontface = "bold") +
        theme_RepSeq()
    }
  } else {
    
    if(scaleChoice == "log"){
      tmp <- GGally::ggpairs(data2plot[, -1], sampleNames,
                           lower = list(continuous = log10_points, alpha = 0.5) ,
                          upper = list(continuous = GGally::wrap("cor", size = 3)),
                            diag = NULL)
    } else {
      tmp <- GGally::ggpairs(data2plot[, -1], sampleNames,
                           lower = list(continuous = GGally::wrap(GGally::ggally_points, alpha = 0.5)),
                            upper = list(continuous = GGally::wrap("cor", size = 3)),
                           diag = NULL)
    }
    
     p <- tmp +
        theme_RepSeq() +
        ggplot2::xlab("") +
        ggplot2::ylab("") 
  
   }
  
  return(p)
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
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param index a character specifying the diversity index to be estimated. Should be one of "shannon","invsimpson","simpson", "bergerparker" or "gini".
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param grouped a character indicating one or multiple groups to be compared. A Wilcoxon test is thus performed and adjusted p-values using the Holm method are shown. Colors are attributed to the different groups within the first column, and a facet is applied on the second column. If not specified, no statistical tests will be performed, and calculated values for each sample_id will be represented. 
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param show_stats whether to statistically compare groups
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' 
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotDiversity(x = RepSeqData, level = "V",  colorBy = "sample_id",
#'               facetBy="cell_subset", index="shannon")
#' 
#' plotDiversity(x = RepSeqData, level = "ntCDR3", colorBy = "cell_subset", facetBy="sex", 
#'               grouped=TRUE, index="simpson")
#'
#' plotDiversity(x = RepSeqData, level = "aaClone", colorBy = "cell_subset", facetBy="sex", 
#'               grouped=TRUE, index="shannon",  show_stats=TRUE)

plotDiversity <- function(x, index=c("shannon","simpson", "invsimpson","bergerparker", "gini"),
                          level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                          grouped = FALSE, 
                          colorBy=NULL,
                          facetBy=NULL,
                          label_colors=NULL,
                          show_stats=FALSE){
  
  unq<-length(mData(x)[, colorBy])
 
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (grouped && is.null(colorBy) || grouped && dplyr::n_distinct(mData(x)[, colorBy])==unq) stop("A valid column for colorBy is required when grouped is TRUE")
  if (is.null(index)) stop("a diversity index is expected.")
  if(is.null(colorBy)) stop("need to specify a group column from mData to assign colors.")
  
  
  levelChoice <- match.arg(level)
  sdata<-mData(x)
  
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }
  
  diversity <- diversityIndices(x, level=levelChoice)
  diversity_m <- diversity %>% dplyr::select(sample_id, paste(index)) %>% 
                   dplyr::rename(method = paste(index))
  diversity_m<- merge(diversity_m, sdata[,c('sample_id',colorBy,facetBy), drop=FALSE], by="sample_id")
  
  
  if(grouped==TRUE){
    
    p <-  ggpubr::ggboxplot(
          data = diversity_m,
          x = colorBy,
          y = "method",
          fill = colorBy,
          outlier.shape = NA,
          add = "jitter",
          shape = 21 ) +
          ggplot2::xlab("") +
          ggplot2::ylab(index) +
          theme_RepSeq() +
          ggplot2::scale_fill_manual(values=label_colors[[colorBy]]) +
          ggplot2::theme(legend.position = "none")
        
    # Add faceting
    if (length(facetBy) == 1) {
      p <- p + ggplot2::facet_grid(as.formula(paste("~", facetBy[1])), scales = "free")
    } else if (length(facetBy) == 2) {
      p <- p + ggplot2::facet_grid(as.formula(paste(facetBy[1], "~", facetBy[2])), scales = "free")
    }
    
    # Add p-value annotation if requested
    if (show_stats == TRUE) {
      
        stat.test <- .safe_kruskal_pairwise(
        df = diversity_m,
        group_var = NULL,
        facet_var = facetBy,
        color_var = colorBy,
        value_var = "method" )
        
      p <- p + ggprism::add_pvalue(stat.test, label = "p.adj.signif",tip.length = 0,
                                   step.increase = 0.1)
    }

  }  else {

  
    p<- ggplot2::ggplot(diversity_m,ggplot2::aes(x = sample_id, y = method, fill =.data[[colorBy]])) +
      ggplot2::geom_bar(stat="identity", linewidth=.5, color="black") +
      ggplot2::xlab("")+
      ggplot2::ylab(paste(index))+
      {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
      {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
      ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
      theme_RepSeq()+
      ggplot2::theme(  axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
      {if(length(label_colors[[colorBy]])>20) ggplot2::theme( legend.position = "none") else ggplot2::theme( legend.position = "right")}
  
  }
  return(p)
}


#' @title Visualization of a repertoire structure through a circular treemap
#'
#' @description This function plots all clones in a single repertoire as a circular treemap, showing the clone distributiion within the repertoire. Each circle represents a unique clone, and the circle size corresponds to the clone count.
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param prop a numeric indicating the proportions of clones to be plotted. It ranges from 0 to 1.

#' @export
#' @examples
#'
#' data(RepSeqData)
#' 
#' sname <- rownames(mData(RepSeqData))[1]
#' 
#' plotIndMap(x=RepSeqData, sampleName=sname, level="aaClone", prop=0.01)

plotIndMap <- function (x, sampleName=NULL,
                               prop = 0.01,
                               level = c("aaClone", "ntClone", "V", "J", "VJ", "ntCDR3", "aaCDR3")){
  if (missing(x)) 
    stop("x is missing. An object of class RepSeqExperiment is expected.")
  if (!is.RepSeqExperiment(x)) 
    stop("an object of class RepSeqExperiment is expected.")
  if (prop > 1) 
    stop("prop should be a proportion between 0 and 1")
  if (is.null(sampleName)) {
    sName <- rownames(mData(x))[1]
    cat("Plot the first sample in the dataset:", sName)
  } else sName <- sampleName
  
  levelChoice <- match.arg(level)
  
  cts <- data.table::copy(assay(x))
  cts_b <- cts[, .(count = sum(count)), by = c("sample_id", 
                                               levelChoice)][, `:=`(ranks, lapply(.SD, frankv, ties.method = "first", 
                                                                                  order = -1L)), by = "sample_id", .SDcols = "count"]
  keep <- cts_b[cts_b[, .I[ranks %in% seq_len(ceiling(.N *  prop))], by = "sample_id"]$V1]
  keep <- keep[sample_id == sName]
  
  packing <- packcircles::circleProgressiveLayout(keep$count, sizetype='area')
  packing$radius <- 0.95*packing$radius
  data <- cbind(keep, packing)
  dat.gg <- packcircles::circleLayoutVertices(packing, npoints=50)
  dat.gg$value <- rep(data$count, each=51)
  
 p<-  ggplot2::ggplot() + 
    ggplot2::geom_polygon(data = dat.gg, ggplot2::aes(x, y, group = id, fill=value), colour = "black", linewidth=0.2, alpha = 0.6) +
    ggplot2::scale_fill_distiller(palette = "Spectral", direction = 1 ) +
    ggplot2::theme_void() + 
    ggplot2::theme(legend.position="none")+ 
    ggplot2::coord_equal()
   
 return(p)                          
}



#' @title Visualization of the rarefaction curves
#'
#' @description This function plots the rarefaction curve for each sample within the dataset.
#'
#' Rarefaction is a measure of species richness. The curves plot the number of clones against the number of sequences in a sample, each being obtained by randomly re-sampling a number of sequences multiple times and representing the mean number of found clones.
#'
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
#' plotRarefaction(x = RepSeqData, colorBy = "sample_id")
#'
#'
plotRarefaction <- function(x, colorBy=NULL, label_colors=NULL){
  if (missing(x)) stop("x is required. An object of class RepSeqExperiment is expected.")
  if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected")
  if(is.null(colorBy)) stop("need to specify a group column from mData to assign colors.")
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }

  cts <- data.table::copy(assay(x))
  sdata <- mData(x)
  raretab <- cts[, round(rarefyDT(count),2), by = sample_id][, (colorBy) := lapply(.SD, function(x) sdata[x, colorBy] ), .SDcols = "sample_id"]

  if(colorBy=="sample_id"){
 
  p <- ggplot2::ggplot(data = raretab, ggplot2::aes(x = x, y = y, color = .data[[colorBy]])) +
    ggplot2::geom_line(ggplot2::aes(group=sample_id)) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      x = "Number of sequences",
      y = "Number of clones") +
    ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
    theme_RepSeq()+
    {if(length(label_colors[[colorBy]])>20) list(ggplot2::theme(legend.position = "none"))}
  
    if(length(label_colors[[colorBy]])<10) 
   p <- p + 
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .1))+
    ggplot2::theme(legend.position = "none")+
    ggrepel::geom_text_repel(data=raretab %>% dplyr::group_by(sample_id) %>% dplyr::slice_max(x),
                             nudge_x = -0.1, direction = "y", hjust = "left", size=3 ,ggplot2::aes(label = sample_id)) 
  
  # if(length(sdata$sample_id)>10) 
  #  p <- p+ ggrepel::geom_text_repel(data=raretab %>% dplyr::group_by(sample_id) %>% dplyr::slice_max(x),
  #               nudge_x = -0.1, direction = "y", hjust = "left", size=3 ,ggplot2::aes(label = sample_id)) 
  } else {
    
    p <- ggplot2::ggplot(data = raretab, ggplot2::aes(x = x, y = y,  color = .data[[colorBy]])) +
      ggplot2::geom_line(ggplot2::aes(group=sample_id)) +
      ggplot2::guides(fill = "none") +
      ggplot2::labs(
        x = "Number of sequences",
        y = "Number of clones") +
      ggplot2::scale_color_manual(values=label_colors[[colorBy]])+
      # ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
      theme_RepSeq()+
      ggplot2::theme(legend.position = "right")+
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .1))
    
  }

  return(p)
}


#' @title Visualization of the clonal distribution per interval
#'
#' @description This function plots two histograms of the clonal distribution per a set of intervals in all the samples within the dataset. 
#'
#' The plot titled "Distribution" calculates the proportion of each interval in the whole repertoire, whereas the one titled "Cumulative frequency" shows the cumulative frequency of the sequences within each interval.
#'
#' This could allow a global view of the repertoire fraction contributing the most to the repertoire. For instance, top sequences belonging to the highest interval often constitute a low fraction of the whole repertoire but contribute more significantly in terms of cumulative frequency in view of their high occurrence.
#'
#' Samples can be statistically compared in each interval using the \code{groupBy} parameter.
#' @param x an object of class  \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire to be taken into account when calculating the clonal distribution. Should be one of "aaClone","ntClone", "ntCDR3" or "aaCDR3".
#' @param grouped a character indicating one or multiple groups to be compared. A Wilcoxon test is thus performed and adjusted p-values using the Holm method are shown. Colors are attributed to the different groups within the first column, and a facet is applied on the second column. If not specified, no statistical tests will be performed, and calculated values for each sample_id will be represented. 
#' @param label_colors a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param interval_scale whether intervals should be determined in count or frequency
#' @param calculation_type a character indicating the type of calculation to plot, either the distribution or the cumulative frequency
#' @param show_stats whether to statistically compare groups using a wilcoxon test
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotIntervals(x = RepSeqData, level="aaCDR3", 
#'                facetBy = 'sex', interval_scale="count",
#'                calculation_type = "distribution")
#'                
#' plotIntervals(x = RepSeqData, level="ntCDR3", colorBy="cell_subset", grouped=TRUE, 
#'               interval_scale="frequency", show_stats=TRUE ,
#'               calculation_type = "cumulative frequency")
#' 

plotIntervals <- function(x,
                          level = c("aaClone","ntClone", "ntCDR3","aaCDR3"),
                          colorBy=NULL,
                          facetBy=NULL,
                          label_colors=NULL, 
                          grouped=FALSE, 
                          show_stats=FALSE, 
                          calculation_type = c("distribution", "cumulative frequency"),
                          interval_scale=c("count", "frequency")){
  
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }
  
  interval=percent=nas=ct <- NULL
  sdata<-mData(x)
  unq<-length(sdata[, colorBy])
  levelChoice <- match.arg(level)
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (!grouped && !is.null(colorBy)) stop("colorBy cannot be used when grouped is set to FALSE. Only facetBy can be use in this case. Groups will be colored based on intervals.")
  if (grouped && is.null(colorBy) || grouped && dplyr::n_distinct(mData(x)[, colorBy])==unq) stop("A valid column for colorBy is required when grouped is TRUE")
  

  data2plot <- data.table::copy(assay(x))
  data2plot<- data2plot[, lapply(.SD, sum), by = c(levelChoice, "sample_id"), .SDcols = "count"]
  
  data2plot <- if (interval_scale == "count") {
    intervals = c("1", "]1, 10[", "[10, 100[", "[100, 1000[", "[1000, 10000[", "[10000, Inf]")
    colorBreaks = c("1" = "#1F77B4B2", "]1, 10[" = "#FF7F0EB2", "[10, 100[" = "#2CA02CB2", 
                    "[100, 1000[" = "#D62728B2", "[1000, 10000[" = "#9467BDB2", "[10000, Inf]" = "#8C564BB2")
    f = function(x) {
      if (x == 1) "1"
      else if (x > 1 & x < 10) "]1, 10["
      else if (x < 100) "[10, 100["
      else if (x < 1000) "[100, 1000["
      else if (x < 10000) "[1000, 10000["
      else "[10000, Inf]"
    }
    data2plot[, `:=`(interval, unlist(lapply(get(interval_scale), f))), by = "sample_id"]
    
  } else {
    intervals = c("[0, 0.00001[", "[0.00001, 0.0001[", "[0.0001, 0.001[",
                  "[0.001, 0.01[","[0.01, 1]")
    colorBreaks = c("[0, 0.00001[" = "#1F77B4B2", "[0.00001, 0.0001[" = "#FF7F0EB2",
                    "[0.0001, 0.001[" = "#2CA02CB2", "[0.001, 0.01[" = "#D62728B2", 
                    "[0.01, 1]" = "#9467BDB2")
    f = function(x) {
      if (x < 0.00001) "[0, 0.00001["
      else if (x < 0.0001) "[0.00001, 0.0001["
      else if (x < 0.001) "[0.0001, 0.001["
      else if (x < 0.01) "[0.001, 0.01["
      else "[0.01, 1]"
    }
    data2plot[,frequency := prop.table(count),by="sample_id"][, `:=`(interval, unlist(lapply(get(interval_scale), f))), by = "sample_id"]
  }
  
  
  if(calculation_type == "distribution"){
    data2plot<- data2plot[, .(ct = .N), by = .(interval, sample_id)][, value := ct / sum(ct), by = sample_id][, ct := NULL]
  } else if (calculation_type == "cumulative frequency"){
    data2plot<- data2plot[, .(ct = sum(count)), by = .(interval, sample_id)][, value := ct / sum(ct), by = sample_id][, ct := NULL]
  } else {
    stop("Invalid type. Choose either 'distribution' or 'cumulative frequency'.")
  }
  
  data2plot<- merge(data2plot, sdata[,c('sample_id',colorBy,facetBy), drop=FALSE], by="sample_id")
  data2plot$interval<- factor(data2plot$interval, levels=intervals[intervals %in% data2plot$interval])
  plotBreaks<- intervals
  
  if(grouped){
    
    if(show_stats==TRUE){
      stat.test <- .safe_kruskal_pairwise(
        df = data2plot,
        group_var = "interval",
        facet_var = facetBy,
        color_var = colorBy,
        value_var = "value"
      )
    }

    p <- ggplot2::ggplot(data = data2plot,
                         ggplot2::aes(x = factor(interval, levels=plotBreaks),
                                      y = value), alpha=.7) +
      ggplot2::geom_boxplot(ggplot2::aes(fill=.data[[colorBy]]),outlier.shape = NA, position=ggplot2::position_dodge(width=.8)) +
      ggplot2::geom_jitter(ggplot2::aes(fill=.data[[colorBy]]),shape=21,position=ggplot2::position_jitterdodge() )+
      {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
      {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])), scales="free"))} +
      ggplot2::xlab(paste(interval_scale,'intervals'))+ggplot2::ylab(calculation_type)+
      ggplot2::scale_fill_manual(values = label_colors[[colorBy]]) +
      theme_RepSeq() +
      ggplot2::theme(
        legend.position = "right",
        plot.subtitle=ggplot2::element_text(size=10),
        axis.text.x = ggplot2::element_text(size=8, angle=45, hjust=1),
        axis.text.y = ggplot2::element_text(size=8))+
      {if(show_stats==TRUE) ggpubr::stat_pvalue_manual(stat.test, label = "p.adj.signif",
                                                       size=4, hide.ns = TRUE)}
    
  } else {
    
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = sample_id, y =value, fill=factor(interval, levels=rev(plotBreaks))) ,  alpha=.7) +
      ggplot2::geom_bar(stat = "identity",  alpha=.8) +
      ggplot2::scale_fill_manual(values=colorBreaks, 
                                 name = paste(interval_scale, 'intervals'))+
      {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
      {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
      theme_RepSeq()+
      ggplot2::theme( axis.text.x = ggplot2::element_text( size=7,angle = 90),
                      axis.text = ggplot2::element_text(size=7),
                      legend.position = "right",
                      strip.text.x = ggplot2::element_text( size = 7),
                      strip.text.y = ggplot2::element_text( size = 7),
                      plot.subtitle=ggplot2::element_text(size=7))+
      ggplot2::xlab("")+ggplot2::ylab(calculation_type)

  }
  return(p)
}



#' @title Visualization of basic or diversity statistics for one sample
#'
#' @description This function plots the statistics in the mData slot and highlights one specific sample amongst all samples.
#'
#' Basic statistics include:
#'
#' - nSequences: the total number of sequences in a sample
#'
#' - ntCDR3: the number of unique nucleotide CDR3s
#'
#' - aaCDR3: the number of unique amino acid CDR3s
#'
#' - V: the total number of V genes expressed in each sample
#'
#' - J: the total number of J genes
#'
#' - VJ: the total number of V-J gene combinations
#'
#' - aaClone: the number of unique aaClones
#'
#' - ntClone: the number of unique ntClones
#' 
#' Diversity statistics include:
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
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3". Only required when plotting diversity indices.
#' @param stat a character specifying whether to plot basic or diversity statistics.
#' @param sampleName a character specifying the sample_id to analyze. Default is NULL, which plots the first sample in the dataset.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotIndStatistics(x = RepSeqData, stat = "metadata", level= "ntCDR3" )
#' 
#' plotIndStatistics(x = RepSeqData, stat = "diversity",  level = "aaClone")
#'

plotIndStatistics <- function(x, sampleName=NULL,
                               stat=c("metadata","diversity"),
                               level = c("aaClone","ntClone", "ntCDR3","aaCDR3")){
  
  if (missing(x)) 
    stop("x is missing. An object of class RepSeqExperiment is expected.")
  if (!is.RepSeqExperiment(x)) 
    stop("an object of class RepSeqExperiment is expected.")
  if (is.null(sampleName)) {
    sName <- rownames(mData(x))[1]
    cat("Plot the first sample in the dataset:", sName)
  } else sName <- sampleName
  
  if(stat=="diversity"){
  if (is.null(level))
    stop("a repertoire level is required.")

  levelChoice <- match.arg(level)
    
  dat <-diversityIndices(x, level=levelChoice)
  dat_s <- dat %>% 
          data.table::melt("sample_id") 

  p<- ggplot2::ggplot(dat_s, ggplot2::aes(x=variable, y=value))+ 
    ggplot2::geom_boxplot( width = 0.35,  outlier.shape=21)+
    ggplot2::geom_point(data=dat_s[dat_s$sample_id == sName,], shape=4, color="red", size=3) +
    ggplot2::facet_wrap(~variable, scales="free", nrow=1)+
    theme_RepSeq()+
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank() )+
    ggplot2::xlab("")
  
  } else {
  
    dat<- mData(x)
    dat_s <- dat %>% 
            dplyr::select(sample_id, nSequences, V,J,VJ, ntCDR3,aaCDR3,aaClone,ntClone) %>%
            setDT() %>%
            data.table::melt(id.vars="sample_id")
    
    p<- ggplot2::ggplot(dat_s, ggplot2::aes(x=variable, y=value))+ 
      ggplot2::geom_boxplot( width = 0.35,  outlier.shape=21)+
      ggplot2::geom_point(data=dat_s[dat_s$sample_id == sName,], shape=4, color="red", size=3) +
      ggplot2::facet_wrap(~variable, scales="free", nrow=2)+
      theme_RepSeq()+
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank() )+
      ggplot2::xlab("")+
      ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=2))
    
  }
  return(p)
}



#' @title Visualization of basic statistics
#'
#' @description This function plots the statistics in the mData slot, calculated for each sample during the building of the RepSeqExperiment object.
#'
#' These statistics include:
#'
#' - nSequences: the total number of sequences in a sample
#'
#' - ntCDR3: the number of unique nucleotide CDR3s
#'
#' - aaCDR3: the number of unique amino acid CDR3s
#'
#' - V: the total number of V genes expressed in each sample
#'
#' - J: the total number of J genes
#'
#' - VJ: the total number of V-J gene combinations
#'
#' - aaClone: the number of unique aaClones
#'
#' - ntClone: the number of unique ntClones
#' 
#' Or any other user-defined numeric statistics in the metadata.
#'
#' They can be compared between groups of samples or simply plotted for each sample.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param grouped a character indicating one or multiple groups to be compared. A Wilcoxon test is thus performed and adjusted p-values using the Holm method are shown. Colors are attributed to the different groups within the first column, and a facet is applied on the second column. If not specified, no statistical tests will be performed, and calculated values for each sample_id will be represented. 
#' @param label_colors a list of colors for each variable in groupBy See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param stat a character specifying any numeric column name in the metadata to plot.
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @param show_stats whether to statistically compare groups
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotStatistics(x = RepSeqData, stat = "V", colorBy = "sample_id")
#' 
#' plotStatistics(x = RepSeqData, colorBy = "cell_subset", 
#'               facetBy="sex", stat = "aaClone", grouped=TRUE)
#'
plotStatistics <- function(x, stat = NULL,
                            colorBy=NULL,
                            facetBy=NULL,
                            label_colors = NULL,
                            grouped=FALSE,  
                            show_stats=FALSE){
  
  unq<-length(mData(x)[, colorBy])
  sdata<-mData(x)
  
  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (grouped && is.null(colorBy) || grouped && dplyr::n_distinct(mData(x)[, colorBy])==unq) stop("A valid column for colorBy is required when grouped is TRUE")
  if (is.null(stat)) stop("a statistic to plot is expected.")
  if(! stat %in% colnames(sdata)) stop("the chosen statistic is not in the metadata slot.")
  if(!is.numeric(sdata[[stat]])) stop("a numeric column is expected.")
  if(is.null(colorBy)) stop("need to specify a group column from mData to assign colors.")
  
  
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }
  
  sdata_m<- sdata[,c('sample_id',colorBy,facetBy, stat)]
  
  if(grouped==TRUE){
    
    if(show_stats==TRUE){
      
      stat.test <- .safe_kruskal_pairwise(
            df = sdata_m,
            group_var = NULL,
            facet_var = facetBy,
            color_var = colorBy,
            value_var = stat )
          }
  
    p<-  ggpubr::ggboxplot(
      data = sdata_m,
      x = colorBy,
      y = stat,
      fill = colorBy,
      outlier.shape = NA,
      add = "jitter",
      shape = 21 ) +
      ggplot2::xlab("") +
      ggplot2::ylab(stat) +
      theme_RepSeq() +
      ggplot2::scale_fill_manual(values=label_colors[[colorBy]]) +
      ggplot2::theme(legend.position = "none")
    # Add faceting
    if (length(facetBy) == 1) {
      p <- p + ggplot2::facet_grid(as.formula(paste("~", facetBy[1])), scales = "free")
    } else if (length(facetBy) == 2) {
      p <- p + ggplot2::facet_grid(as.formula(paste(facetBy[1], "~", facetBy[2])), scales = "free")
    }
    
    # Add p-value annotation if requested
    if (show_stats == TRUE) {
      p <- p + ggprism::add_pvalue(stat.test, label = "p.adj.signif",tip.length = 0,
                                   step.increase = 0.1)
    }

  }  else {
   
    if(colorBy=="sample_id")  sdata_m<- sdata[,c(colorBy,facetBy, stat)]

    p<- ggplot2::ggplot(sdata_m,ggplot2::aes(x = sample_id, y = .data[[stat]], fill = .data[[colorBy]])) +
      ggplot2::geom_bar(stat="identity", linewidth=.5, color="black") +
      ggplot2::xlab("")+
      ggplot2::ylab(paste(stat))+
      {if(length(facetBy)==1)list(ggplot2::facet_grid(~.data[[facetBy[1]]], scales="free"))} +
      {if(length(facetBy)==2)list(ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facetBy[1])),
                                                      cols = ggplot2::vars(!!rlang::sym(facetBy[2])),scales="free"))} +
      ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
      theme_RepSeq()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))+
      {if(length(label_colors[[colorBy]])>20) ggplot2::theme( legend.position = "none") else ggplot2::theme( legend.position = "right")}

  }
  return(p)
}


#' @title Visualization of basic statistics in a scatter plot
#'
#' @description This function plots the values of two statistics in the mData slot along two axes for all the samples in the datatset. 
#' It shows the relationship between them. 
#'
#' These statistics include:
#'
#' - nSequences: the total number of sequences in a sample
#'
#' - ntCDR3: the number of unique nucleotide CDR3s
#'
#' - aaCDR3: the number of unique amino acid CDR3s
#'
#' - V: the total number of V genes expressed in each sample
#'
#' - J: the total number of J genes
#'
#' - VJ: the total number of V-J gene combinations
#'
#' - aaClone: the number of unique aaClones
#'
#' - ntClone: the number of unique ntClones
#' 
#' - Chao1: Estimates undetected species using the information on the rarest species (the numbers of singletons and doubletons) (Chao, 1984).
#'  
#' - Improved Chao1: An extension of Chao1 which uses additional information, namely, the numbers of tripletons and quadrupletons  (Chiu et al., 2014).
#'
#'
#'Or any other user-defined numeric statistics in the metadata.
#'
#' @param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' @param stat1 a character specifying the first numeric column name in the metadata to plot.
#' @param stat2 a character specifying the second numeric column name in the metadata to plot.
#' @param label_colors a list of colors for each variable in groupBy See \code{\link{plotColors}}. If NULL, default colors are used.
#' @param facetBy a vector of character indicating one or two column names in mData to apply a facet on.
#' @param colorBy a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.
#' @export
#' @examples
#'
#' data(RepSeqData)
#'
#' plotStatScatter(x = RepSeqData, stat1 = "nSequences", stat2 = "aaClone", colorBy = "sample_id")
#' 
#' plotStatScatter(x = RepSeqData, stat1 = "ntClone", stat2 = "aaClone", colorBy = "cell_subset")
#' plotStatScatter(x = RepSeqData, stat1 = "chao1", stat2 = "aaClone", colorBy = "cell_subset")

plotStatScatter <- function(x, stat1 = NULL, stat2=NULL,
                            colorBy=NULL,
                            facetBy=NULL,
                            label_colors = NULL){
  
  sdata<-mData(x)

  if (missing(x)) stop("x is missing.")
  if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
  if (is.null(stat1) | is.null(stat2)) stop("two statistics are required")
  if(any(!(c(stat1,stat2) %in% colnames(sdata)))) stop("the chosen statistics are not in the metadata slot.")
  if(any(!is.numeric(sdata[[stat1]]) | !is.numeric(sdata[[stat2]]))) stop("two numeric columns are expected.")
  if(is.null(colorBy)) stop("need to specify a group column from mData to assign colors.")
  
  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
  }

  if(colorBy=="sample_id"){
    sdata_m <-sdata[,c(colorBy,facetBy, stat1, stat2)]
  } else {
    sdata_m <- sdata[,c('sample_id',colorBy,facetBy,stat1, stat2)]
  }
  
  p<- ggplot2::ggplot(sdata_m, ggplot2::aes(x = .data[[stat1]], y =.data[[stat2]], fill = .data[[colorBy]])) +
    ggplot2::geom_point(size = 2, shape=21) +
    ggplot2::xlab(stat1)+
    ggplot2::ylab(stat2)+
    {if(length(facetBy)==1)list(ggplot2::facet_wrap(~.data[[facetBy[1]]], scales="free"))} +
    {if(length(facetBy)==2)list(ggplot2::facet_wrap(.data[[facetBy[1]]]~.data[[facetBy[2]]],scales="free"))} +
    ggplot2::scale_fill_manual(values=label_colors[[colorBy]])+
    ggplot2::geom_abline(linetype="dashed",linewidth=0.5, color="red") + 
    theme_RepSeq()+
    {if(length(label_colors[[colorBy]])>20) ggplot2::theme( legend.position = "none") else ggplot2::theme( legend.position = "right")}
  
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
#' @keywords internal
#' @examples
#'
#' data(RepSeqData)
#' plotPerturbationScore(x = RepSeqData,
#'                       ctrl.names = c("tripod-30-813" ,"tripod-31-846"  ,"tripod-35-970"),
#'                       distance = "manhattan",
#'                       order= "cell_subset")
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
  if(is.null(distance)) stop(" a distance method is expected.")
  if (missing(ctrl.names))
    stop("ctrl.names is missing. A vector of characters containing the names of control samples is expected.")

  sdata <- mData(x)
  groups <- sdata[, unlist(lapply(sdata, is.factor)), drop = FALSE][-1]

  if (is.null(label_colors)) {
    label_colors= oData(x)$label_colors 
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
#' @param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' @param group a vector of character indicating the column name in the mData slot, as well as the two groups to be compared.
#' @param top an integer indicating the top n significant labels to be shown on the volcano plot. Default is 10. 
#' @param FC.TH an integer indicating the log2FoldChange threshold. Default is 2. 
#' @param PV.TH an integer indicating the adjusted pvalue threshold. Default is 0.05. 
#' @export
#' @examples
#' 
#' plotDiffExp(x = RepSeqData,
#'             level = "V",
#'             group = c("cell_subset", "amTreg", "nTreg"),
#'             top = 10,
#'             FC.TH = 1,
#'             PV.TH = 0.05)
#'             
plotDiffExp <- function(x,
                        level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
                        group = c("cell_subset", "amTreg", "nTreg"),
                        FC.TH=2, 
                        PV.TH=0.05, 
                        top=10){

  if (missing(x))
    stop("x is missing, an object of class RepSeqExperiment is expected.")
  if (!is.RepSeqExperiment(x))
    stop("an object of class RepSeqExperiment is expected.")
  if (any(grepl(paste(c("\\+","\\-"), collapse="|"),group[-1]))) 
    stop("Subgroups should not contain (-) or (+).")
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The DESeq2 package is required but not installed. Please install it using BiocManager::install('DESeq2').")
  }
  
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

  fc_limits <- ceiling(max(abs(degTab$log2FoldChange)))

  p <- ggplot2::ggplot()+
        ggplot2::geom_point(data = degTab, 
                            ggplot2::aes(x = log2FoldChange,
                                         y = BHpvalue, 
                                         fill = group, group=rn), shape = 21, size = 2)+
        ggplot2::labs(x = x_legend,
                      y = y_legend)+
        ggrepel::geom_text_repel(data = degTab[degTab$group!="Other",][seq_len(top)], ggplot2::aes(x = log2FoldChange,
                                                           y = BHpvalue,
                                                           label = labels), size = 3,
                                                           max.overlaps = 9999999999)+
        ggplot2::geom_hline(yintercept = -log10(PV.TH), linetype = "dashed")+
        ggplot2::geom_vline(xintercept = c(-FC.TH, FC.TH), linetype = "dashed")+
        ggplot2::scale_x_continuous(limits = c(-fc_limits, fc_limits))+
               
        ggplot2::scale_fill_manual(values = c("Other"="gray",
                                              "Over-expression"="red",
                                              "Down-expression"="#6495ED"))+
        theme_RepSeq()
 
  return(p)
}


#' title Multidimensional Scaling analysis and Principal Component Analysis
#' 
#' description This function can be used to visualize:
#' 
#' - repertoire dissimilarities by performing a multidimensional scaling (MDS) on a pairwise distance matrix calculated between samples on a selected repertoire level, using a specific dissimilarity method.
#' 
#' - differentially expressed repertoire levels calculated using the \code{\link{diffExpGroup}} function by performing a principal component analysis (PCA).
#' 
#' param x an object of class \code{\linkS4class{RepSeqExperiment}}
#' param level a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".
#' param method a character specifying the distance method to be computed. Should be specified only if the dim_method parameter is set to "MDS", and should be one of the following: "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis."
#' param groupBy a character specifying up to three column names in mData to be used to attribute group colors.
#' param label_colors a list of colors for each factor column in metaData. See \code{\link{plotColors}}. If NULL, default colors are used.
#' param dim_method a character indicating the dimensional reduction method to be performed. Should be one of "PCA" or "MDS".
#' details Details on the proposed dissimilarity indices can be found in the vegan package:
#' 
#' https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist
#' 
#' export
#' keywords internal
#' examples
#' 
#' data(RepSeqData)
#' 
#' plotDimReduction(x = RepSeqData,
#'                  level = "V",
#'                  method = "euclidean",
#'                  groupBy = "cell_subset",
#'                  dim_method="MDS")
#' 
#' plotDimReduction(x = RepSeqData,
#'                 level = "J",
#'                 groupBy = "cell_subset",
#'                 dim_method="PCA")

# plotDimReduction <- function(x, level = c("aaClone","ntClone", "V", "J", "VJ", "ntCDR3","aaCDR3"),
#                     method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
#                                "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#                                "binomial", "chao", "cao", "mahalanobis"), groupBy = NULL, label_colors=NULL, dim_method=c("MDS", "PCA")) {
# 
#   if (missing(x)) stop("x is missing.")
#   if(missing(dim_method)) stop("dim_method is missing.")
#   if (missing(method) && dim_method=="MDS") stop("method is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
# 
#   sdata <- mData(x)
#   sNames <- rownames(sdata)
# 
#   if (is.null(groupBy)) stop("a group column from the metadata must be specified.")
#   
#    # if (length(grep(colorBy, colnames(sdata))) == 0) colorBy <- "sample_id" else colorBy <- colorBy
#   
#   levelChoice <- match.arg(level)
#   if(dim_method=="MDS"){
# 
#     variable <- NULL
#     methodChoice <- match.arg(method)
#     cols <- c("sample_id", levelChoice, "count")
#     tmp <- data.table::copy(assay(x))[, ..cols]
#     groups <- sdata[,unlist(lapply(sdata, is.factor)), drop = FALSE]
#     dat <- data.table::dcast(data = tmp, paste(levelChoice, "~sample_id"), value.var = "count", fun.aggregate = sum)
#     simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag=TRUE, upper=TRUE), .SDcols = sNames]
#     fit <- simmat %>% stats::cmdscale( k = 2) %>%  tidyr::as_tibble(.name_repair="unique" )
#     colnames(fit) <- c("D1" ,"D2")
# 
#     if(length(groupBy)==1){
#     fit <- fit %>%
#       dplyr::mutate(groups=x@metaData[[groupBy]])
#     } else if (length(groupBy)==2){
#       fit <- fit %>%
#       dplyr::mutate(groups=paste(x@metaData[,groupBy[1]], x@metaData[,groupBy[2]]))
#       mData(x)=mData(x) %>% mutate(groups=as.factor(fit$groups))
#     } else if (length(groupBy)==3){
#       fit <- fit %>%
#         dplyr::mutate(groups=paste(x@metaData[,groupBy[1]], x@metaData[,groupBy[2]], x@metaData[,groupBy[3]]))
#       mData(x)=mData(x) %>% mutate(groups=as.factor(fit$groups))
#       }
#       
#     if (is.null(label_colors)) {
#       label_colors = plotColors(x, samplenames = FALSE)
#     }
#     
#     mnmx<- list()
#     for(i in unique(fit$groups)){
#       outi <- car::dataEllipse(fit$D1[fit$groups==i],fit$D2[fit$groups==i], levels=c(0.95,0.95), draw = FALSE)
#       rng <- do.call(rbind, lapply( outi, function(mtx) apply(mtx,2,range)))
#       mnmx[[i]] <- apply(rng, 2, range)
#     }
#     mnmx <- plyr::ldply(mnmx, data.frame)
# 
#     ref= mnmx %>% 
#       select(x,y) %>%
#       abs(.) %>% max()
# 
#     p<-ggpubr::ggscatter(fit, x = "D1", y = "D2",
#                       color = "black",
#                       fill = "groups",
#                       palette = if(length(groupBy)==1) label_colors[[groupBy]] else unique(label_colors$groups),
#                       shape = 21,
#                       conf.int = TRUE,
#                       #size = 1.5,
#                       conf.int.level=0.95,
#                       ellipse = TRUE, 
#                       repel = TRUE,
#                       ggtheme = theme_RepSeq())+
#         ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
#         ggplot2::geom_vline(xintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
#         ggplot2::scale_x_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))+
#         ggplot2::scale_y_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))
#   } else if(dim_method == "PCA"){
# 
#     dds <- .toDESeq2(x, colGrp = groupBy, level = levelChoice)
#     dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
#     dds <- DESeq2::DESeq(dds, fitType = 'local')
#     rsd <- DESeq2::rlog(dds)
#     datapca <- DESeq2::plotPCA(rsd, intgroup = groupBy, returnData = TRUE)
#     percentVar <- round(100 * attr(datapca, "percentVar"))
#     
#     mnmx<- list()
#     for(i in unique(datapca$group)){
#       outi <- car::dataEllipse(datapca$PC1[datapca$group==i],datapca$PC2[datapca$group==i], levels=c(0.95,0.95), draw = FALSE)
#       rng <- do.call(rbind, lapply( outi, function(mtx) apply(mtx,2,range)))
#       mnmx[[i]] <- apply(rng, 2, range)
#     }
#     mnmx <- plyr::ldply(mnmx, data.frame)
#     
#     ref= mnmx %>% 
#       select(x,y) %>%
#       abs(.) %>% max()
#     
# 
#     p <- ggpubr::ggscatter(datapca, x = "PC1", y = "PC2",
#                            color = "black",
#                            fill = groupBy,
#                            palette = unique(label_colors[[groupBy]]),
#                            shape = 21,
#                            conf.int = TRUE,
#                            size = 1.5,
#                            ellipse = TRUE,
#                            repel = TRUE,
#                            ggtheme = theme_RepSeq())+
#       ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#       ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#       ggplot2::geom_hline(yintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
#       ggplot2::geom_vline(xintercept = 0, color = "gray", size = 0.1, linetype = "dashed")+
#       theme_RepSeq()+
#       ggplot2::scale_x_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))+
#       ggplot2::scale_y_continuous(limits = c(-ref-(ref/10), ref+(ref/10)))
# 
#   }
# 
#   return(p)
# 
# }



#' @title Calculate wilcoxon
#'
#' @description calculate wilcoxon tests in all plots
#'
#' @details function calculates wilcoxon tests in all plots 
#'
#' @param df dataframe
#' @param group_var x axis
#' @param facet_var facet group(s)
#' @param color_var color group
#' 
#' @return dataframe of statistical test results
#' @export
#' @keywords internal




.safe_kruskal_pairwise <- function(df, group_var, facet_var, color_var, value_var) {
  if (is.null(group_var)) {
    group_var <- 'group'
    df[[group_var]] <- "group"
  }
  
  if (!is.null(facet_var)) {
    df <- df %>% dplyr::group_by(!!rlang::sym(group_var), !!rlang::sym(facet_var))
  } else {
    df <- df %>% dplyr::group_by(!!rlang::sym(group_var))
  }
  

 test <- df %>%
    dplyr::group_modify(~{
      subdf <- .x
      color_levels <- length(unique(subdf[[color_var]]))
      n_unique <- dplyr::n_distinct(subdf[[value_var]])
      # If only one group, return NA
      if (color_levels < 2 || n_unique < 4) {
        return(tibble::tibble(
          test = NA_character_,
          n_groups = color_levels
        ))
      }
      
      if (color_levels == 2) {
        return(tibble::tibble(
          test = "Wilcoxon",
          n_groups = color_levels
        ))
        
      } else {
        return(tibble::tibble(
          test = "Pairwise Wilcoxon",
          n_groups = color_levels
        ))
      }
    })
 
 test=stats::na.omit(test)
 if (unique(test$test) ==  "Wilcoxon") {
   print('Performing Wilcoxon test with Bonferroni correction for 2 groups')
   
   res<- df %>% 
     merge(., test, by = c(group_var, facet_var)) %>%
     {
       if(!is.null(facet_var)) {
         dplyr::group_by(., !!rlang::sym(group_var), !!rlang::sym(facet_var))
       } else {
         dplyr::group_by(., !!rlang::sym(group_var))
       }
     } %>%
     rstatix::wilcox_test(., as.formula(paste(value_var, "~", color_var))) %>%
     rstatix::adjust_pvalue(method = "bonferroni") %>%
     rstatix::add_significance() %>%
     rstatix::add_xy_position()
     
 } else if (unique(test$test) == "Pairwise Wilcoxon"){
   print('Performing Pairwise Wilcoxon test with Bonferroni correction for 3 groups')
   
   res<- df %>% 
     merge(., test, by = c(group_var, facet_var)) %>%
     {
       if(!is.null(facet_var)) {
         dplyr::group_by(., !!rlang::sym(group_var), !!rlang::sym(facet_var))
       } else {
         dplyr::group_by(., !!rlang::sym(group_var))
       }
     } %>%
     rstatix::pairwise_wilcox_test(
     .,
     formula = as.formula(paste(value_var, "~", color_var)),
     p.adjust.method = "bonferroni"  ) %>%
     add_xy_position(x = group_var, dodge = .8)
   
 }
   
 res<- res %>%
    dplyr::ungroup()
}


#' @title Calculate wilcoxon safe
#'
#' @description calculate wilcoxon tests in all plots
#'
#' @details function calculates wilcoxon tests in all plots 
#'
#' @param df dataframe
#' @param group_var x axis
#' @param facet_var facet group(s)
#' @param color_var color group
#' 
#' @return dataframe of statistical test results
#' @export
#' @keywords internal

.safe_wilcox<-function(df, group_var, facet_var, color_var, value_var){
  if (is.null(group_var)) {
    group_var = 'group'
    
    df[[group_var]]= "group"
  }
  # if (is.null(group_var)) {
  #   if (!is.null(facet_var)) {
  #     df <- df %>% dplyr::group_by(!!sym(facet_var))
  #   } else {
  #     df <- df %>% dplyr::group_by(!!sym(color_var))
  #   }
  # } else {
  if (!is.null(facet_var)) {
    df <- df %>% dplyr::group_by(!!rlang::sym(group_var), !!rlang::sym(facet_var))
  } else {
    df <- df %>% dplyr::group_by(!!rlang::sym(group_var))
  }
  # }
  
  
  df %>%
    dplyr::group_modify(~{
      subdf <- .x
      color_levels <- unique(subdf[[color_var]])
      # Check if color_var has exactly 2 levels
      if (length(color_levels) != 2) {
        return(tibble::tibble(
          p.value = NA_real_,
          n1 = sum(subdf[[color_var]] == color_levels[1]),
          n2 = ifelse(length(color_levels) > 1, sum(subdf[[color_var]] == color_levels[2]), NA_integer_)
        ))
      }
      n1 <- sum(subdf[[color_var]] == color_levels[1])
      n2 <- sum(subdf[[color_var]] == color_levels[2])
      # Check if both groups have at least one observation
      if (n1 < 1 | n2 < 1) {
        return(tibble::tibble(p.value = NA_real_, n1 = n1, n2 = n2))
      }
      # Perform Wilcoxon test
      
      res <- rstatix::wilcox_test(subdf, as.formula(paste(value_var, "~", color_var)))
      
      res<- res %>% rstatix::adjust_pvalue(method = "bonferroni") %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position()
      
      
    }) %>%
    dplyr::ungroup()
}
