#' @title AnalyzAIRR - AIRR-compliant R package for the analysis of bulk Ig/TCR repertoire datasets
#' 
#' @description AnalyzAIRR is a AIRR-compliant R package developed to analyze bulk Ig/TCR repertoire datasets.
# 
#' @details 
#' AnalyzAIRR allows a general data exploration to evaluate the homogeneity 
#' within defined groups, identify outliers and filter them out. AnalyzAIRR also proposes
#' a set of diversity measures and statistical metrics applicable at any level of granularity. Thus, single-sample repertoire explorations or in-depth cross-comparisons of AIRR datasets 
#' can be conducted leading to ready-to-publish visualization graphics. 
#' AnalyzAIRR is complemented with a guided workflow to help users in their analytical strategy, and a Shiny web application making it user-friendly for biologists with little or no background in bioinformatics.
#
#' @author 
#' \itemize{
#' \item Vanessa Mhanna
#' \item Gabriel Pires
#' \item Grégoire Bohl
#' \item Karim El Soufi
#' \item Nicolas Tchitchek
#' \item David Klatzmann
#' \item Adrien Six
#' \item Hang-Phuong Pham
#' \item Encarnita Mariotti-Ferrandiz 
#' }
#
# 
#' @seealso 
#' Useful links:
#' \itemize{
#'  \item \href{https://github.com/i3-unit/AnalyzAIRR}{https://github.com/i3-unit/AnalyzAIRR}
#'  \item \href{https://github.com/vanessajmh/Shiny-AnalyzAIRR}{https://github.com/vanessajmh/Shiny-AnalyzAIRR}
#' }
#
#
# 
# Imports

# @importClassesFrom data.table data.table
#' @import data.table
#' @import utils
#' @import graphics
#' @import pbapply
#' @import parallel
#' @import methods
#' @import ggplot2
#' @import stringr
#' @import GGally
#' @import grDevices
#' @import ade4
#' @import rstatix
#' @import ggprism
#' @import ggpubr
#' @import tibble
#' @import tidyverse
#' @import smplot2
#' @import vegan
#' @import naturalsort
#' @import RColorBrewer
#' @import packcircles
#' @import ComplexHeatmap
#' @importFrom testthat test_that
#' @importFrom car dataEllipse
#' @importFrom rlang sym syms
#' @importFrom circlize chordDiagram circos.track circos.text 
#' @importFrom stats median sd as.formula var frequency cmdscale lag lm cor.test na.omit
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr group_by summarize n mutate filter select rename
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom R.utils decompressFile
#' @name AnalyzAIRR
NULL


#"_PACKAGE"
#> [1] "_PACKAGE"
