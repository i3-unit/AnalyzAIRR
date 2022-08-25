##### variable cleaning and environment setup #####
rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)

##### load libraries #####
library("roxygen2")
library("devtools")
library("knitr")
library("rmarkdown")
library("BiocCheck")
library("Rcpp")
library("BiocStyle")

##### set working directory #####
if(R.version$os == "linux-gnu"){
  setwd("/mnt/mukkuri/RepSeq/RS_Analysis/GPI/AnalyzAIRR/")
} else {
  setwd("/Users/vanessamhanna/Documents/PostDoc/AnalyzAIRR/")
}

##### remove previous package folders
unlink("./package", recursive = TRUE)
unlink("./package.zip", recursive = TRUE)
unlink("./package.tar.gz", recursive = TRUE)
unlink("./package.pdf", recursive = TRUE)
unlink("./package.doc", recursive = TRUE)
Sys.sleep(2)


##### create package template #####
dir.create("package/")
if(.Platform$OS.type == "windows"){
  dir.create("package.zip/")
} else {
  dir.create("package.tar.gz/")
}

dir.create("package.doc/")
dir.create("package/R")
dir.create("package/data")
dir.create("package/vignettes")
#dir.create("package/src")
dir.create("package/man")
dir.create("package/inst")
dir.create("package/inst/extdata")
dir.create("package/inst/extdata/mixcr")
dir.create("package/inst/extdata/informat")



##### copy sources files #####
file.copy("./sources_files/DESCRIPTION", "./package/DESCRIPTION", overwrite=TRUE)
file.copy("./sources_files/NEWS","./package/NEWS", overwrite = TRUE)
file.copy("./sources_files/LICENSE","./package/LICENSE", overwrite = TRUE)
file.copy("./data/RepSeqData.rda", "./package/data/RepSeqData.rda", overwrite = TRUE)
file.copy("./vignettes/my-vignette.Rmd", "./package/vignettes/my-vignette.Rmd", overwrite = TRUE)


##### copy sources scripts #####
source.file <- list.files("./sources_scripts","\\.R$")
file.copy(paste0("./sources_scripts/", source.file), paste0("./package/R/",source.file), overwrite = TRUE)

#
# ##### copy src files #####
# source.file <- list.files("./sources_scripts.src")
# file.copy(paste0("./sources_scripts.src/", source.file), paste0("./package/src/",source.file), overwrite = TRUE)
#

##### copy extdata #####
source.file <- list.files("./extdata/", recursive = TRUE)
file.copy(paste0('./extdata/', source.file), paste0("./package/inst/extdata/", source.file), overwrite = TRUE)


##### build package and create doc #####
setwd("./package/")
usethis::use_namespace(roxygen = TRUE)
#Rcpp::compileAttributes()
devtools::document(roclets = c("namespace", "rd", "collate"))
devtools::build(binary = TRUE, args = c('--preclean'))
devtools::build()
setwd("../")

#devtools::document("./package")
if(.Platform$OS.type == "windows"){
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Program Files/MiKTeX/miktex/bin/x64", sep = .Platform$path.sep))
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/rtools40/usr/bin/", sep  = .Platform$path.sep))
}
system("R CMD Rd2pdf --no-preview package")


##### move zip file #####
pdffile = list.files("./", pattern = "*.pdf", full.names = TRUE)
pdffile = pdffile[grepl("package",pdffile)]
file.rename(pdffile, paste0("./package.doc/", basename(pdffile)))


##### move compressed files #####
if(.Platform$OS.type == "windows"){
	file = list.files("./", pattern = "*.zip", full.names = TRUE)
	file = file[!grepl("package",file)]
	file.rename(file,paste0("./package.zip/", basename(file)))
} else {
	file = list.files("./", pattern = "*.gz", full.names = TRUE)
	file = file[!grepl("package",file)]
	file.rename(file, paste0("./package.tar.gz/", basename(file)))
}

##### check #####
# devtools::check("./package", cran = FALSE)
# 
# if(.Platform$OS.type == "windows"){
# 	zipfile = list.files("./package.zip/", pattern = "*.zip", full.names = TRUE)
# } else {
# 	zipfile = list.files("./package.tar.gz/", pattern = "*.gz", full.names = TRUE)
# }
# 
# BiocCheck::BiocCheck(zipfile[1], "no-check-vignettes" = TRUE)
# BiocCheck::BiocCheck(zipfile[2], "no-check-vignettes" = TRUE)

###### install package #####
detach("package:AnalyzAIRR", unload = TRUE)
remove.packages("AnalyzAIRR")
Sys.sleep(2)
print("Installing package AnalyzAIRR")
if(.Platform$OS.type == "windows"){
	install.packages("./package.zip/AnalyzAIRR_1.0.0.zip", repos = NULL, type = "source")
} else{
	install.packages("./package.tar.gz/AnalyzAIRR_1.0.0.tar.gz", repos = NULL, type = "source")
}
print("Uploading AnalyzAIRR library to environment")
library(AnalyzAIRR)
print("Package AnalyzAIRR successfully uploaded")
