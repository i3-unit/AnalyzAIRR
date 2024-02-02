FROM rocker/r-ver:4.3.2
RUN apt update
RUN apt-get install -y libssl-dev libxml2-dev libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev zlib1g-dev
RUN apt-get install -y cmake
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_github("vanessajmh/AnalyzAIRR", INSTALL_opts= "--install-tests")'
#RUN Rscript -e 'testthat::test_dir(path= paste0(.libPaths()[1], "/AnalyzAIRR/tests/"))'
#RUN Rscript -e 'library(AnalyzAIRR)'
#RUN Rscript -e 'testthat::test_package("AnalyzAIRR")'