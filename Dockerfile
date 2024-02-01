FROM rocker/r-ver:4.3.2
RUN apt update
RUN apt-get install -y libssl-dev libxml2-dev libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev zlib1g-dev
RUN apt-get install -y cmake
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_github("vanessajmh/AnalyzAIRR")'
RUN Rscript -e 'install.packages('devtools',repos = 'http://cran.us.r-project.org')'
RUN Rscript -e 'devtools::test()'

#RUN Rscript -e 'testthat::test_package("AnalyzAIRR")'