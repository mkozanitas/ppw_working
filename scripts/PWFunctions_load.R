## LOAD PPW FUNCTIONS
# clear workspace
rm(list=ls())

# The following packaged need to be installed the following functions to work: "RCurl", "data.table","picante"
library("RCurl")
library("data.table")
library("picante")

# this code is from some delightful human on the internet who figured out how to source from GitHub; run this function in order to source in the script with all the Pepperwood Functions 
source_https <- function(url, ...) { 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

# sources in all functions (described below) that allow access to the PPW Vegetation Plot Data
source_https('https://raw.githubusercontent.com/dackerly/PepperwoodVegPlots/master/Analyses/PWfunctions_GitHub.R')
