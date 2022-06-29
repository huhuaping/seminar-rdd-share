# pkg needed
library('tidyverse')
library("rmarkdown")
library('bookdown')
library('knitr')
library("magrittr")  # for pipe

# package management
require('remotes')

# project tools
require("RCurl")    
require("openxlsx")
require("here")

# plot tool
library('ggplot2')
library('svglite')
require("latex2exp")

# data clean tool
require("janitor")
library("scales")

# data analysis
require("foreign")
require("haven")

# tool for publish
require("DT")

# my custom pkg
# remotes::install_github("KWB-R/kwb.utils")
# remotes::install_github("huhuaping/xmerit")
require("xmerit")

