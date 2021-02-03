conflictRules("MASS", exclude = "select")
#install.packages('hapsim',exclude="select")
library('hapsim')

library('AlphaSimR')
library('tidyverse')
library('tidygraph')
library('ggpubr')
library('ggridges')
library('ggraph')
library(purrr)        # Functional programming
library(dplyr)        # Data wrangling
library(tidyr)        # Tidy-ing data
library(broom)        # List columns within tibbles
library(corrplot)
library(gridExtra)
library(data.table) #For quickly collapsing a list of data.frames
library(ggplot2)
require('Ghat')
require('rrBLUP')
library(latex2exp)
library(patchwork)


source("~/epistasis_in_breeding/code/Functions.R")
options(scipen = 1000000)
