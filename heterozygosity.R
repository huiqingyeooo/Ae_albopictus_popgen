# Calculate basic popgen statistics

rm(list=ls())
library(PopGenome)
library(diveRsity)
library(HWxtest)
library(Rcpp)
library(adegenet)
library(poppr)
library(hierfstat)

genind <- import2genind("albo_structure.str")

###########################################
# Expected and observed heterozygosity
############################################
# basic stats across entire population

basic_stat <- basic.stats(genind, diploid=TRUE) #calculate basic stats
basic_stat #mean of Ho, Hs, Ht and FIS

Ho <- basic_stat$Ho
Ho <- as.data.frame(Ho)
stderror <- function(x) sd(x)/sqrt(length(x))
stderror(Ho[,1]) #stderror

Hs <- basic_stat$Hs
Hs<-as.data.frame(Hs)
stderror(Hs[,1]) #stderror

###########################################
# FIS (Inbreeding coefficient)
############################################
Fis <- basic_stat$Fis
Fis <-as.data.frame(Fis)
stderror(Fis[,1]) #stderror