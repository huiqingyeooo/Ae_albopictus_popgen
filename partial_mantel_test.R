# Partial Mantel Test

rm(list = ls())
library(adegenet)
library(vcfR)
library(poppr)
library(ncf)

genind <- import2genind("albo_structure.str")

coords <- read.csv("b5_sorted.csv") #ensure that order of samples is same as genind
other(genind)$xy <- coords
genind$other$xy
Dgeo <- dist(other(genind)$xy)
Dgen <- provesti.dist(genind$tab)

lai <- read.csv("lai.csv")
lai_values <- lai[,c(2)]
Dlai <- dist(lai_values)

MDgen <- as.matrix(Dgen)
MDgeo <- as.matrix(Dgeo)
MDlai <- as.matrix(Dlai)

dim(MDgen) #check dimensions are the same
dim(MDgeo)
dim(MDlai)
partial.mantel.test(MDgen, MDlai, MDgeo, resamp = 1000, method = "pearson", quiet = FALSE)


