## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE,warning = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(np)
library(stats)
library(MASS)
library(kdml)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  library(devtools)
#  install_github("jrjthompson/R-package-kdml")
#  library(kdml)

## -----------------------------------------------------------------------------
df <- data.frame(
  x1 = runif(100, 0, 100),
  x2 = factor(sample(c("A", "B", "C"), 100, replace = TRUE)),
  x3 = factor(sample(c("A", "B", "C"), 100, replace = TRUE)),
  x4 = rnorm(100, 10, 3),
  x5 = ordered(sample(c("Low", "Medium", "High"), 100, replace = TRUE), 
               levels = c("Low", "Medium", "High")),
  x6 = ordered(sample(c("Low", "Medium", "High"), 100, replace = TRUE), 
               levels = c("Low", "Medium", "High"))
)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # DKPS distance (Equation 1 and 3)
#  dis_dkps <- dkps(df = df)
#  
#  # DKSS distance (Equations 2 and 3)
#  dis_kdss <- kdss(df = df)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # DKPS distance (Equation 1 and 3)
#  dis_dkps_np <- dkps(df = df, bw = "np")
#  
#  # DKSS distance (Equations 2 and 3)
#  dis_kdss_np <- kdss(df = df, bw = "np")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  dis_dkps_custom_kernels <- dkss(df = df, bw = "mscv",
#      cFUN = "c_epanechnikov", uFUN = "u_aitken", oFUN = "o_habbema")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # MSCV bandwidth specification using the similarity function in Equation (1)
#  mscv.dkps(df, nstart = NULL, ckernel = "c_gaussian",
#             ukernel = "u_aitken", okernel = "o_wangvanryzin", verbose = TRUE)
#  
#  # MSCV bandwidth specification using the similarity function in Equation (2)
#  mscv.dkss(df, nstart = NULL, ckernel = "c_gaussian",
#            ukernel = "u_aitken", okernel = "o_wangvanryzin", verbose = TRUE)

