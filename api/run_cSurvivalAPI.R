#!/usr/bin/env Rscript

library(stringr)
library(hash)
library(plumber)

pr("data_file_api.R") %>%
  pr_run(port=4000, host="0.0.0.0")
# Setting the host option on a VM instance ensures the application can be accessed externally.
# (This may be only true for Linux users.)