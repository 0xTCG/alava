#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
print(args)

# 3. This will generate a file which contains F coefficient estimates for assessing heterozygosity.
# We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean,
# which can be performed using the following R commands:

dat <- read.table(args[1], header = TRUE, comment.char = "")

m <- mean(dat$F) # Calculate the mean
s <- sd(dat$F) # Calculate the SD
print(c(m, s))

# Get any samples with F coefficient within 3 SD of the population mean
valid <- subset(dat, F <= m + 3 * s & F >= m - 3 * s)
write.table(valid[, c(1, 2)], args[2], quote = FALSE, row.names = FALSE)
