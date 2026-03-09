library(ape)
library(magrittr)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)

args <- commandArgs(TRUE)
infile <- args[1]
outname <- args[2]

data <- read.table(infile, header = TRUE, sep = "\t", row.names = 1)

## ----------------------------
## 1. set outgroup
## ----------------------------
if (any(grepl("_W", colnames(data)))) {
  outgroups <- colnames(data)[grepl("_W", colnames(data))]
} else {
  stop("No outgroup (_W) samples found in column names")
}

## ----------------------------
## 2. build tree
## ----------------------------
tree <- nj(as.dist(1 - cor(data, method = "spearman", use = "everything")))
tree <- root(tree, outgroup = outgroups, resolve.root = TRUE)

pdf(paste0(outname, ".tree.pdf"))
plot(tree, cex = 1, edge.width = 2, label.offset = 0.01)
dev.off()

write.tree(tree, paste0(outname, ".tree.nwk"))

## ----------------------------
## 3. Bootstrap + treelength calculation（Normalized）
## ----------------------------
bootstrap.reps <- 1000
set.seed(123)

mytrees <- vector("list", bootstrap.reps)

branch_lengths_raw <- numeric(bootstrap.reps)
branch_lengths_norm <- numeric(bootstrap.reps)
ntips_vec <- numeric(bootstrap.reps)

for (i in 1:bootstrap.reps) {

  sampled.data <- data[
    sample(1:nrow(data), size = nrow(data), replace = TRUE),
  ]

  this.tree <- nj(as.dist(1 - cor(sampled.data, method = "spearman")))

  this.tree <- root(
    this.tree,
    outgroup = outgroups,
    resolve.root = TRUE
  )

  mytrees[[i]] <- this.tree

  total_branch_length <- sum(this.tree$edge.length)
  n_tips <- Ntip(this.tree)

  branch_lengths_raw[i] <- total_branch_length
  branch_lengths_norm[i] <- total_branch_length / (n_tips - 1)
  ntips_vec[i] <- n_tips
}

branch_lengths_df <- data.frame(
  Tree_ID = 1:bootstrap.reps,
  N_tips = ntips_vec,
  Branch_Length_Raw = branch_lengths_raw,
  Branch_Length_Normalized = branch_lengths_norm
)

write.table(
  branch_lengths_df,
  paste0(outname, "_branch_lengths_raw_and_normalized.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

## ----------------------------
## 4. Choose bootstrap tree
## ----------------------------
parts <- prop.part(mytrees, check.labels = TRUE)

Sum <- numeric(bootstrap.reps)

for (i in 1:bootstrap.reps) {
  boot <- prop.clades(
    mytrees[[i]],
    part = parts,
    rooted = TRUE
  )
  Sum[i] <- sum(boot, na.rm = TRUE)
}

Sum_ID <- data.frame(
  Tree_ID = 1:bootstrap.reps,
  Sum = Sum
) %>% arrange(desc(Sum))

best_id <- Sum_ID$Tree_ID[1]
trw <- mytrees[[best_id]]

## ----------------------------
## 5. Calculate and color the bootstrap rate
## ----------------------------
boot <- prop.clades(
  trw,
  part = parts,
  rooted = TRUE
)

trw$boot.color <- rep("red", length(boot))

for (i in seq_along(boot)) {

  trw$label_[i] <- boot[i] / bootstrap.reps

  if (!is.na(boot[i]) && boot[i] / bootstrap.reps > 0.9) {
    trw$boot.color[i] <- "white"
  } else if (boot[i] / bootstrap.reps > 0.7) {
    trw$boot.color[i] <- "yellow"
  } else if (boot[i] / bootstrap.reps > 0.5) {
    trw$boot.color[i] <- "goldenrod1"
  } else {
    trw$boot.color[i] <- "red"
  }
}

## ----------------------------
## 6. Plots
## ----------------------------
pdf(paste0(outname, ".spearman.rooted.boot.pdf"))

plot(
  trw,
  main = "Rooted Phylogenetic Tree Based on Spearman Correlation",
  cex = 0.5
)

add.scale.bar(x = 0.0, y = 1.3, length = 0.1)

nodelabels(
  pch = 21,
  col = "black",
  bg = trw$boot.color,
  cex = 1
)

dev.off()

