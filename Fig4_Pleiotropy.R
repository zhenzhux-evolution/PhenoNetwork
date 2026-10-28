rm(list=ls())
library(data.table)
library(dplyr)
library(ComplexHeatmap) #install.packages("BiocManager") BiocManager::install("ComplexHeatmap")
library(vegan)
library(grid)


setwd("/home/zhenzhux/pheno/gwas/PolygScore/")
out_dir = "/home/zhenzhux/pheno"


# make the list from all significant snps for all traits
files <- list.files(path = "./traitPlot_28cov", pattern = "^sig_Trait[0-9]{2}\\.tsv$", full.names = TRUE)
trait_list <- setNames(
  lapply(files, \(f) fread(f)$rs),
  sub("^sig_(Trait[0-9]{2})\\.tsv$", "\\1", basename(files))
)
# binary matrix + pleiotropy snps
biMatrix <- list_to_matrix(trait_list)
biMatrix_pleio <- biMatrix[rowSums(biMatrix) >= 2, , drop = FALSE]
biMatrix_pleio <- biMatrix_pleio[, colSums(biMatrix_pleio) >= 1, drop = FALSE]
#biMatrix_pleio <- biMatrix_pleio[, colnames(biMatrix_pleio) != "Trait19"]
write.table(biMatrix_pleio, "biMatrix_pleio.txt", quote = F, sep = "\t", col.names = T, row.names = T)

colnames(biMatrix_pleio)
traits_morph <- sprintf("Trait%02d", 25:40); 
traits_voc <- sprintf("Trait%02d", 11:24)

cols_morph <- intersect(colnames(biMatrix_pleio), traits_morph)
cols_voc   <- intersect(colnames(biMatrix_pleio), traits_voc)
n_morph <- length(cols_morph)
n_voc   <- length(cols_voc)

## Row sums per group
combined_matrix <- data.table(
  SNP        = rownames(biMatrix_pleio),
  sum_morph  = rowSums(biMatrix_pleio[, cols_morph, drop = FALSE]),
  sum_voc    = rowSums(biMatrix_pleio[, cols_voc,   drop = FALSE])
)

## Fractions + observed difference
combined_matrix[, `:=`(
  Frac_morph = sum_morph / n_morph,
  Frac_voc   = sum_voc   / n_voc,
  obs_diff   = abs((sum_morph / n_morph) - (sum_voc / n_voc))
)]
quantile(combined_matrix$obs_diff)


######### shuffle the trait labels 
# function to compute the fraction difference
compute_FracDiff <- function(mat, cols_morph, cols_voc, n_morph, n_voc) {
  abs(
    rowSums(mat[, cols_morph, drop = FALSE]) / n_morph - 
      rowSums(mat[, cols_voc, drop = FALSE]) / n_voc
  )
}
obs_diff <- combined_matrix$obs_diff
names(obs_diff) <- combined_matrix$SNP
obs_diff_mean <- mean(obs_diff)

n_perm <- 1000
perm_res <- permatfull(
  biMatrix_pleio,
  fixedmar = "rows", mtype = "prab", #for presence-absence type incidence data
  times = n_perm
)
perm_mats <- perm_res$perm
head(perm_mats[[3]])
# apply to all the simulated lists
null_diff <- sapply(
  perm_mats,
  compute_FracDiff,
  cols_morph = cols_morph,
  cols_voc   = cols_voc,
  n_morph    = n_morph,
  n_voc      = n_voc
)
null_diff[1:2, 1:5]

# mean of null diff across all SNPs 
null_diff_mean <- colMeans(null_diff)
p_mean <- sum(null_diff_mean >= obs_diff_mean) / length(null_diff_mean)
mean(null_diff_mean)
summary(null_diff_mean)
obs_diff_mean  # 0.447805


colSums(biMatrix)
pdf(file = file.path(out_dir, "hist_meandiff_1k.pdf"), width = 3, height = 4, bg = "white")
#png(filename = file.path(out_dir, "hist_meandiff_1k.png"),
#    width = 3, height = 4, units = "in", res = 600) #, bg = "white"
hist(null_diff_mean, breaks = 4, main = NA, xlab = "Mean difference", xlim = c(0, 0.5))
abline(v = mean(obs_diff_mean), col = "blue", lty = 2, lwd = 2)
dev.off()



biMatrix_trait_pleiotropy <- biMatrix_pleio 
m <- make_comb_mat(biMatrix_trait_pleiotropy)
m_pleiotropy = m[comb_degree(m) > 1]
set_order_trait <- sort(set_name(m_pleiotropy))

pdf(file = file.path(out_dir, "intersectPlot_pleiotropy.pdf"), 
    width = 12, height = 8, bg = "white")
#png(filename = file.path(out_dir, "intersectPlot_pleiotropy.png"),
#    width = 18, height = 12, units = "cm", res = 600) #, bg = "white"
ht <- UpSet(
  m_pleiotropy,
  set_order = set_order_trait,
  comb_order = order(comb_size(m_pleiotropy), decreasing = TRUE),
  pt_size = unit(5, "mm"),
  lwd = 3,
  comb_col = c("grey50", "#39a339ff", "#47d7e7ff",
               rep("#fa841cff", 7))[comb_degree(m_pleiotropy)],
  top_annotation = upset_top_annotation(
    m_pleiotropy,
    add_numbers = TRUE,
    annotation_name_rot = 90,
    labels_gp = gpar(fontsize = 50)
  ),
  right_annotation = upset_right_annotation(
    m_pleiotropy,
    add_numbers = TRUE
  )
)

draw(ht)

grid.lines(
  x = unit(c(0, 1), "npc"),
  y = unit(c(0.375, 0.375), "npc"),   # middle of the page
  gp = gpar(col = "grey20", lwd = 3, lty = 2)
)
dev.off()
