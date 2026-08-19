args <- commandArgs(trailing = TRUE)
name <- tools::file_path_sans_ext((basename(args[1])))
frq_path <- args[1]
idepth <- read.table(args[2], header = TRUE)
imiss <- read.table(args[3], header = TRUE)
ldepth_mean <- read.table(args[4], header = TRUE)
lqual <- read.table(args[5], header = TRUE)
lmiss <- read.table(args[6], header = TRUE)
het <- read.table(args[7], header = TRUE)
hwe <- read.table(args[8], header = TRUE)

# Count greatest allele number to format frq table for MAF
frq_header_index <- 1
frq_field_counts <- count.fields(frq_path)[-frq_header_index]
frq_n_metadata_columns <- 4
frq_max_alleles <- max(frq_field_counts, na.rm = TRUE) - frq_n_metadata_columns
frq_col_names <- c(
  "CHROM",
  "POS",
  "N_ALLELES",
  "N_CHROMS",
  paste0("A", seq_len(frq_max_alleles))
)

frq <- read.table(
  frq_path,
  skip = frq_header_index,
  header = FALSE,
  fill = TRUE,
  col.names = frq_col_names,
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

frq$MAF <- (frq[grep("^A", names(frq), value = TRUE)] |> apply(1, \(x) min(x)))
frq$MAC <- round(frq$MAF * frq$N_CHROMS) # Not plotted. Account for missingness

png(
  filename = paste(name, ".png", sep = ""),
  width = 6.75, height = 9, res = 600,
  units = "in"
)
par(mfrow = c(4, 2))

hist(het$F, main = "Inbreeding Coefficient (F)", xlab = "")
hist(idepth$MEAN_DEPTH, main = "Mean Depth (Ind.)", xlab = "")
hist(imiss$F_MISS, main = "Missing Data (Ind.)", xlab = "Fraction Missing")

cutoff <- ldepth_mean$MEAN_DEPTH <= quantile(ldepth_mean$MEAN_DEPTH, 0.99)
hist(
  ldepth_mean$MEAN_DEPTH[cutoff],
  main = "Mean Depth (Site; 0th - 99th Percentile)",
  xlab = ""
)

cutoff <- lqual$QUAL <= quantile(lqual$QUAL, 0.99)
hist(
  lqual$QUAL[cutoff],
  main = "Quality (Site; 0th - 99th Percentile)",
  xlab = ""
)

hist(lmiss$F_MISS, main = "Missing Data (Site)", xlab = "Fraction Missing")
hist(frq$MAF, main = "Minor Allele Frequency (Site)", xlab = "")

hist(
  -log10(hwe$P_HWE),
  main = "Hardy-Weinberg Test (Site)",
  xlab = "-log10(P)"
)

dev.off()
