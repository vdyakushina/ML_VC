# Inititate
args <- commandArgs(T)
id <- args[1]
path <- c("samtools" = "/pipeline/tools/samtools-1.11/bin/samtools")

# TODO: comment after debug
# bam <- "/mnt/atlas/home/gkhvorykh/proc/20200428_130220/k1.bam"
# path <- c("samtools" = "/home/gennady/tools/samtools-1.12/bin/samtools")
# id <- "k1"

# Check input
bam <- paste0(id, ".bam")
if (! file.exists(bam)) stop(bam, "doesn't exist", .call = F)

# Estimate the coverage by chromosomw
out <- system2(path['samtools'], c("idxstats", bam), stdout = T)

# Convert in data frame
data <- strsplit(out, "\t")
data <- do.call(what = "rbind", args = data)
df <- data.frame(data, stringsAsFactors = F)
df <- df[-26, ]

# Plot
png(paste0(id, ".chr.cov.png"), width = 1200, height = 900, units = "px")
barplot(as.numeric(df$X3)/1000, names.arg = df$X1, 
        las = 2, ylab = "Number of reads, K")
title = "Chromosome coverage"
subtitle = paste("Sample:", id)
mtext(side = 3, line = 3, cex = 1.2, title)
mtext(side = 3, line = 2, cex = 1, subtitle)
invisible(dev.off())
