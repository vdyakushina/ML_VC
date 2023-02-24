# Get arguments
args <- commandArgs(T)
input <- args[1]
cov <- args[2]
pref <- args[3]

## TODO: comment after debug 
#input <- "/mnt/atlas/proc/20200428_130220/k1.mosdepth.region.dist.txt"
#cov <- 120

# Load data
d1 <- read.table(input)

# Subset 
d2 <- d1[rev(rownames(d1[which(d1$V1 == "total"),])),]

# Save to the file
png(paste0(pref, ".cov.png"), h = 800, w = 800, pointsize = 20)

# Create plot area. Add gridlines and axis labels.
plot(d2$V2, d2$V3, 
     type = 'n', xlab = "Depth", 
     ylab = "Fraction of capture target bases \u2265 depth", 
     xlim= c(0, 800), ylim = c(0, 1.0))
mtext(side = 3, line = 3, "Target Region Coverage")
mtext(side = 3, line = 2, paste("Average coverage:", cov))
abline(v = cov, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(2, at = c(0.90), labels = c(0.90))
axis(2, at = c(0.50), labels = c(0.50))

# Plot the data 
points(d2$V2, d2$V3, type = 'l', lwd = 2, col = "red")

# Add a legend using the sample names
legend("topright", legend = pref, col = "red", lty = 1, lwd = 1)

invisible(dev.off())


