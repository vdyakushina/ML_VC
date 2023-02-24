## Estimate the eveness score
## The reference: https://www.nature.com/articles/jhg201621
# Initiate
args <- commandArgs(T)
fn <- args[1]
# Check file exist
if(!file.exists(fn)) stop(fn, "doesn't exist", .call = F)
# Load data
data <- read.table(fn, stringsAsFactors = F)
# Estimate the eveness score
D <- data$V7
C <- round(mean(D))
D2 <- D[D<=C]
E <- 1 - (length(D2)-sum(D2)/C)/length(D)
# Output the result
cat(round(E, 2))
