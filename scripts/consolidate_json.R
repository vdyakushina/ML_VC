# Initiate
## TODO: uncomment for debug
# input <- "/mnt/atlas/home/gkhvorykh/proc/20200428_130220"
# project <- "20200428_130220"

args <- commandArgs(T)

input <- args[1]
project <- args[2]

json_to_dataframe <- function(x){
  # Read all json files
  res <- rjson::fromJSON(file = x)
  df <- as.data.frame(res, stringsAsFactors = F)
  id <- df$sample[1]
  colnames(df) <- gsub(paste0(id, "_"), "", colnames(df))
  # Correct colnames
  ind <- grepl(".1$", colnames(df))
  colnames(df)[ind] <- gsub("pool", "pool2", colnames(df)[ind])
  colnames(df)[ind] <- gsub(".1$", "", colnames(df)[ind])
  colnames(df) <- gsub("pool\\.", "pool1.", colnames(df))
  # Drop redundent columns
  df <- df[!grepl("pool1.number|pool2.number", colnames(df))]
  df
}

# Get list of json files
lf <- list.files(input, "*.json", recursive = T, full.names = T)
message("json files found: ", length(lf))
out <- lapply(lf, json_to_dataframe)
out <- do.call(what = "rbind", args = out)

# Save as csv
write.csv(out, sprintf("%s/%s.csv", input, project), quote = F,
          row.names = F)


