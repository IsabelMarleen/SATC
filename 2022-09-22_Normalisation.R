# Read in data
#filename <- "/home/genis/projects/satc2/leopards/test_leopardwin100kb.txt"
filename <- "examples/test_leopardwin100kb.txt"
raw <- read.table(filename, header=T, as.is=T, row.names = 1)

# First normalisation step - normalise for scaffold length
norm1 <- apply(raw, 2, function(x) x/raw[, "length"])
stopifnot(all(norm1[, "length"] ==1))

# Second normalisation step - normalise for library size
norm2 <- t(apply(norm1[,-1], 1, function(x) x/median(x)))
