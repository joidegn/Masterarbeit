data <- read.csv('data/1959-2014_normalized.csv')
data <- data[, -1]
y <- data[, 1]
x <- data[5:dim(data)[1], -1]

res <- princomp(x)
loadings <- res$loadings
