library("forecast")

monthly.data <- read.csv("data/1959-2014.csv")

monthly.data <- monthly.data[1:nrow(monthly.data)-1,]  # last line is NA for almost all series
nas <- apply(monthly.data, 2, function(col) sum(is.na(col)))
nas[1] <- 1 # first column is data column and we dont want that

monthly.noNas <- monthly.data[, nas==0]

demean <- function(series) series-mean(series)
rescale <- function(series) series/sd(series)
make.stationary <- function(series) {
    temp <- series
    ndiffs <- ndiffs(series, alpha=0.05, test="kpss")
    if (ndiffs > 0)
        for (i in 1:ndiffs)  # TODO: ADF seems to be more conservative than KPSS
            temp <- diff(temp)
    return(temp)
}

monthly.stationary <- apply(monthly.noNas, 2, make.stationary)
min.length <- Reduce(min, lapply(monthly.stationary, length))
# TODO: what to do with the longer data series? Cut off last elements? cut-off first? make missing points median of the other few?
monthly.clean <- as.data.frame(lapply(monthly.stationary, function(series) rescale(demean(series[1:min.length]))))

write.csv(monthly.clean, "data/1959-2014_normalized.csv")
