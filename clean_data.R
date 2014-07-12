library(forecast)
library(zoo)

clean.data <- function(data) {  # takes csv file which will be demeaned, normalized and made stationary
    demean <- function(series) series-mean(series, na.rm=T)
    rescale <- function(series) series/sd(series, na.rm=T)
    make.stationary <- function(series, times) {
        ndiffs <- ndiffs(series, alpha=0.05, test="kpss")
        if (ndiffs > 0)
            for (i in 1:ndiffs)  # TODO: ADF seems to be more conservative than KPSS?
                series <- diff(series)
        return(zoo(series, times[(ndiffs+1):length(times)]))  # return a zoo object and keep the time index intactrtaking account for the number of differences
    }

    stationary.data <- do.call(merge, apply(data, 2, function(column) rescale(demean(make.stationary(column, names(column))))))
    # TODO: choose data range!
    return(stationary.data)
}

make.data.csv <- function(files, output.file=F) {  # merges all csv files into one by date, interpolates, demeans and normalized the data
    raw.data.sets <- lapply(files, function(file) read.csv(file))
    # naive filter: checks if any dates are duplicates (this includes NAs and similar shenanigans)
    keep.data.sets.idx <- sapply(raw.data.sets, function(data.set) anyDuplicated(as.Date(paste0(as.character(data.set[5:dim(data.set)[1], 1]), "-01"), format="%Y-%m-%d"))==0)
    valid.data.sets = raw.data.sets[keep.data.sets.idx]
    #titles <- sapply(valid.data.sets, function(data.set) as.character(data.set[1,2])) # works only on bundesbank
    ids <- sapply(files, function(file) tail(strsplit(strsplit(file, ".csv")[[1]], "/")[[1]], n=1))
    # the following probably creates warnings because there are NAs
    data.sets <- lapply(valid.data.sets, function(data.set) zoo(as.numeric(levels(data.set[5:dim(data.set)[1], 2])[data.set[5:dim(data.set)[1], 2]]), as.Date(paste0(as.character(data.set[5:dim(data.set)[1], 1]), "-01"), format="%Y-%m-%d")))
    # filter again, this time for any NAs. If there are Nas, cut them out of the series. Yes exactly as brutally as it sounds: with a knife! This is ok because many csv have an empty line at the end
    data.sets <- lapply(data.sets, function(data.set) na.omit(data.set))

    data <- do.call(merge, data.sets)
    approximated.data <- na.approx(data, na.rm=F)  # interpolate interior NAs (this refers especially to columns which were quarterly)
    cat(sum(is.na(data)) - sum(is.na(approximated.data)), " NA values have been approximated\n")
    data <- clean.data(approximated.data)  # demean, rescale, make stationary, interpolate missings
    colnames(data) <- ids
    if (is.character(output.file))
        write.zoo(data, output.file)  # unfortunately whitespace delimited
    return(data)
}

split.csv <- function(file, folder=F) { # this splits a file into different files named after the column
    data = read.zoo(file, format="%Y-%m-%d", header=T, sep=",")
    datas <- lapply(data, function(dta) dta)
    if (is.character(folder))
        for (series.idx in 1:length(datas))
            write.zoo(datas[[series.idx]], paste0(folder, "/", names(datas)[series.idx], ".csv"))
}



# the following might be useful for debugging
#col.is.equal <- function(col, matr)
#    sapply(1:(dim(matr)[2]), function(idx) all(col %in% matr[, idx]))
#
#sapply(data.sets, function(data.set) which(col.is.equal(data.set, data)))

files <- Sys.glob("data/bundesbank_selected/*.csv")
files2 <- Sys.glob("data/FRED_selected/*.csv")

all.files <- c(files, files2)
