library(forecast)
library(xts)
library(zoo)

normalize.data <- function(data) {  # takes csv file which will be demeaned, normalized and made stationary
    demean <- function(series) series-mean(series, na.rm=T)
    rescale <- function(series) series/sd(series, na.rm=T)
    make.stationary <- function(series, times) {
        ndiffs <- ndiffs(series, alpha=0.05, test="kpss")
        if (ndiffs > 0)
            for (i in 1:ndiffs)  # TODO: ADF seems to be more conservative than KPSS?
                series <- diff(series)
        return(zoo(series, times[(ndiffs+1):length(times)]))  # return a zoo object and keep the time index intactrtaking account for the number of differences
    }

    if (dim(data)[2] > 1)
        stationary.data <- do.call(merge, apply(data, 2, function(column) rescale(demean(make.stationary(column, names(column))))))
        #library(xtable)
        #sink("../data/differences.tex")
        #diffs <- matrix(sapply(data, function(col) ndiffs(col, alpha=0.05, test="kpss")), ncol=1)
        #colnames(diffs) <- colnames(data)
        #xtable(diffs)
        #sink(NULL)
    else
        stationary.data <- rescale(demean(make.stationary(data[,1], index(data))))
    # TODO: choose data range!
    return(stationary.data)
}

sanitize.dates <- function(data.frm) {  # checks if dates are in format "2010-10" (i.e. %Y-%m) or "2010-10-13" (i.e. %Y-%m-%d)
    # first drop rows which dont appear to hold dates (testing is too simple and might fail miserably)
    valid.date.rows <- sapply(as.character(data.frm[,1]), function(datestr) ifelse(length(strsplit(datestr, "-")[[1]]) %in% c(2,3), T, F))
    data.frm <- data.frm[valid.date.rows, ]
    data.frm[, 1] <- sapply(as.character(data.frm[,1]), function(datestr) ifelse(length(strsplit(datestr, "-")[[1]])==2, paste0(datestr, "-01"), datestr))
    return(data.frm)
}

remove.short.series <- function(data, min.length=100) {  # removes columns which have too few non NA values
    lengths <- apply(data, 2, function(column) length(column) - sum(is.na(column)))
    return(data[, which(lengths>min.length)])
}

cut.to.size <- function(data) {  # cuts whole data.frame or zoo thing to drop NA margins
    start.index <- max(apply(data, 2, function(column) for (i in 1:length(column)) if (!is.na(column[i])) return(i-1)))  # count from bottom up until NA
    end.index <- min(apply(data, 2, function(column) for (i in length(column):1) if (!is.na(column[i])) return(i+1)))  # count from top down until NA
    if (end.index-start.index <= 0) {
        warning(" no observations left!")
        return(data.frame())
    }
    return(data[start.index:end.index, ])
}

remove.series.that.cost.observations <- function(data, max.margin.costs=c(10,10)) {  # some series have a bad time frame making it such that cut.to.size removes almost all observations
    start.indices <- apply(data, 2, function(column) for (i in 1:length(column)) if (!is.na(column[i])) return(i-1))  # count from bottom up until NA
    end.indices <- apply(data, 2, function(column) for (i in length(column):1) if (!is.na(column[i])) return(i+1))  # count from top down until NA
    leftover <- data[, -unique(c(which(start.indices > max.margin.costs[2]+1), which(dim(data)[1] - end.indices > max.margin.costs[2])))]
    return(leftover)
}

make.data.csv <- function(files, output.file=F) {  # merges all csv files into one by date, interpolates, demeans and normalized the data
    raw.data.sets <- lapply(files, function(file) sanitize.dates(read.csv(file, sep=" ", header=F)))  # read.zoo needs data everywhere and is oppinionated about csv structure
    # naive filter: checks if any dates are duplicates (this includes NAs and similar shenanigans)
    # update date column to be of the format 2010-10-13 i.e. %Y-%m-%d
    keep.data.sets.idx <- sapply(raw.data.sets, function(data.set) anyDuplicated(as.Date(as.character(data.set[, 1]), format="%Y-%m-%d"))==0)
    valid.data.sets = raw.data.sets[keep.data.sets.idx]
    #titles <- sapply(valid.data.sets, function(data.set) as.character(data.set[1,2])) # works only on bundesbank
    ids <- sapply(files[keep.data.sets.idx], function(file) tail(strsplit(strsplit(file, ".csv")[[1]], "/")[[1]], n=1))
    # the following probably creates warnings because there are NAs
    data.sets <- lapply(valid.data.sets, function(data.set) zoo(data.set[, 2], as.Date(as.character(data.set[, 1]), format="%Y-%m-%d")))
    na.sets <- sapply(data.sets, function(data.set) all(is.na(data.set)))  # if there is one, it should be inspected by hand why this happened
    if (sum(na.sets) >0 ) warning(paste0("There are ", length(na.sets), " sets with no data."))
    # filter again, this time for any NAs. If there are Nas, cut them out of the series. Yes exactly as brutally as it sounds: with a knife! This is ok because many csv have an empty line at the end
    data.sets <- lapply(data.sets, function(data.set) na.omit(data.set))

    data <- do.call(merge, data.sets)
    approximated.data <- na.approx(data, na.rm=F)  # interpolate interior NAs (this refers especially to columns which were quarterly)
    cat(sum(is.na(data)) - sum(is.na(approximated.data)), " NA values have been approximated\n")
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



# the following was useful for debugging
#col.is.equal <- function(col, matr)
#    sapply(1:(dim(matr)[2]), function(idx) all(col %in% matr[, idx]))
#
#sapply(data.sets, function(data.set) which(col.is.equal(data.set, data)))

bb.files <- Sys.glob("../data/bundesbank_selected/*.csv")
fred.files <- Sys.glob("../data/FRED_selected/*.csv")

all.files <- c(bb.files, fred.files)

# I am only interested in the dates.
data <- make.data.csv(all.files)
selected.data <- window(data, start = as.Date("1991-01-01"))  # arbitrarily
long.enough <- remove.series.that.cost.observations(remove.short.series(selected.data, 200), c(20,20))
final.data <- cut.to.size(long.enough)
final.data <- xts(final.data, as.Date(as.yearmon(as.character(index(final.data)))))

gdp <- read.csv("../data/BB-GDP_quarterly.csv", sep=" ", header=F)
gdp.time.index <- as.Date(as.yearmon(as.character(gdp[,1])))
gdp <- xts(gdp[,2], gdp.time.index)
colnames(gdp) <- c("BBK01.JQA000")
final.data <- na.trim(normalize.data(na.omit(merge(zoo(gdp), zoo(final.data)))))
#write.zoo(final.data, "../data/final_data.csv")

series <- names(final.data)
bb.descriptions.mat <- read.csv("../data/series_bundesbank_descriptions", as.is=T, sep="\t", header=F)
bb.descriptions <- bb.descriptions.mat[,2]
names(bb.descriptions) <- bb.descriptions.mat[,1]
fred.descriptions.mat <- read.csv("../data/series_FRED_descriptions", as.is=T, header=F)
fred.descriptions = fred.descriptions.mat[,2]
names(fred.descriptions) <- fred.descriptions.mat[,1]
descriptions <- c(bb.descriptions, fred.descriptions)

variable.descriptions <- na.omit(descriptions[series])
nas = sapply(series[which(is.na(descriptions[series]))], function(id) system(paste0("cat ../data/items.grouped.json | jq '.[] | .[] | select(.id==\"", id, "\") | .title' | sed 's/\"//g'"), intern=T))
# not enough time to solve this issue properly, we just grep the missing descriptions from the file. This is clumsy, slow and require jq (the command line json parser to be installed)
variable.descriptions <- c(variable.descriptions, nas)
#library(xtable)
#sink("variable_table.tex")
#xtable(matrix(c(names(variable.descriptions), variable.descriptions), ncol=2))
#sink(NULL)
 
