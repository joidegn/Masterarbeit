library(xtable)

lines <- readLines("../data/bundesbank_series")
lines <- lines[lines!=""]
matr <- matrix(lines, ncol=2, byrow=T)
colnames(matr) <- c("Series Id", "Series Title")
sink("../data/bundesbank_series.tex")
print(xtable(matr))
sink(NULL)
