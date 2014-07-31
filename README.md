Replication files for "Forecasting Economic Time Series Using Dynamic Factor Models under Structural Breaks"
====
### Requirements
The calculations for this master thesis have been done mostly in [R] (http://cran.r-project.org/) and [Julia] (http://julialang.org/). For R the packages that have been used are 
* forecast
* xts
* zoo
* xtable (only used to extract latex tables)

For Julia the list of required Packages is:
* FactorModels (written for this thesis)
* Distributions
* Gadfly (only used for graphs)
* DataFrames
* DateTime
* DataStructures


### Folder structure
The  main folder holds all program code required for replication. Namely the files are `prediction.jl`, `bootstrap.jl` and `break_removal_montecarlo.jl`  which can best be run in a  [Julia REPL] (http://julialang.org/). The data folder holds all data files including the raw data from the bundesbank and from FRED.

### Downloading, choosing and cleaning the data
Some of the data series have been downloaded from the bundesbank website using Python and [Scrapy] (http://scrapy.org/). The program code to download the data from the bundesbank website can be found in the folder data/scrape/bundesbank. The results are on the one hand the csv files in data/bundesbank and on the other hand a list of all series descriptions in data/scrape/bundesbank/items.json which has been prettyfied in /data/scrape/bundesbank/items.pretty.json.

The rest of the data comes from the Federal Reserve of St. Louis (specifically the FRED2 API) and has been downloaded straightforwardly using the code in the [FRED] (https://github.com/joidegn/FRED.jl) Julia package.  The result can be found in data/FRED\_selected/.

Shell scripts, R and Julia code to clean, normalize, difference and merge the data can be found in the data\_scripts folder. The cleaned and normalized data set on which the calculations are run can be found at data/final\_data.csv.


### Replication of the Results
Most easily the results can be replicated by installing the necessary dependencies for julia, checking out this repository and running the main `prediction.jl`, `bootstrap.jl` and `break_removal_montecarlo.jl` in a julia REPL e.g. by running `include("predictions.jl")`, `include("bootstrap.jl")` or `include("break_removal_montecarlo.jl")`.

Note that some time consuming tasks have been commented out and it might be necessary to walk through the code line by line and uncomment manually.
