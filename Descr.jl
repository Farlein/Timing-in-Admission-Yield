



```
Load packages
```

#using Random, Distributions

using Plots, StatsPlots

#using PlotlyJS

using CSV

using DataFrames

using FreqTables, StatsBase


```
Load raw data
```

raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2208_20220816.csv", DataFrame)
raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2218_20220816.csv", DataFrame)
raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2228_20220816.csv", DataFrame)

hcat(1:ncol(raw), names(raw))

freqtable(raw, :Period)

freqtable(raw, :Period, :Dpst_Ind)