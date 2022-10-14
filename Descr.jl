



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

raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2208_20220816.csv", DataFrame)
raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2218_20220816.csv", DataFrame)
raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2228_20220816.csv", DataFrame)


hcat(1:ncol(raw), names(raw))


```
Table 1
```
freqtable(raw, :Period)

freqtable(raw, :Period, :Dpst_Ind)

```
Table 3
```

mean(raw.inst_grant_rate)
std(raw.inst_grant_rate)
minimum(raw.inst_grant_rate)
maximum(raw.inst_grant_rate)