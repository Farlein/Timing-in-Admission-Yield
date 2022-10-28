



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

raw_chain_2228 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2228_20220819.csv", DataFrame)
raw_chain_2218 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2218_20220819.csv", DataFrame)
raw_chain_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20220819.csv", DataFrame)


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

```
Results
```


n_Period = 8
n_var = 19
n_σ = 2

raw_chain = raw_chain_2208

#=
raw_chain[(1*n_Period+1):(1*n_Period+n_Period),:]
chain_mtx = transpose(Matrix(raw_chain[(1*n_Period+1):(1*n_Period+n_Period),:]))
percentile(chain_mtx[:,1], 95)
percentile(chain_mtx[:,1], 5)
percentile(chain_mtx[:,1], 50)

percentile(chain_mtx[:,7], 95)
percentile(chain_mtx[:,7], 5)
percentile(chain_mtx[:,7], 50)

effect_1 = [percentile(chain_mtx[:,1], 5), mean(chain_mtx[:,1]), percentile(chain_mtx[:,1], 50), percentile(chain_mtx[:,1], 95)]
=#

var_mtx = transpose(Matrix(raw_chain[(1*n_Period+1):(1*n_Period+n_Period),:]))
effect_summ = [1, 1, percentile(var_mtx[:,1], 5), mean(var_mtx[:,1]), percentile(var_mtx[:,1], 50), percentile(var_mtx[:,1], 95)]

for j in 2:n_Period
    effect_temp = [1, j, percentile(var_mtx[:,j], 5), mean(var_mtx[:,j]), percentile(var_mtx[:,j], 50), percentile(var_mtx[:,j], 95)]
    effect_summ = hcat(effect_summ, effect_temp)
end

for i in 2:n_var
    var_mtx = transpose(Matrix(raw_chain[(i*n_Period+1):(i*n_Period+n_Period),:]))
    for j in 1:n_Period
        effect_temp = [i, j, percentile(var_mtx[:,j], 5), mean(var_mtx[:,j]), percentile(var_mtx[:,j], 50), percentile(var_mtx[:,j], 95)]
        effect_summ = hcat(effect_summ, effect_temp)
    end
end

effect = transpose(effect_summ)
effect_df = DataFrame(hcat( effect, (effect[:,3] .<0) .* (effect[:,6] .>0)) , :auto )
rename!(effect_df, :x1 => :Var_Seq)
rename!(effect_df, :x2 => :Period)
rename!(effect_df, :x3 => :Percentile_5)
rename!(effect_df, :x4 => :Mean)
rename!(effect_df, :x5 => :Median)
rename!(effect_df, :x6 => :Percentile_95)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Effect_2208_20221018.csv", effect_df)

effect_gdf = groupby(effect_df, :Var_Seq)
combine(effect_gdf, :x7 => minimum)

