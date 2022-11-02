



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
raw_chain_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20221028.csv", DataFrame)


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
n_var = 20
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
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Effect_2208_20221031.csv", effect_df)

effect_gdf = groupby(effect_df, :Var_Seq)
combine(effect_gdf, :x7 => minimum)

#=
para_df = DataFrame(
    β0 = value.(β0)
    , β_Admit = value.(β)[1]
    , β_home_distance = value.(β)[2]
    , β_Admit_Honor = value.(β)[3]
    , β_Diff_Major = value.(β)[4]
    , β_Gender = value.(β)[5]
    , β_inst_grant = value.(β)[6]
    , β_loan = value.(β)[7]
    , β_fed_efc = value.(β)[8]
    , β_Pell = value.(β)[9]
    , β_ASIAN = value.(β)[10]
    , β_BLACK = value.(β)[11]
    , β_HISPA = value.(β)[12]
    , β_WHITE = value.(β)[13]
    , β_Multi = value.(β)[14]
    , β_Postcard = value.(β)[15]
    , β_Pros_Event = value.(β)[16]
    , β_CampusTour = value.(β)[17]
    , β_DecisionDay = value.(β)[18]
    , β_Financing = value.(β)[19]
    , β_FinAid = value.(β)[20]
)
=#


```
whether variables are important
```

p_mtx = ones(n_var+1,n_Period)
p_intval_mtx = Array{Tuple{Float64 ,Float64}, 2}(undef, n_var+1, 8)
p_star_mtx = Array{String, 2}(undef, 21, 8)
for i in 1:(n_var+1)
    for j in 1:n_Period
        p_mtx[i,j] = pvalue(OneSampleZTest(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))))
        p_intval_mtx[i,j] = round.(confint(OneSampleZTest(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))));digits=3)
        if p_mtx[i,j] <0.01
            p_star_mtx[i,j] = "***"
        elseif p_mtx[i,j] <0.05
            p_star_mtx[i,j] = "**"
        elseif p_mtx[i,j] <0.1
            p_star_mtx[i,j] = "*"
        else p_star_mtx[i,j] = " "
        end
    end
end

p_mtx[4,3]
p_intval_mtx[4,3]
p_star_mtx[4,3]

p_df = DataFrame(hcat(1:(n_var+1),p_mtx), :auto)
p_intval_df = DataFrame(hcat(1:(n_var+1),p_intval_mtx), :auto)
p_star_df = DataFrame(hcat(1:(n_var+1),p_star_mtx), :auto)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/p_df_2208_20221101.csv", p_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/p_intval_df_2208_20221101.csv", p_intval_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/p_star_df_2208_20221101.csv", p_star_df)

```
whether a variable has time-varying effect
```

Comp_intval = Array{Tuple{Float64 ,Float64}, 2}(undef, n_var, Int(n_Period*(n_Period-1)/2))
for i in 1:n_var
    for j in 1:(n_Period-1)
        for k in (j+1):n_Period
            x = vec(Matrix(raw_chain[(i*n_Period+j):(i*n_Period+j),:]))
            y = vec(Matrix(raw_chain[(i*n_Period+k):(i*n_Period+k),:]))
            Comp_intval[i, Int(((j-1)*n_Period-j*(j-1)/2+(k-j)))] = round.(confint(UnequalVarianceZTest(x,y); level = 0.95, tail=:both); digits=3)
        end
    end
    
end

Comp_intval_df = DataFrame(hcat(1:n_var,Comp_intval), :auto)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Comp_intval_df_2208_20221101.csv", Comp_intval_df)

Need to get max difference for each variable