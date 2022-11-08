



```
Load packages
```

#using Random, Distributions

using Plots, StatsPlots, ColorSchemes

#using PlotlyJS

using CSV

using DataFrames

using FreqTables, StatsBase

using HypothesisTests

```
Load raw data
```

raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2208_20220816.csv", DataFrame)
raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2218_20220816.csv", DataFrame)
raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2228_20220816.csv", DataFrame)

raw_chain_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20221102.csv", DataFrame)
raw_chain_2218 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2218_20221102.csv", DataFrame)
raw_chain_2228 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2228_20221102.csv", DataFrame)




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


```
whether variables are important
```

β_mean = ones(n_var+1,n_Period)
β_low = ones(n_var+1,n_Period)
β_high = ones(n_var+1,n_Period)
β_tupl = Array{Tuple{Float64, Float64, Float64}, 2}(undef, n_var+1, n_Period)
#p_mtx = ones(n_var+1,n_Period)
#intval_mtx = Array{Tuple{Float64 ,Float64}, 2}(undef, n_var+1, 8)
p_star_mtx = Array{String, 2}(undef, n_var+1, n_Period)
for i in 1:(n_var+1)
    for j in 1:n_Period
        β_mean[i,j] = round(mean(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))); digits=2)
        β_low[i,j] = round(percentile(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:])), 2.5); digits=2)
        β_high[i,j] = round(percentile(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:])), 97.5); digits=2)
        β_tupl[i,j] = (β_low[i,j], β_mean[i,j], β_high[i,j])
 #       p_mtx[i,j] = pvalue(OneSampleZTest(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))))
 #       intval_mtx[i,j] = round.(confint(OneSampleZTest(vec(Matrix(raw_chain[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))));digits=3)
        if (β_low[i,j] <0.0) && (β_high[i,j] > 0.0)
            p_star_mtx[i,j] = " "
        else p_star_mtx[i,j] = "*"
        end
    end
end

β_mean[1,:]
β_mean[10,:]
p_mtx[4,3]
#intval_mtx[4,3]
#p_star_mtx[4,3]

#p_df = DataFrame(hcat(1:(n_var+1),p_mtx), :auto)
#intval_df = DataFrame(hcat(1:(n_var+1),intval_mtx), :auto)
β_tupl_df = DataFrame(hcat(1:(n_var+1),β_tupl), :auto)
p_star_df = DataFrame(hcat(1:(n_var+1),p_star_mtx), :auto)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/tupl_df_2208_20221107.csv", β_tupl_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/p_star_df_2208_20221107.csv", p_star_df)

###CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/intval_df_2208_20221107.csv", intval_df)

### Baseline
β_test = -5.93-(-7.47)
θ_test = 0.05
1-exp(exp(β_test)log(1-θ_test))

```
whether a variable has time-varying effect
```

Comp_mtx = Array{String, 2}(undef, n_var+1, Int(n_Period*(n_Period-1)/2)) 
#Comp_intval = Array{Tuple{Float64 ,Float64}, 2}(undef, n_var, Int(n_Period*(n_Period-1)/2))

for i in 1:(n_var+1)
    for j in 1:(n_Period-1)
        for k in (j+1):n_Period
            if (β_high[i,j] < β_low[i,k]) || (β_high[i,k] < β_low[i,j])
                Comp_mtx[i, Int(((j-1)*n_Period-j*(j-1)/2+(k-j)))] = "**"
            else Comp_mtx[i, Int(((j-1)*n_Period-j*(j-1)/2+(k-j)))] = " "
            end
#            x = vec(Matrix(raw_chain[(i*n_Period+j):(i*n_Period+j),:]))
#            y = β_low[i,k]
#            Comp_intval[i, Int(((j-1)*n_Period-j*(j-1)/2+(k-j)))] = round.(confint(UnequalVarianceZTest(x,y); level = 0.95, tail=:both); digits=3)
#            z_mtx[i, Int(((j-1)*n_Period-j*(j-1)/2+(k-j)))] = round(abs(UnequalVarianceZTest(x,y).z); digits=1)
        end
    end
    
end

findmax(Comp_mtx, dims=2)[1]
Comp_df = DataFrame(hcat(1:(n_var+1),findmax(Comp_mtx, dims=2)[1], Comp_mtx), :auto)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Comp_df_2208_20221107.csv", Comp_df)

###CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Comp_intval_df_2208_20221101.csv", Comp_intval_df)
###CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/z_df_2208_20221101.csv", z_df)

```
Examples of time-varying effects
```
sample_var = 9
plot_1 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Pell", yaxis = "", label = "Period 8")

sample_var = 5
plot_2 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Gender", yaxis = "", label = "Period 8",legend=false)            

sample_var = 16
plot_3 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Early Events", yaxis = "", label = "Period 8")   

sample_var = 18
plot_4 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Decision Day Event", yaxis = "", label = "Period 8")  

plot(plot_1,plot_2,plot_3,plot_4,layout = (2, 2), legend=false)

```
Examples of time-independent effects
```
sample_var = 7
plot_1 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Loan", yaxis = "", label = "Period 8")

sample_var = 5
plot_2 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Gender", yaxis = "", label = "Period 8",legend=false)            

sample_var = 16
plot_3 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Early Events", yaxis = "", label = "Period 8")   

sample_var = 18
plot_4 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Decision Day Event", yaxis = "", label = "Period 8")  

plot(plot_1,plot_2,plot_3,plot_4,layout = (2, 2), legend=false)

#####
# BoxPlot
#####

var_name = ["β_Admit", "β_home_distance", "β_Admit_Honor", "β_Diff_Major", "β_Gender", "β_inst_grant", "β_loan", "β_fed_efc", "β_Pell", 
            "β_ASIAN", "β_BLACK", "β_HISPA", "β_WHITE", "β_Multi", "β_Postcard", "β_Pros_Event", "β_CampusTour", "β_DecisionDay", "β_Delay_Review"]


for i in 1:n_var
    plot_mtx = ones(1000, n_Period)
    for j in 1:n_Period
        plot_mtx[:,j] = vec(Matrix(raw_chain[(i*n_Period+j):(i*n_Period+j),:]))
    end
    display(boxplot(["1" "2" "3" "4" "5" "6" "7" "8"], plot_mtx
                , legend=false, xaxis="Period", yaxis=var_name[i]; palette = :grayC) )
end


