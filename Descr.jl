



```
Load packages
```

#using Random, Distributions

using Plots, StatsPlots, ColorSchemes

#using PlotlyJS

using CSV

using DataFrames

using FreqTables, StatsBase

#using HypothesisTests

using Statistics
```
Load raw data
```

raw_2208 = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2208_20221123.csv", DataFrame)
raw_2218 = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2218_20221123.csv", DataFrame)
raw_2228 = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2228_20221123.csv", DataFrame)

raw_chain_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20221123.csv", DataFrame)
raw_chain_2218 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2218_20221123.csv", DataFrame)
raw_chain_2228 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2228_20221123.csv", DataFrame)




```
Data and Variables
```
print(hcat(1:ncol(raw_2208), names(raw_2208)))
unique(raw_2208[raw_2208.Resid_Ind .==0 ,:].Student_ID)
unique(raw_2218[raw_2218.Resid_Ind .==0 ,:].Student_ID)
unique(raw_2228[raw_2228.Resid_Ind .==0 ,:].Student_ID)

# non-resident
18450+20619+21216

```
Table 1
```
freqtable(raw_2208[raw_2208.Resid_Ind .==0 ,:], :Period)

freqtable(raw_2208[raw_2208.Resid_Ind .==0 ,:], :Period, :Dpst_Ind)
freqtable(raw_2218[raw_2218.Resid_Ind .==0 ,:], :Period, :Dpst_Ind)
freqtable(raw_2228[raw_2228.Resid_Ind .==0 ,:], :Period, :Dpst_Ind)

```
Table 2
```
total_dpst = sum(raw_2208[raw_2208.Resid_Ind .==0 ,:].Dpst_Ind)+sum(raw_2218[raw_2218.Resid_Ind .==0 ,:].Dpst_Ind)+sum(raw_2228[raw_2228.Resid_Ind .==0 ,:].Dpst_Ind)
total_obs = length(unique(raw_2208[raw_2208.Resid_Ind .==0 ,:].Student_ID)) + length(unique(raw_2218[raw_2218.Resid_Ind .==0 ,:].Student_ID)) + length(unique(raw_2228[raw_2228.Resid_Ind .==0 ,:].Student_ID))
dpst_pct = total_dpst/total_obs

for i in 1:nrow(raw_2208)
    if ismissing(raw_2208.FinAid_Rate[i])
        raw_2208.FinAid_Rate[i] = 0
    end
end

for i in 1:nrow(raw_2218)
    if ismissing(raw_2218.FinAid_Rate[i])
        raw_2218.FinAid_Rate[i] = 0
    end
end

for i in 1:nrow(raw_2228)
    if ismissing(raw_2228.FinAid_Rate[i])
        raw_2228.FinAid_Rate[i] = 0
    end
end

var_list = [:FinAid_Rate,:fed_efc_rate, :Pell_Ind
                , :home_distance, :Gender_Ind
                , :Eth_BLACK_Ind, :Eth_ASIAN_Ind, :Eth_HISPA_Ind, :Eth_WHITE_Ind, :Eth_Multi_Ind
                , :Pros_Event_Ind, :Admit_Honor_Ind, :Diff_Major_Ind, :CampusTour_Ever_Ind, :DecisionDay_Ever_Ind, :Delay_Review_Ind]

gdf = groupby(raw_2208[raw_2208.Resid_Ind .==0 ,:], :Student_ID)
FinAid_1 = combine(gdf, var_list .=> mean; renamecols=false)
gdf = groupby(raw_2218[raw_2218.Resid_Ind .==0 ,:], :Student_ID)
FinAid_2 = combine(gdf, var_list .=> mean; renamecols=false)
gdf = groupby(raw_2228[raw_2228.Resid_Ind .==0 ,:], :Student_ID)
FinAid_3 = combine(gdf, var_list .=> mean; renamecols=false)

FinAid = vcat(FinAid_1,FinAid_2,FinAid_3)

mean(FinAid.FinAid_Rate)
std(FinAid.FinAid_Rate)

mean(FinAid.fed_efc_rate)
std(FinAid.fed_efc_rate)

println(round.(mean(Matrix(FinAid), dims=1); digits=3))
println(round.(sum(Matrix(FinAid), dims=1); digits=3))

std(FinAid.home_distance)
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
n_var = 16
n_σ = 2

raw_chain = raw_chain_2208
#raw_chain = raw_chain_2218
#raw_chain = raw_chain_2228

var_name = ["β_FinAid", "β_Pell", "β_efc", "β_home", "β_Gender", "β_ASIAN", "β_BLACK", "β_HISPA", "β_WHITE", "β_Multi",
           "β_Pros_Event", "β_Admit_Honor", "β_Diff_Major", "β_CampusTour", "β_DecisionDay","β_Delay_Review"]

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
#p_mtx[4,3]
#intval_mtx[4,3]
#p_star_mtx[4,3]

#p_df = DataFrame(hcat(1:(n_var+1),p_mtx), :auto)
#intval_df = DataFrame(hcat(1:(n_var+1),intval_mtx), :auto)
β_tupl_df = DataFrame(hcat(1:(n_var+1),β_tupl), :auto)
p_star_df = DataFrame(hcat(1:(n_var+1),p_star_mtx), :auto)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/tupl_df_2228_20221124.csv", β_tupl_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/p_star_df_2228_20221124.csv", p_star_df)

###CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/intval_df_2208_20221107.csv", intval_df)

findmax(p_star_mtx, dims=2)[1]


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
Comp_df = DataFrame(hcat(1:(n_var+1), findmax(p_star_mtx, dims=2)[1], findmax(Comp_mtx, dims=2)[1], Comp_mtx), :auto)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Comp_df_2208_20221123.csv", Comp_df)

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

sample_var = 4
plot_2 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Major Change", yaxis = "", label = "Period 8",legend=false)            

sample_var = 17
plot_3 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Campus Tour", yaxis = "", label = "Period 8")   

sample_var = 2
plot_4 = density(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7, label = "Period 1")
density!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+8):(sample_var*n_Period+8),:])), :, 1), alpha=0.7
            , xaxis = "Home Distance", yaxis = "", label = "Period 8")  

plot(plot_1,plot_2,plot_3,plot_4,layout = (2, 2), legend=false)


#####
# BoxPlot
#####

plot_mtx = ones(1000, n_Period)
for j in 1:n_Period
    plot_mtx[:,j] = vec(Matrix(raw_chain_2228[(j):(j),:]))
end
display(boxplot(["1" "2" "3" "4" "5" "6" "7" "8"], plot_mtx
            , legend=false, xaxis="Period", yaxis="Baseline"; palette = :grayC) )

raw_chain = raw_chain_2228


for i in 1:n_var
    plot_mtx = ones(1000, n_Period)
    for j in 1:n_Period
        plot_mtx[:,j] = vec(Matrix(raw_chain[(i*n_Period+j):(i*n_Period+j),:]))
    end
    display(boxplot(["1" "2" "3" "4" "5" "6" "7" "8"], plot_mtx
                , legend=false, xaxis="Period", yaxis=var_name[i]; palette = :grayC) )
end

```
Comparison among three falls
```

plot_mtx_1 = ones(n_Period, 1000)
plot_mtx_2 = ones(n_Period, 1000)
plot_mtx_3 = ones(n_Period, 1000)
for j in 1:n_Period
    plot_mtx_1[j,:] = vec(Matrix(raw_chain_2208[(j):(j),:]))
    plot_mtx_2[j,:] = vec(Matrix(raw_chain_2218[(j):(j),:]))
    plot_mtx_3[j,:] = vec(Matrix(raw_chain_2228[(j):(j),:]))
end
boxplot([0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8], vec(plot_mtx_1)
            , xaxis="Period", yaxis="Baseline", label="2020", bar_width = 0.2, seriescolor=:gray50)
boxplot!([1, 2, 3, 4, 5, 6, 7, 8], vec(plot_mtx_2)
            , xaxis="Period", yaxis="Baseline", label="2021", bar_width = 0.2, seriescolor=:gray75)            
display(boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis="Baseline", label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=:topleft, xticks=0:1:9)  )


for i in 1:n_var
    plot_mtx_1 = ones(n_Period, 1000)
    plot_mtx_2 = ones(n_Period, 1000)
    plot_mtx_3 = ones(n_Period, 1000)
    for j in 1:n_Period
        plot_mtx_1[j,:] = vec(Matrix(raw_chain_2208[(i*n_Period+j):(i*n_Period+j),:]))
        plot_mtx_2[j,:] = vec(Matrix(raw_chain_2218[(i*n_Period+j):(i*n_Period+j),:]))
        plot_mtx_3[j,:] = vec(Matrix(raw_chain_2228[(i*n_Period+j):(i*n_Period+j),:]))
    end

    boxplot([0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8], vec(plot_mtx_1)
            , xaxis="Period", yaxis=var_name[i], label="2020", bar_width = 0.2, seriescolor=:gray50) 
    boxplot!([1, 2, 3, 4, 5, 6, 7, 8], vec(plot_mtx_2)
            , xaxis="Period", yaxis=var_name[i], label="2021", bar_width = 0.2, seriescolor=:gray75) 
    display(boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis=var_name[i], label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=:outertopright, xticks=0:1:9) )

end

```
Figure 1
```

plot_mtx_1 = ones(n_Period, 1000)
plot_mtx_2 = ones(n_Period, 1000)
plot_mtx_3 = ones(n_Period, 1000)
for j in 1:n_Period
    plot_mtx_1[j,:] = vec(Matrix(raw_chain_2208[(j):(j),:]))
    plot_mtx_2[j,:] = vec(Matrix(raw_chain_2218[(j):(j),:]))
    plot_mtx_3[j,:] = vec(Matrix(raw_chain_2228[(j):(j),:]))
end
boxplot([0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8], vec(plot_mtx_1)
            , xaxis="Period", yaxis="Baseline", label="2020", bar_width = 0.2, seriescolor=:gray50)
boxplot!([1, 2, 3, 4, 5, 6, 7, 8], vec(plot_mtx_2)
            , xaxis="Period", yaxis="Baseline", label="2021", bar_width = 0.2, seriescolor=:gray75)            
plot_1 = boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis="Baseline", label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=:topleft, xticks=0:1:9, size=(600,400))  

var_num = 1
plot_mtx_1 = ones(n_Period, 1000)
plot_mtx_2 = ones(n_Period, 1000)
plot_mtx_3 = ones(n_Period, 1000)
for j in 1:n_Period
    plot_mtx_1[j,:] = vec(Matrix(raw_chain_2208[(var_num*n_Period+j):(var_num*n_Period+j),:]))
    plot_mtx_2[j,:] = vec(Matrix(raw_chain_2218[(var_num*n_Period+j):(var_num*n_Period+j),:]))
    plot_mtx_3[j,:] = vec(Matrix(raw_chain_2228[(var_num*n_Period+j):(var_num*n_Period+j),:]))
end
boxplot([0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8], vec(plot_mtx_1)
            , xaxis="Period", yaxis="Financial Aid", label="2020", bar_width = 0.2, seriescolor=:gray50)
boxplot!([1, 2, 3, 4, 5, 6, 7, 8], vec(plot_mtx_2)
            , xaxis="Period", yaxis="Financial Aid", label="2021", bar_width = 0.2, seriescolor=:gray75)            
plot_2 = boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis="Financial Aid", label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=false, xticks=0:1:9, size=(1200,400))  

var_num = 5
plot_mtx_1 = ones(n_Period, 1000)
plot_mtx_2 = ones(n_Period, 1000)
plot_mtx_3 = ones(n_Period, 1000)
for j in 1:n_Period
    plot_mtx_1[j,:] = vec(Matrix(raw_chain_2208[(var_num*n_Period+j):(var_num*n_Period+j),:]))
    plot_mtx_2[j,:] = vec(Matrix(raw_chain_2218[(var_num*n_Period+j):(var_num*n_Period+j),:]))
    plot_mtx_3[j,:] = vec(Matrix(raw_chain_2228[(var_num*n_Period+j):(var_num*n_Period+j),:]))
end
boxplot([0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8], vec(plot_mtx_1)
            , xaxis="Period", yaxis="Gender", label="2020", bar_width = 0.2, seriescolor=:gray50)
boxplot!([1, 2, 3, 4, 5, 6, 7, 8], vec(plot_mtx_2)
            , xaxis="Period", yaxis="Gender", label="2021", bar_width = 0.2, seriescolor=:gray75)            
plot_3 = boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis="Gender", label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=false, xticks=0:1:9, size=(600,400)) 

var_num = 15
plot_mtx_1 = ones(n_Period, 1000)
plot_mtx_2 = ones(n_Period, 1000)
plot_mtx_3 = ones(n_Period, 1000)
for j in 1:n_Period
    plot_mtx_1[j,:] = vec(Matrix(raw_chain_2208[(var_num*n_Period+j):(var_num*n_Period+j),:]))
    plot_mtx_2[j,:] = vec(Matrix(raw_chain_2218[(var_num*n_Period+j):(var_num*n_Period+j),:]))
    plot_mtx_3[j,:] = vec(Matrix(raw_chain_2228[(var_num*n_Period+j):(var_num*n_Period+j),:]))
end
boxplot([0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8], vec(plot_mtx_1)
            , xaxis="Period", yaxis="Decision Day", label="2020", bar_width = 0.2, seriescolor=:gray50)
boxplot!([1, 2, 3, 4, 5, 6, 7, 8], vec(plot_mtx_2)
            , xaxis="Period", yaxis="Decision Day", label="2021", bar_width = 0.2, seriescolor=:gray75)            
plot_4 = boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis="Decision Day", label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=false, xticks=0:1:9, size=(600,400)) 

plot(plot_1,plot_2,plot_3,plot_4,layout = (2, 2))