


```
Load packages
```

#using Random, Distributions

using Plots, StatsPlots

#using PlotlyJS

using CSV

using DataFrames

using FreqTables, StatsBase

using HypothesisTests

using Random, Distributions

```
Load raw data

```

###
raw = CSV.read("C:\\DATA\\FSAN\\Adm_Yield\\for_julia_all_20220315.csv", DataFrame)
names(raw)

raw.headcount = ones(size(raw)[1])

### Number of admits and deposits
combine(groupby(raw[(raw.week .== 12) .& (raw.Resid_Ind .== 0), :], :apply_term), :headcount => sum)
combine(groupby(raw[(raw.week .== 13) .& (raw.Resid_Ind .== 0), :], :apply_term), :headcount => sum)
combine(groupby(raw[(raw.week .== 22) .& (raw.Resid_Ind .== 0), :], :apply_term), :headcount => sum)

combine(groupby(raw[(raw.week .== 22) .& (raw.Resid_Ind .== 0), :], :apply_term), :dpst_flag => sum)

### Yield on May 1
yield_May = combine(groupby(raw[(raw.week .== 22) .& (raw.Resid_Ind .== 0), :], :apply_term), :dpst_flag => mean)

plot_x = string.("20", string.(SubString.(string.(yield_May.apply_term[5:8]),2,3)))
plot_y = round.(yield_May.dpst_flag_mean[5:8],digits=3)

#plotlyjs()
#Plots.PlotlyJSBackend()



#gr()

plot(
    bar(plot_x, plot_y, series_annotations = plot_y
    , label = ""
    , xaxis = "Year", yaxis = "Yield", guidefontsize = 20
    , tickfontsize = 16
    , title = "Admission Yield on May 1", titlefontsize = 20
    )
     )

plot(
    plot_x, plot_y, st= :bar
    , label = ""
    , xaxis = "Year", yaxis = "Yield", guidefontsize = 20
    , tickfontsize = 16
    , title = "Admission Yield on May 1", titlefontsize = 20
     )

plot(
    plot_x, plot_y;
    st= :bar, txt = plot_y
    , label = ""
    , xaxis = "Year", yaxis = "Yield", guidefontsize = 20
    , tickfontsize = 16
    , title = "Admission Yield on May 1", titlefontsize = 20
     )


### Yield trend


### Yield trend
trd_2198 = combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2198), :], :week), :dpst_flag => mean)
trd_2188 = combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2188), :], :week), :dpst_flag => mean)
trd_2178 = combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2178), :], :week), :dpst_flag => mean)
trd_2168 = combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2168), :], :week), :dpst_flag => mean)

scatter(1:11, trd_2198.dpst_flag_mean, label = "2019")
scatter!(1:11, trd_2188.dpst_flag_mean, label = "2018")
scatter!(1:11, trd_2178.dpst_flag_mean, label = "2017")
scatter!(1:11, trd_2168.dpst_flag_mean, label = "2016"
            , legend = :topleft
            , xaxis = "Week", yaxis = "Yield"
            , guidefontsize = 20
            , tickfontsize = 16)

### Yield trend
combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2198), :], :week), :dpst_prop => mean)
combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2188), :], :week), :dpst_prop => mean)
combine(groupby(raw[(raw.Resid_Ind .== 0) .& (raw.week .>= 12) .& (raw.apply_term .== 2178), :], :week), :dpst_prop => mean)

histogram(raw_target.fed_efc_rate)
histogram(raw_source.fed_efc_rate)

mean(raw_target.fed_efc_rate)
mean(raw_source.fed_efc_rate)

histogram(raw_target.inst_grant_rate)
histogram(raw_source.inst_grant_rate)

mean(raw_target.inst_grant_rate)
mean(raw_source.inst_grant_rate)

mean(raw_target.Feeder_HS)
mean(raw_source.Feeder_HS)

mean(raw_target.Eth_WHITE_Ind)
mean(raw_source.Eth_WHITE_Ind)

mean(raw_target.Pell_Ind)
mean(raw_source.Pell_Ind)

mean(raw_target.Fst_Gen_Ind)
mean(raw_source.Fst_Gen_Ind)

mean(raw_target.Eth_BLACK_Ind)
mean(raw_source.Eth_BLACK_Ind)


```

Compare parameters

```

```
Load raw parameters
```

###
raw_para_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_20220719.csv", DataFrame)
raw_para_2218 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2218_20220719.csv", DataFrame)
raw_para_2228 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2228_20220719.csv", DataFrame)

para_2208 = raw_para_2208[2:9, :]
para_2208.term .= "Fall 2020"
para_2208.period = 1:8

para_2218 = raw_para_2218[2:9, :]
para_2218.term .= "Fall 2021"
para_2218.period = 1:8

para_2228 = raw_para_2228[2:9, :]
para_2228.term .= "Fall 2022"
para_2228.period = 1:8

#para_union = vcat(para_2208, para_2218, para_2228)

para_list = names(para_2208)

plot(para_2208.β_Delay_Review, label = "Fall 2020")
plot!(para_2218.β_Delay_Review, label = "Fall 2021")
plot!(para_2228.β_Delay_Review, label = "Fall 2022"
        , xaxis = ("Period")
        , yaxis = names(para_2228)[21]
        , legend=:topleft)

function cc_plot_para(para_num, para_name)
    plot(para_2208[:, para_num] , label = "Fall 2020")
    plot!(para_2218[:, para_num] , label = "Fall 2021")
    plot!(para_2228[:, para_num] , label = "Fall 2022"
            , xaxis = ("Period")
            , yaxis = para_name
            , legend=:topleft)
end

cc_plot_para(1,"β0")
cc_plot_para(2,"β_Admit")
cc_plot_para(3,"β_hs_gpa")
cc_plot_para(4,"β_Admit_Honor")
cc_plot_para(5,"β_Diff_Major")
cc_plot_para(6,"β_Gender")
cc_plot_para(7,"β_inst_grant")
cc_plot_para(8,"β_loan")
cc_plot_para(9,"β_fed_efc")
cc_plot_para(10,"β_Pell")
cc_plot_para(11,"β_ASIAN")
cc_plot_para(12,"β_BLACK")
cc_plot_para(13,"β_HISPA")
cc_plot_para(14,"β_WHITE")
cc_plot_para(15,"β_Slate_Time")
cc_plot_para(16,"β_Slate_Absent")
cc_plot_para(17,"β_Postcard")
cc_plot_para(18,"β_Blue_Gold")
cc_plot_para(19,"β_Early_Event")
cc_plot_para(20,"β_Event")
cc_plot_para(21,"β_Delay_Review")




cc_plot_para(21,para_list[21])


```

Draw Trend Plots

```

```
Load raw parameters
```
raw_chain_2228 = CSV.read("I:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2228_20220809_2.csv", DataFrame)


plot(Vector(raw_chain_2228[15,:]))
plot(Vector(raw_chain_2228[15,:]), label = "β_Diff_Major_5")



var_name = ["β_Diff_Major_1", "β_Diff_Major_5"]
plot(Vector(raw_chain_2228[15,:]), label = var_name[2])

var_name = []
for i in 1:10
    var_temp = "β0_"*string(i)
    if i == 1
        var_name = var_temp
    else var_name = vcat(var_name,var_temp)
    end
end

for i in 1:10
    var_temp = "β_Diff_Major_"*string(i)
    var_name = vcat(var_name,var_temp)
end

for i in 1:10
    var_temp = "β_inst_grant_"*string(i)
    var_name = vcat(var_name,var_temp)
end

for i in 3:30
    display(plot(Vector(raw_chain_2228[i,:]), label = var_name[i]))
end



```

Draw Trend Plots

```

```
Load raw parameters
```
raw_chain_2228 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2228_20220819.csv", DataFrame)
raw_chain_2218 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2218_20220819.csv", DataFrame)
raw_chain_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20220819.csv", DataFrame)

n_Period = 8
n_var = 19
n_σ = 2

raw_chain = raw_chain_2208
μ_chain_df = mean(Matrix(raw_chain), dims = 2)

n_β = (n_var+1)*n_Period
β0_fit = μ_chain_df[1:n_Period]
β_fit = reshape(μ_chain_df[(n_Period+1):n_β], (n_Period,n_var))
display(β_fit)
β_h_fit = μ_chain_df[(n_β+1):(n_β+n_var)]
σ_h_fit = exp(μ_chain_df[n_β+n_var+1])
σ_l_fit = exp(μ_chain_df[n_β+n_var+2])

### Trend Plot

#var_name = ["β_Admit", "β_Admit_Honor", "β_Diff_Major", "β_Gender", "β_inst_grant", "β_loan", "β_fed_efc", "β_Pell", 
#            "β_ASIAN", "β_BLACK", "β_HISPA", "β_WHITE", "β_Slate_Time", "β_Slate_Absent", "β_Postcard", "β_Blue_Gold", "β_Early_Event", "β_Event", "β_Delay_Review"]

var_name = ["β_Admit", "β_hs_gpa", "β_Admit_Honor", "β_Diff_Major", "β_Gender", "β_inst_grant", "β_loan", "β_fed_efc", "β_Pell", 
            "β_ASIAN", "β_BLACK", "β_HISPA", "β_WHITE", "β_Multi", "β_Major_Finder", "β_SFS", "β_Blue_Gold", "β_Early_Event", "β_Delay_Review"]

print(hcat(Vector(1:n_var), var_name))

display(plot(reshape(transpose(Array(raw_chain[1:n_Period,:])), :, 1), label = "", yaxis = ("β0"), legend=:topleft))

#plot(reshape(transpose(Array(raw_chain[(1*n_Period+2):(1*n_Period+9),:])), :, 1), label = "", yaxis = var_name[1], legend=:topleft)
#plot(Vector(raw_chain[15,:]), label = "β_Diff_Major_5")

for i in 1:n_var
   display(plot(reshape(transpose(Array(raw_chain[(i*n_Period+1):(i*n_Period+n_Period),:])), :, 1), label = "", yaxis = var_name[i], legend=:topleft)) 
end

#=
for i in 1:n_var
    display(plot(reshape(transpose(Array(raw_chain[(n_var*n_Period+n_Period+i),:])), :, 1), label = "", yaxis = var_name[i], legend=:topleft)) 
end
=#

σ_name = ["σ_h", "σ_l"]

for i in 1:n_σ
    display(plot(reshape(transpose(Array(raw_chain[(n_var*n_Period+n_Period+n_var+i),:])), :, 1), label = "", yaxis = σ_name[i], legend=:topleft)) 
end



```

Draw Interaction Plots

```

var_name = ["β_Admit", "β_home_distance", "β_Admit_Honor", "β_Diff_Major", "β_Gender", "β_inst_grant", "β_loan", "β_fed_efc", "β_Pell", 
            "β_ASIAN", "β_BLACK", "β_HISPA", "β_WHITE", "β_Multi", "β_Postcard", "β_Pros_Event", "β_CampusTour", "β_DecisionDay", "β_Financing","β_FinAid"]

plot(reshape(transpose(Array(raw_chain[(6*n_Period+1):(6*n_Period+n_Period),:])), :, 1), reshape(transpose(Array(raw_chain[(7*n_Period+1):(7*n_Period+n_Period),:])), :, 1))

plot(reshape(transpose(Array(raw_chain[(6*n_Period+n_Period):(6*n_Period+n_Period),:])), :, 1), reshape(transpose(Array(raw_chain[(7*n_Period+n_Period):(7*n_Period+n_Period),:])), :, 1))

for i in 1:n_var
    for j in i:n_var
        display(plot(reshape(transpose(Array(raw_chain[(i*n_Period+1):(i*n_Period+n_Period),:])), :, 1), reshape(transpose(Array(raw_chain[(j*n_Period+1):(j*n_Period+n_Period),:])), :, 1),
                        label = "", xaxis = var_name[i], yaxis = var_name[j]))
    end
end

for i in 8:8
    for j in 16:16
        for k in 1:n_Period
            display(plot(reshape(transpose(Array(raw_chain[(i*n_Period+k):(i*n_Period+k),:])), :, 1), reshape(transpose(Array(raw_chain[(j*n_Period+k):(j*n_Period+k),:])), :, 1),
                        label = k, xaxis = var_name[i], yaxis = var_name[j]))
        end
    end
end



```

Draw Trend Plots
Three years together

```

```
Load raw chains
```
raw_chain_2208 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20221123.csv", DataFrame)
raw_chain_2218 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2218_20221123.csv", DataFrame)
raw_chain_2228 = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2228_20221123.csv", DataFrame)



n_Period = 8
n_var = 16
n_σ = 2

plot(reshape(transpose(Array(raw_chain_2208[1:n_Period,:])), :, 1), label = "Fall 2020", yaxis = ("β0"), legend=:topleft, alpha=0.2) 
plot!(reshape(transpose(Array(raw_chain_2218[1:n_Period,:])), :, 1), label = "Fall 2021", alpha=0.2)
display(plot!(reshape(transpose(Array(raw_chain_2228[1:n_Period,:])), :, 1), label = "Fall 2022", alpha=0.2))

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

var_name = ["β_Admit", "β_home_distance", "β_Admit_Honor", "β_Diff_Major", "β_Gender", "β_inst_grant", "β_loan", "β_fed_efc", "β_Pell", 
            "β_ASIAN", "β_BLACK", "β_HISPA", "β_WHITE", "β_Multi", "β_Postcard", "β_Pros_Event", "β_CampusTour", "β_DecisionDay", "β_Delay_Review"]

for i in 1:n_var
    plot(reshape(transpose(Array(raw_chain_2208[(i*n_Period+1):(i*n_Period+n_Period),:])), :, 1), label = "Fall 2020", alpha=0.2)
    plot!(reshape(transpose(Array(raw_chain_2218[(i*n_Period+1):(i*n_Period+n_Period),:])), :, 1), label = "Fall 2021", alpha=0.2)
    display(plot!(reshape(transpose(Array(raw_chain_2228[(i*n_Period+1):(i*n_Period+n_Period),:])), :, 1), label = "Fall 2022", yaxis = var_name[i], legend=:topleft, alpha=0.2)) 
 end

##############
sample_var = 6
quantile(vec(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), 0.025)
quantile(vec(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), 0.975)
OneSampleZTest(vec(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), 0.0)

histogram(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1), alpha=0.7)
histogram!(reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+5):(sample_var*n_Period+5),:])), :, 1), alpha=0.7)

x = reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])), :, 1)
y = reshape(transpose(Array(raw_chain_2208[(sample_var*n_Period+5):(sample_var*n_Period+5),:])), :, 1)
x = vec(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])))
y = vec(transpose(Array(raw_chain_2208[(sample_var*n_Period+5):(sample_var*n_Period+5),:])))

z = UnequalVarianceZTest(x,y)

x = vec(transpose(Array(raw_chain_2208[(sample_var*n_Period+1):(sample_var*n_Period+1),:])))
y = vec(transpose(Array(raw_chain_2208[(sample_var*n_Period+7):(sample_var*n_Period+7),:])))
UnequalVarianceZTest(x,y)
confint(UnequalVarianceZTest(x,y); level = 0.95, tail=:both)

transpose(Matrix(raw_chain[(i*n_Period+1):(i*n_Period+n_Period),:]))

OneSampleZTest(vec(Matrix(raw_chain[((4-1)*n_Period+3):((4-1)*n_Period+3),:])))
confint(OneSampleZTest(vec(Matrix(raw_chain[((4-1)*n_Period+3):((4-1)*n_Period+3),:]))); level = 0.90, tail=:both)

boxplot(x)


###

plot_mtx = ones(1000, n_Period)
for j in 1:n_Period
    plot_mtx[:,j] = vec(Matrix(raw_chain_2228[(j):(j),:]))
end
display(boxplot(["1" "2" "3" "4" "5" "6" "7" "8"], plot_mtx
            , legend=false, xaxis="Period", yaxis="Baseline"; palette = :grayC) )

d = Normal()
test_1 = rand(d, 100)
test_2 = rand(d, 100)
test_3 = rand(d, 100)

test_df = DataFrame(y=vcat(test_1,test_2,test_3), z=repeat(["2020"; "2021"; "2022"], inner=100))
gr()
@df test_df boxplot(
    :y,
    #group=:z
)
@df test_df boxplot(
    :z,
    :y
)

test_df_1 = DataFrame(x=repeat(1:8,inner=100), y=repeat(test_1, 8))
test_df_2 = DataFrame(x=repeat([1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1],inner=100), y=repeat(test_2, 8))
@df test_df_1 boxplot(
    :x, :y,
    bar_width = 0.1
)
@df test_df_2 boxplot!(
    :x, :y,
    bar_width = 0.1
)

boxplot([0.9, 1.9, 2.9, 3.9, 4.9, 5.9, 6.9, 7.9], repeat(test_1,8), xlims = (0,9), bar_width = 0.1)
boxplot!([1, 2, 3, 4, 5, 6, 7, 8], repeat(test_1,8), xlims = (0,9), bar_width = 0.1)
boxplot!([1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1], repeat(test_1,8), xlims = (0,9), bar_width = 0.1)
boxplot!(["1","2"],)

boxplot(["0.5"], test_1)
boxplot!(["1"], test_2)
boxplot!(["1.5"], test_3)

boxplot(["1.1","2.1"], repeat(test_1,2)
        , xlims = (0,9), bar_width = 0.1)
##############

i_test = 10
"β0_"*string(i_test)


plot(reshape(transpose(Array(raw_chain_2218[1:n_Period,:])), :, 1), label = "Fall 2021")

p_star_mtx = Array{String, 2}(undef, 21, 8)


d = Normal()
quantile(d, 0.5)
quantile(d, 0.95)
pdf(d, 0)
cdf(d, 0)
d_rdm = rand(d, 1000)

quantile(d_rdm, 0.025)
quantile(d_rdm, 0.975)
OneSampleZTest(d_rdm, 0.0)

percentile(d_rdm,2.5)


```
Comparison among three falls
```

β_mean_1 = ones(n_var+1,n_Period)
β_mean_2 = ones(n_var+1,n_Period)
β_mean_3 = ones(n_var+1,n_Period)

for i in 1:(n_var+1)
    for j in 1:n_Period
        β_mean_1[i,j] = round(mean(vec(Matrix(raw_chain_2208[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))); digits=2)
        β_mean_2[i,j] = round(mean(vec(Matrix(raw_chain_2218[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))); digits=2)
        β_mean_3[i,j] = round(mean(vec(Matrix(raw_chain_2228[((i-1)*n_Period+j):((i-1)*n_Period+j),:]))); digits=2)
    end
end

plot(1:n_Period, β_mean_1[1,:], label="Fall 2020")
plot!(1:n_Period, β_mean_2[1,:], label="Fall 2021")
plot!(1:n_Period, β_mean_3[1,:], label="Fall 2022"
        , legend = :topleft, xaxis="Period", yaxis = "Baseline Force")




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
plot_0 = boxplot!([1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2], vec(plot_mtx_3)
            , xaxis="Period", yaxis="Baseline", label="2022", bar_width = 0.2, seriescolor=:gray100
            , legend=:outertopright, xticks=0:1:9, size=(600,400))          

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

l = @layout [a{0.3w} grid(2,2)]
plot(plot_0, plot_1,plot_2,plot_3,plot_4, layout = l)

plot(plot_1,plot_2,plot_3,plot_4, layout = 4)