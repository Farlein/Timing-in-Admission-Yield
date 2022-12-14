

```
Load packages
```
# using Turing, Random 

# using DifferentialEquations, DiffEqBayes

# using MCMCChains, Distributions

using Random, Distributions

using Plots, StatsPlots

using CSV

using DataFrames

using FreqTables, StatsBase

#using JuMP, Ipopt

using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics
using MCMCDiagnostics
using Parameters, Statistics
#import ForwardDiff  
#import ReverseDiff
import Zygote



```
Load raw data

```

###
raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2208_20221123.csv", DataFrame)
#raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2218_20221123.csv", DataFrame)
#raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2228_20221123.csv", DataFrame)

#raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2228_20220615.csv", DataFrame)


hcat(1:ncol(raw), names(raw))
#print(hcat(1:ncol(raw), names(raw)))
#display(hcat(1:ncol(raw), names(raw)))
unique(raw.Apply_Term)

#raw_source = raw[ (raw.Apply_Term .== "Fall 2020") .& (raw.Resid_Ind .== 0), :]
#raw_source = raw[ (raw.Apply_Term .== "Fall 2021") .& (raw.Resid_Ind .== 0), :]
raw_source = raw[ (raw.Resid_Ind .== 0), :]

dpst_source = combine(groupby(raw_source, :Period), :Dpst_Ind => sum)
bar(dpst_source.Period, dpst_source.Dpst_Ind_sum, label = "")


#=
Periods:

1: Feb 1 to Feb 28/29
2: Mar 1 to Mar 15
3: Mar 15 to Mar 31
4: Apr 1 to Apr 7
5: Apr 8 to Apr 14
6: Apr 15 to Apr 21
7: Apr 22 to Apr 28
8: Apr 29 to May 1

=#


```
Data Management
```

#data_fit = raw_source[( 2 .<= raw_source.Period .<= 9 ), :]
data_fit = raw_source
unique(data_fit.Period)

# numeric variable

n_Period = 8
n_fit = size(data_fit)[1]

#=
data_fit.Delay_Review_Ind = zeros(nrow(data_fit))
for i in eachindex(data_fit.Delay_Review_Decision)
    if data_fit.Delay_Review_Decision[i] >= 2
        data_fit.Delay_Review_Ind[i] = 1
        else data_fit.Delay_Review_Ind[i] = 0
    end
end
=#

for i in 1:nrow(data_fit)
    if ismissing(data_fit.FinAid_Rate[i])
        data_fit.FinAid_Rate[i] = 0
    end
end

mean(data_fit.FinAid_Rate)
gdf = groupby(data_fit, :Period)
combine(gdf, [:FinAid_Rate] .=> mean; renamecols=false)

data_fit.Non_IG_Rate = data_fit.FinAid_Rate .- data_fit.inst_grant_rate


maximum(data_fit.Delay_Review_Ind)
minimum(data_fit.Delay_Review_Ind)
median(data_fit.Delay_Review_Ind)
mean(data_fit.Delay_Review_Ind)

function cc_standardize(var_seq, var_name)

    for i in 1:n_Period
        var_temp = (data_fit[ data_fit.Period .==i , var_seq] .- mean(data_fit[ data_fit.Period .==i , var_seq]) ) ./ std(data_fit[ data_fit.Period .==i , var_seq])
        join_key = data_fit[ data_fit.Period .==i , [2,3]]
        if i == 1
            var_name = hcat(join_key, var_temp)
            else var_name = vcat(var_name, hcat(join_key, var_temp))
        end
    end
    
    var_name
end

print(hcat(1:43, names(data_fit)))

#=
hs_gpa_std = []
hs_gpa_std = cc_standardize(10, hs_gpa_std)
rename!(hs_gpa_std,:x1 => :hs_gpa_std)
data_fit = innerjoin(data_fit, hs_gpa_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :hs_gpa_std => mean)
combine(groupby(data_fit, :Period), :hs_gpa_std => std)
=#

home_distance_std = []
home_distance_std = cc_standardize(11, home_distance_std)
rename!(home_distance_std,:x1 => :home_distance_std)
data_fit = innerjoin(data_fit, home_distance_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :home_distance_std => mean)
combine(groupby(data_fit, :Period), :home_distance_std => std)

#=
major_finder_std = []
major_finder_std = cc_standardize(16, major_finder_std)
rename!(major_finder_std,:x1 => :major_finder_std)
data_fit = innerjoin(data_fit, major_finder_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :major_finder_std => mean)
combine(groupby(data_fit, :Period), :major_finder_std => std)

SFS_std = []
SFS_std = cc_standardize(17, SFS_std)
rename!(SFS_std,:x1 => :SFS_std)
data_fit = innerjoin(data_fit, SFS_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :SFS_std => mean)
combine(groupby(data_fit, :Period), :SFS_std => std)
=#

#=
mean(data_fit.Diff_inst_grant_rate)
std(data_fit.Diff_inst_grant_rate)
diff_grant_std = []
diff_grant_std = cc_standardize(27, diff_grant_std)
rename!(diff_grant_std,:x1 => :diff_grant_std)
data_fit = innerjoin(data_fit, diff_grant_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :diff_grant_std => mean)
combine(groupby(data_fit, :Period), :diff_grant_std => std)
=#


###############################

```
Period 1 to Period 8
Only non-resident
Piece-wise exponential
MLE model without period
```

surv_fit = Model(Ipopt.Optimizer)

l = -10 .* ones(n_Period)  
u = 10 .* ones(n_Period)  
@variable(surv_fit, l[i] <= ??0[i = 1:n_Period] <= u[i])
#@variable(surv_fit, ??0[i = 1:n_Period] )

n_var = 16
l_?? = -10 .* ones(n_var)  
u_?? = 10 .* ones(n_var)  
@variable(surv_fit, l_??[i] <= ??[i = 1:n_var] <= u_??[i])


@NLexpression(
    surv_fit,
    X??[i = 1:n_fit], ??0[data_fit.Period[i]] 
                        + ??[1] * data_fit.FinAid_Rate[i]
                        + ??[2] * data_fit.Pell_Ind[i]
                        + ??[3] * data_fit.fed_efc_rate[i]
                        + ??[4] * data_fit.home_distance_std[i]
                        + ??[5] * data_fit.Gender_Ind[i]
                        + ??[6] * data_fit.Eth_ASIAN_Ind[i]
                        + ??[7] * data_fit.Eth_BLACK_Ind[i]
                        + ??[8] * data_fit.Eth_HISPA_Ind[i]
                        + ??[9] * data_fit.Eth_WHITE_Ind[i]
                        + ??[10] * data_fit.Eth_Multi_Ind[i]
                        + ??[11] * data_fit.Pros_Event_Ind[i]
                        + ??[12] * data_fit.Admit_Honor_Ind[i]
                        + ??[13] * data_fit.Diff_Major_Ind[i]
                        + ??[14] * data_fit.CampusTour_Ever_Ind[i]
                        + ??[15] * data_fit.DecisionDay_Ever_Ind[i]
                        + ??[16] * data_fit.Delay_Review_Ind[i]

    )

#=       
@NLconstraint(
    surv_fit, [i = 1:n_fit], 
    X??[i] >= 1e-4
)

### ??[i] = X??[i]
@NLobjective(
    surv_fit,
        Max,
        sum( data_fit.Dpst_Ind[i]*log(X??[i]) - X??[i]*data_fit.Period_length[i]  for i = 1:n_fit)
    )
=#

### ??[i] = exp(X??[i])


w = 1.05
#r = 0.00001
r = 0
@NLobjective(
    surv_fit,
        Min,
        -sum( w*data_fit.Dpst_Ind[i]*X??[i] - exp(X??[i])*data_fit.Period_length[i]  for i = 1:n_fit) / n_fit
        + r * sum( (??[j])^2  for j = 1:n_var)
        
             
    )

#print(surv_test)

optimize!(surv_fit)

termination_status(surv_fit)
objective_value(surv_fit)

value.(??0)
value.(??)

print(value.(??0))

para_df = DataFrame(
    ??0 = value.(??0)
    , ??_FinAid = value.(??)[1]
    , ??_Pell = value.(??)[2]
    , ??_efc = value.(??)[3]
    , ??_home = value.(??)[4]
    , ??_Gender = value.(??)[5]
    , ??_ASIAN = value.(??)[6]
    , ??_BLACK = value.(??)[7]
    , ??_HISPA = value.(??)[8]
    , ??_WHITE = value.(??)[9]
    , ??_Multi = value.(??)[10]
    , ??_Pros_Event = value.(??)[11]
    , ??_Admit_Honor = value.(??)[12]
    , ??_Diff_Major = value.(??)[13]
    , ??_CampusTour = value.(??)[14]
    , ??_DecisionDay = value.(??)[15]
    , ??_Delay_Review = value.(??)[16]
)

print(para_df)

#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_high_20221123.csv", para_df)

# para_df = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_Baseline_20220628.csv", DataFrame)

??0_fit = value.(??0)
??_fit = value.(??)

X??_fit = (??0_fit[data_fit.Period] 
.+ ??_fit[1] .* data_fit.FinAid_Rate
.+ ??_fit[2] .* data_fit.Pell_Ind
.+ ??_fit[3] .* data_fit.fed_efc_rate
.+ ??_fit[4] .* data_fit.home_distance_std
.+ ??_fit[5] .* data_fit.Gender_Ind
.+ ??_fit[6] .* data_fit.Eth_ASIAN_Ind
.+ ??_fit[7] .* data_fit.Eth_BLACK_Ind
.+ ??_fit[8] .* data_fit.Eth_HISPA_Ind
.+ ??_fit[9] .* data_fit.Eth_WHITE_Ind
.+ ??_fit[10] .* data_fit.Eth_Multi_Ind
.+ ??_fit[11] .* data_fit.Pros_Event_Ind
.+ ??_fit[12] .* data_fit.Admit_Honor_Ind
.+ ??_fit[13] .* data_fit.Diff_Major_Ind
.+ ??_fit[14] .* data_fit.CampusTour_Ever_Ind
.+ ??_fit[15] .* data_fit.DecisionDay_Ever_Ind
.+ ??_fit[16] .* data_fit.Delay_Review_Ind
)

??_fit = exp.(X??_fit)
S_fit = exp.(-??_fit .* data_fit.Period_length)
histogram(S_fit)
minimum(S_fit)

data_fit.S_fit = exp.(-??_fit .* data_fit.Period_length)
data_fit.??_fit = 1 .- data_fit.S_fit

sum(data_fit.??_fit)
sum(data_fit.Dpst_Ind)


combine(groupby(data_fit, :Period), :??_fit => sum)
combine(groupby(data_fit, :Period), :Dpst_Ind => sum)

Actual_vs_Fit = DataFrame(
                    Period = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Period
                    , Actual = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Dpst_Ind_sum
                    , Fit = combine(groupby(data_fit, :Period), :??_fit => sum).??_fit_sum )

#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Performance_2228_20220817.csv", Actual_vs_Fit)

###############################

```
Period 1 to Period 8
Only non-resident
Piece-wise exponential
MLE model without regulation
```

surv_fit = Model(Ipopt.Optimizer)

l = -10 .* ones(n_Period)  
u = 10 .* ones(n_Period)  
@variable(surv_fit, l[i] <= ??0[i = 1:n_Period] <= u[i])
#@variable(surv_fit, ??0[i = 1:n_Period] )

n_var = 16
l_?? = -10 .* ones(n_var)  
u_?? = 10 .* ones(n_var)  
@variable(surv_fit, l_??[i] <= ??[i = 1:n_Period, 1:n_var] <= u_??[i])

@NLexpression(
    surv_fit,
    X??[i = 1:n_fit], ??0[data_fit.Period[i]] 
                        + ??[data_fit.Period[i],1] * data_fit.FinAid_Rate[i]
                        + ??[data_fit.Period[i],2] * data_fit.Pell_Ind[i]
                        + ??[data_fit.Period[i],3] * data_fit.fed_efc_rate[i]
                        + ??[data_fit.Period[i],4] * data_fit.home_distance_std[i]
                        + ??[data_fit.Period[i],5] * data_fit.Gender_Ind[i]
                        + ??[data_fit.Period[i],6] * data_fit.Eth_ASIAN_Ind[i]
                        + ??[data_fit.Period[i],7] * data_fit.Eth_BLACK_Ind[i]
                        + ??[data_fit.Period[i],8] * data_fit.Eth_HISPA_Ind[i]
                        + ??[data_fit.Period[i],9] * data_fit.Eth_WHITE_Ind[i]
                        + ??[data_fit.Period[i],10] * data_fit.Eth_Multi_Ind[i]
                        + ??[data_fit.Period[i],11] * data_fit.Pros_Event_Ind[i]
                        + ??[data_fit.Period[i],12] * data_fit.Admit_Honor_Ind[i]
                        + ??[data_fit.Period[i],13] * data_fit.Diff_Major_Ind[i]
                        + ??[data_fit.Period[i],14] * data_fit.CampusTour_Ever_Ind[i]
                        + ??[data_fit.Period[i],15] * data_fit.DecisionDay_Ever_Ind[i]
                        + ??[data_fit.Period[i],16] * data_fit.Delay_Review_Ind[i]


    )

#=       
@NLconstraint(
    surv_fit, [i = 1:n_fit], 
    X??[i] >= 1e-4
)

### ??[i] = X??[i]
@NLobjective(
    surv_fit,
        Max,
        sum( data_fit.Dpst_Ind[i]*log(X??[i]) - X??[i]*data_fit.Period_length[i]  for i = 1:n_fit)
    )
=#

### ??[i] = exp(X??[i])

w = 1.05
#r = 0.00001
r = 0
@NLobjective(
    surv_fit,
        Min,
        -sum( w*data_fit.Dpst_Ind[i]*X??[i] - exp(X??[i])*data_fit.Period_length[i]  for i = 1:n_fit) / n_fit
        + r * sum( (??[j,k])^2  for j = 1:n_Period, k = 1:n_var)
        
             
    )

#print(surv_test)

optimize!(surv_fit)

termination_status(surv_fit)
objective_value(surv_fit)

value.(??0)
value.(??)

print(value.(??0))

para_df = DataFrame(
    ??0 = value.(??0)
    , ??_FinAid = value.(??)[:,1]
    , ??_Pell = value.(??)[:,2]
    , ??_efc = value.(??)[:,3]
    , ??_home = value.(??)[:,4]
    , ??_Gender = value.(??)[:,5]
    , ??_ASIAN = value.(??)[:,6]
    , ??_BLACK = value.(??)[:,7]
    , ??_HISPA = value.(??)[:,8]
    , ??_WHITE = value.(??)[:,9]
    , ??_Multi = value.(??)[:,10]
    , ??_Pros_Event = value.(??)[:,11]
    , ??_Admit_Honor = value.(??)[:,12]
    , ??_Diff_Major = value.(??)[:,13]
    , ??_CampusTour = value.(??)[:,14]
    , ??_DecisionDay = value.(??)[:,15]
    , ??_Delay_Review = value.(??)[:,16]
)

print(para_df)

#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_20221123.csv", para_df)

# para_df = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_Baseline_20220628.csv", DataFrame)

??0_fit = value.(??0)
??_fit = value.(??)

X??_fit = (??0_fit[data_fit.Period] 
.+ ??_fit[data_fit.Period,1] .* data_fit.FinAid_Rate
.+ ??_fit[data_fit.Period,2] .* data_fit.Pell_Ind
.+ ??_fit[data_fit.Period,3] .* data_fit.fed_efc_rate
.+ ??_fit[data_fit.Period,4] .* data_fit.home_distance_std
.+ ??_fit[data_fit.Period,5] .* data_fit.Gender_Ind
.+ ??_fit[data_fit.Period,6] .* data_fit.Eth_ASIAN_Ind
.+ ??_fit[data_fit.Period,7] .* data_fit.Eth_BLACK_Ind
.+ ??_fit[data_fit.Period,8] .* data_fit.Eth_HISPA_Ind
.+ ??_fit[data_fit.Period,9] .* data_fit.Eth_WHITE_Ind
.+ ??_fit[data_fit.Period,10] .* data_fit.Eth_Multi_Ind
.+ ??_fit[data_fit.Period,11] .* data_fit.Pros_Event_Ind
.+ ??_fit[data_fit.Period,12] .* data_fit.Admit_Honor_Ind
.+ ??_fit[data_fit.Period,13] .* data_fit.Diff_Major_Ind
.+ ??_fit[data_fit.Period,14] .* data_fit.CampusTour_Ever_Ind
.+ ??_fit[data_fit.Period,15] .* data_fit.DecisionDay_Ever_Ind
.+ ??_fit[data_fit.Period,16] .* data_fit.Delay_Review_Ind

)

??_fit = exp.(X??_fit)
S_fit = exp.(-??_fit .* data_fit.Period_length)
histogram(S_fit)
minimum(S_fit)

data_fit.S_fit = exp.(-??_fit .* data_fit.Period_length)
data_fit.??_fit = 1 .- data_fit.S_fit

sum(data_fit.??_fit)
sum(data_fit.Dpst_Ind)


combine(groupby(data_fit, :Period), :??_fit => sum)
combine(groupby(data_fit, :Period), :Dpst_Ind => sum)

Actual_vs_Fit = DataFrame(
                    Period = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Period
                    , Actual = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Dpst_Ind_sum
                    , Fit = combine(groupby(data_fit, :Period), :??_fit => sum).??_fit_sum )

#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Performance_2228_20220817.csv", Actual_vs_Fit)



#################################


```
Period 1 to Period 8
Only non-resident
Hierarchical
Dynamic ??
```

raw_MLE_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_20221123.csv", DataFrame)
raw_MLE_high_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_high_20221123.csv", DataFrame)
#raw_MLE_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2218_20221123.csv", DataFrame)
#raw_MLE_high_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2218_high_20221123.csv", DataFrame)
#raw_MLE_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2228_20221123.csv", DataFrame)
#raw_MLE_high_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2228_high_20221123.csv", DataFrame)

print(raw_MLE_para)

n_var = 16 ### This n_var contains ??0
n_Period

struct SurvivalProblem{Ty, TX}
    y::Ty
    X::TX
end

??0_MLE = raw_MLE_para.??0

function (problem::SurvivalProblem)(??)
  @unpack y, X = problem   # extract the data           
  @unpack ??0, ??, ??_high, ??_h, ??_l  = ??      # extract the parameters 
  

  ??_high = zeros(n_var)
  ??_0 = 0.1
  ??_h = 0.5
  ??_l = 0.5
  #??_h = 10.0
  #??_l = 10.0
  
  X?? = (
    ??0[X] 
    .+ ??[0*n_Period .+ X] .* data_fit.FinAid_Rate
    .+ ??[1*n_Period .+ X] .* data_fit.Pell_Ind
    .+ ??[2*n_Period .+ X] .* data_fit.fed_efc_rate
    .+ ??[3*n_Period .+ X] .* data_fit.home_distance_std
    .+ ??[4*n_Period .+ X] .* data_fit.Gender_Ind
    .+ ??[5*n_Period .+ X] .* data_fit.Eth_ASIAN_Ind
    .+ ??[6*n_Period .+ X] .* data_fit.Eth_BLACK_Ind
    .+ ??[7*n_Period .+ X] .* data_fit.Eth_HISPA_Ind
    .+ ??[8*n_Period .+ X] .* data_fit.Eth_WHITE_Ind
    .+ ??[9*n_Period .+ X] .* data_fit.Eth_Multi_Ind
    .+ ??[10*n_Period .+ X] .* data_fit.Pros_Event_Ind
    .+ ??[11*n_Period .+ X] .* data_fit.Admit_Honor_Ind
    .+ ??[12*n_Period .+ X] .* data_fit.Diff_Major_Ind
    .+ ??[13*n_Period .+ X] .* data_fit.CampusTour_Ever_Ind
    .+ ??[14*n_Period .+ X] .* data_fit.DecisionDay_Ever_Ind
    .+ ??[15*n_Period .+ X] .* data_fit.Delay_Review_Ind
    )

  w = 1.05  
  loglike = sum(w .*y .* X?? .- exp.(X??) .* data_fit.Period_length)
  logpri_??0 = sum(logpdf(MultivariateNormal(??0_MLE, ??_0), ??0))
  logpri_high = sum(logpdf(MultivariateNormal(??_high, ??_h), ??_high))
  logpri_??_h = sum(logpdf(Exponential(??_h), ??_h))
  
  ??_?? = []
  for i in 1:n_var
    ??_temp = ??_high[i]*ones(n_Period)
    if i == 1
        ??_?? = ??_temp
        else ??_?? = vcat(??_??, ??_temp)
    end
  end

  logpri_low = sum(logpdf(MultivariateNormal(??_??, ??_l), ??))
  logpri_??_l = sum(logpdf(Exponential(??_l), ??_l))

  loglike + logpri_??0 + logpri_high + logpri_low + logpri_??_h + logpri_??_l
end


p = SurvivalProblem(data_fit.Dpst_Ind, data_fit.Period)

p((??0 = zeros(n_Period), ?? = zeros(n_var*n_Period), ??_high = zeros(n_var), ??_h = 1.0, ??_l = 1.0 ))
p((??0 = zeros(n_Period), ?? = zeros(n_var*n_Period), ??_high = zeros(n_var), ??_h = 10.0, ??_l = 10.0 ))

p((??0 = ??0_MLE, ?? = zeros(n_var*n_Period), ??_high = zeros(n_var), ??_h = 1.0, ??_l = 1.0  ))
p((??0 = ??0_MLE, ?? = zeros(n_var*n_Period), ??_high = ones(n_var), ??_h = 1.0, ??_l = 1.0  ))

p((??0 = ??0_MLE, ?? = zeros(n_var*n_Period), ??_high = zeros(n_var), ??_h = 10.0, ??_l = 10.0  ))
p((??0 = ??0_MLE, ?? = zeros(n_var*n_Period), ??_high = ones(n_var), ??_h = 10.0, ??_l = 10.0  ))


??0_init = ones(n_Period)
??_init = zeros(n_var*n_Period)
??_h_init = zeros(n_var)

t = as((??0 = as(Array, length(??0_init)), ?? = as(Array, length(??_init)), ??_high = as(Array, length(??_h_init)), ??_h = as??????, ??_l = as??????))
P = TransformedLogDensity(t, p)
#???P = ADgradient(:ForwardDiff, P);
#???P = ADgradient(:ReverseDiff, P);
???P = ADgradient(:Zygote, P);

q??? = vcat(raw_MLE_para.??0
            , raw_MLE_para.??_FinAid
            , raw_MLE_para.??_Pell
            , raw_MLE_para.??_efc
            , raw_MLE_para.??_home
            , raw_MLE_para.??_Gender
            , raw_MLE_para.??_ASIAN
            , raw_MLE_para.??_BLACK
            , raw_MLE_para.??_HISPA
            , raw_MLE_para.??_WHITE 
            , raw_MLE_para.??_Multi
            , raw_MLE_para.??_Pros_Event
            , raw_MLE_para.??_Admit_Honor
            , raw_MLE_para.??_Diff_Major
            , raw_MLE_para.??_CampusTour
            , raw_MLE_para.??_DecisionDay
            , raw_MLE_para.??_Delay_Review

            , raw_MLE_high_para.??_FinAid[1]
            , raw_MLE_high_para.??_Pell[1]
            , raw_MLE_high_para.??_efc[1]
            , raw_MLE_high_para.??_home[1]
            , raw_MLE_high_para.??_Gender[1]
            , raw_MLE_high_para.??_ASIAN[1]
            , raw_MLE_high_para.??_BLACK[1]
            , raw_MLE_high_para.??_HISPA[1]
            , raw_MLE_high_para.??_WHITE[1] 
            , raw_MLE_high_para.??_Multi[1]
            , raw_MLE_high_para.??_Pros_Event[1]
            , raw_MLE_high_para.??_Admit_Honor[1]
            , raw_MLE_high_para.??_Diff_Major[1]
            , raw_MLE_high_para.??_CampusTour[1]
            , raw_MLE_high_para.??_DecisionDay[1]
            , raw_MLE_high_para.??_Delay_Review[1]

            , log(1.0)
            , log(0.2) )

results = mcmc_with_warmup(Random.GLOBAL_RNG, ???P, 1000;
                            initialization = (q=q???, ) )


# To get the posterior for ``??``, we need to use `get_position` and
# then transform

#posterior = transform.(t, results.chain);
posterior = results.chain

posterior_?? = first.(posterior);

# check the mean

mean(posterior_??)

# check the effective sample size

ess_?? = effective_sample_size(posterior_??)

# NUTS-specific statistics

summarize_tree_statistics(results.tree_statistics)

#=
Array(results.chain)
mean(results.chain, dims = 1)
mean(DataFrame(results.chain,:auto)[1,:])
=#

chain_df = DataFrame(results.chain,:auto)
plot(Vector(chain_df[15,:]))
plot(Vector(chain_df[25,:]))
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2208_20221123.csv", chain_df)

print(mean(Matrix(chain_df), dims = 2))
??_chain_df = mean(Matrix(chain_df), dims = 2)

n_?? = (n_var+1)*n_Period
??0_fit = ??_chain_df[1:n_Period]
??_fit = reshape(??_chain_df[(n_Period+1):n_??], (n_Period,n_var))
??_h_fit = ??_chain_df[(n_??+1):(n_??+n_var)]
??_h_fit = exp(??_chain_df[n_??+n_var+1])
??_l_fit = exp(??_chain_df[n_??+n_var+2])
# exp(??_h_init)
# exp(??_l_init)


X??_fit = (??0_fit[data_fit.Period] 
.+ ??_fit[data_fit.Period,1] .* data_fit.FinAid_Rate
.+ ??_fit[data_fit.Period,2] .* data_fit.Pell_Ind
.+ ??_fit[data_fit.Period,3] .* data_fit.fed_efc_rate
.+ ??_fit[data_fit.Period,4] .* data_fit.home_distance_std
.+ ??_fit[data_fit.Period,5] .* data_fit.Gender_Ind
.+ ??_fit[data_fit.Period,6] .* data_fit.Eth_ASIAN_Ind
.+ ??_fit[data_fit.Period,7] .* data_fit.Eth_BLACK_Ind
.+ ??_fit[data_fit.Period,8] .* data_fit.Eth_HISPA_Ind
.+ ??_fit[data_fit.Period,9] .* data_fit.Eth_WHITE_Ind
.+ ??_fit[data_fit.Period,10] .* data_fit.Eth_Multi_Ind
.+ ??_fit[data_fit.Period,11] .* data_fit.Pros_Event_Ind
.+ ??_fit[data_fit.Period,12] .* data_fit.Admit_Honor_Ind
.+ ??_fit[data_fit.Period,13] .* data_fit.Diff_Major_Ind
.+ ??_fit[data_fit.Period,14] .* data_fit.CampusTour_Ever_Ind
.+ ??_fit[data_fit.Period,15] .* data_fit.DecisionDay_Ever_Ind
.+ ??_fit[data_fit.Period,16] .* data_fit.Delay_Review_Ind
)


??_fit = exp.(X??_fit)
S_fit = exp.(-??_fit .* data_fit.Period_length)
histogram(S_fit)
minimum(S_fit)

data_fit.S_fit = exp.(-??_fit .* data_fit.Period_length)
data_fit.??_fit = 1 .- data_fit.S_fit

sum(data_fit.??_fit)
sum(data_fit.Dpst_Ind)


combine(groupby(data_fit, :Period), :??_fit => sum)
combine(groupby(data_fit, :Period), :Dpst_Ind => sum)

Actual_vs_Fit = DataFrame(
                    Period = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Period
                    , Actual = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Dpst_Ind_sum
                    , Fit = combine(groupby(data_fit, :Period), :??_fit => sum).??_fit_sum )



#################################