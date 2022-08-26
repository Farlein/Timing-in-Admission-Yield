

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

# using JuMP, Ipopt


using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics
using MCMCDiagnostics
using Parameters, Statistics
import ForwardDiff  
import ReverseDiff
import Zygote


```
Load raw data

```

###
#raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2208_20220816.csv", DataFrame)
raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2218_20220816.csv", DataFrame)
#raw = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Data/for_julia_2228_20220816.csv", DataFrame)

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

data_fit.Delay_Review_Ind = zeros(nrow(data_fit))
for i in eachindex(data_fit.Delay_Review_Decision)
    if data_fit.Delay_Review_Decision[i] >= 2
        data_fit.Delay_Review_Ind[i] = 1
        else data_fit.Delay_Review_Ind[i] = 0
    end
end

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


hs_gpa_std = []
hs_gpa_std = cc_standardize(10, hs_gpa_std)
rename!(hs_gpa_std,:x1 => :hs_gpa_std)
data_fit = innerjoin(data_fit, hs_gpa_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :hs_gpa_std => mean)
combine(groupby(data_fit, :Period), :hs_gpa_std => std)

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

mean(data_fit.Diff_inst_grant_rate)
std(data_fit.Diff_inst_grant_rate)
diff_grant_std = []
diff_grant_std = cc_standardize(27, diff_grant_std)
rename!(diff_grant_std,:x1 => :diff_grant_std)
data_fit = innerjoin(data_fit, diff_grant_std, on = [:Student_ID => :Student_ID, :Period => :Period])
combine(groupby(data_fit, :Period), :diff_grant_std => mean)
combine(groupby(data_fit, :Period), :diff_grant_std => std)



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
@variable(surv_fit, l[i] <= β0[i = 1:n_Period] <= u[i])
#@variable(surv_fit, β0[i = 1:n_Period] )

n_var = 19
l_β = -10 .* ones(n_var)  
u_β = 10 .* ones(n_var)  
@variable(surv_fit, l_β[i] <= β[i = 1:n_var] <= u_β[i])

@NLexpression(
    surv_fit,
    Xβ[i = 1:n_fit], β0[data_fit.Period[i]] 
                        + β[1] * data_fit.Admit_Ind[i]
                        + β[2] * data_fit.hs_gpa_std[i]
                        + β[3] * data_fit.Admit_Honor_Ind[i]
                        + β[4] * data_fit.Diff_Major_Ind[i]
                        + β[5] * data_fit.Gender_Ind[i]
                        + β[6] * data_fit.inst_grant_rate[i]
                        + β[7] * data_fit.student_loan_rate[i]
                        + β[8] * data_fit.fed_efc_rate[i]
                        + β[9] * data_fit.Pell_Ind[i]
                        + β[10] * data_fit.Eth_ASIAN_Ind[i]
                        + β[11] * data_fit.Eth_BLACK_Ind[i]
                        + β[12] * data_fit.Eth_HISPA_Ind[i]
                        + β[13] * data_fit.Eth_WHITE_Ind[i]
                        + β[14] * data_fit.Eth_Multi_Ind[i]
                        + β[15] * data_fit.major_finder_std[i]
                        + β[16] * data_fit.SFS_std[i]
                        + β[17] * data_fit.Blue_Gold_Ind[i]
                        + β[18] * data_fit.Early_Event_Ind[i]
                        + β[19] * data_fit.Delay_Review_Ind[i]

    )

#=       
@NLconstraint(
    surv_fit, [i = 1:n_fit], 
    Xβ[i] >= 1e-4
)

### λ[i] = Xβ[i]
@NLobjective(
    surv_fit,
        Max,
        sum( data_fit.Dpst_Ind[i]*log(Xβ[i]) - Xβ[i]*data_fit.Period_length[i]  for i = 1:n_fit)
    )
=#

### λ[i] = exp(Xβ[i])

#r = 0.00001
r = 0
@NLobjective(
    surv_fit,
        Min,
        -sum( data_fit.Dpst_Ind[i]*Xβ[i] - exp(Xβ[i])*data_fit.Period_length[i]  for i = 1:n_fit) / n_fit
        + r * sum( (β[j])^2  for j = 1:n_var)
        
             
    )

#print(surv_test)

optimize!(surv_fit)

termination_status(surv_fit)
objective_value(surv_fit)

value.(β0)
value.(β)

print(value.(β0))

para_df = DataFrame(
    β0 = value.(β0)
    , β_Admit = value.(β)[1]
    , β_hs_gpa = value.(β)[2]
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
    , β_Major_Finder = value.(β)[15]
    , β_SFS = value.(β)[16]
    , β_Blue_Gold = value.(β)[17]
    , β_Early_Event = value.(β)[18]
    , β_Delay_Review = value.(β)[19]
)

print(para_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/results/para_2208_Adm_Bio_20220628.csv", para_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_high_20220817.csv", para_df)

# para_df = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_Baseline_20220628.csv", DataFrame)

β0_fit = value.(β0)
β_fit = value.(β)

#=
    Xβ[i = 1:n_fit], β0[data_fit.Period[i]] 
                        + β[data_fit.Period[i],1] * data_fit.Admit_Ind[i]
                        + β[data_fit.Period[i],2] * data_fit.hs_gpa_std[i]
                        + β[data_fit.Period[i],3] * data_fit.Admit_Honor_Ind[i]
                        + β[data_fit.Period[i],4] * data_fit.Diff_Major_Ind[i]
                        + β[data_fit.Period[i],5] * data_fit.Gender_Ind[i]
                        + β[data_fit.Period[i],6] * data_fit.inst_grant_rate[i]
                        + β[data_fit.Period[i],7] * data_fit.student_loan_rate[i]
                        + β[data_fit.Period[i],8] * data_fit.fed_efc_rate[i]
                        + β[data_fit.Period[i],9] * data_fit.Pell_Ind[i]
                        + β[data_fit.Period[i],10] * data_fit.Eth_ASIAN_Ind[i]
                        + β[data_fit.Period[i],11] * data_fit.Eth_BLACK_Ind[i]
                        + β[data_fit.Period[i],12] * data_fit.Eth_HISPA_Ind[i]
                        + β[data_fit.Period[i],13] * data_fit.Eth_WHITE_Ind[i]
                        + β[data_fit.Period[i],14] * data_fit.Eth_Multi_Ind[i]
                        + β[data_fit.Period[i],15] * data_fit.major_finder_std[i]
                        + β[data_fit.Period[i],16] * data_fit.SFS_std[i]
                        + β[data_fit.Period[i],17] * data_fit.Blue_Gold_Ind[i]
                        + β[data_fit.Period[i],18] * data_fit.Early_Event_Ind[i]
                        + β[data_fit.Period[i],19] * data_fit.Delay_Review_Ind[i]
=#

Xβ_fit = (β0_fit[data_fit.Period] 
.+ β_fit[1] .* data_fit.Admit_Ind
.+ β_fit[2] .* data_fit.hs_gpa_std
.+ β_fit[3] .* data_fit.Admit_Honor_Ind
.+ β_fit[4] .* data_fit.Diff_Major_Ind
.+ β_fit[5] .* data_fit.Gender_Ind
.+ β_fit[6] .* data_fit.inst_grant_rate
.+ β_fit[7] .* data_fit.student_loan_rate
.+ β_fit[8] .* data_fit.fed_efc_rate
.+ β_fit[9] .* data_fit.Pell_Ind
.+ β_fit[10] .* data_fit.Eth_ASIAN_Ind
.+ β_fit[11] .* data_fit.Eth_BLACK_Ind
.+ β_fit[12] .* data_fit.Eth_HISPA_Ind
.+ β_fit[13] .* data_fit.Eth_WHITE_Ind
.+ β_fit[14] .* data_fit.Eth_Multi_Ind
.+ β_fit[15] .* data_fit.major_finder_std
.+ β_fit[16] .* data_fit.SFS_std
.+ β_fit[17] .* data_fit.Blue_Gold_Ind
.+ β_fit[18] .* data_fit.Early_Event_Ind
.+ β_fit[19] .* data_fit.Delay_Review_Ind

)

λ_fit = exp.(Xβ_fit)
S_fit = exp.(-λ_fit .* data_fit.Period_length)
histogram(S_fit)
minimum(S_fit)

data_fit.S_fit = exp.(-λ_fit .* data_fit.Period_length)
data_fit.θ_fit = 1 .- data_fit.S_fit

sum(data_fit.θ_fit)
sum(data_fit.Dpst_Ind)


combine(groupby(data_fit, :Period), :θ_fit => sum)
combine(groupby(data_fit, :Period), :Dpst_Ind => sum)

Actual_vs_Fit = DataFrame(
                    Period = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Period
                    , Actual = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Dpst_Ind_sum
                    , Fit = combine(groupby(data_fit, :Period), :θ_fit => sum).θ_fit_sum )

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
@variable(surv_fit, l[i] <= β0[i = 1:n_Period] <= u[i])
#@variable(surv_fit, β0[i = 1:n_Period] )

n_var = 19
l_β = -10 .* ones(n_var)  
u_β = 10 .* ones(n_var)  
@variable(surv_fit, l_β[i] <= β[i = 1:n_Period, 1:n_var] <= u_β[i])

@NLexpression(
    surv_fit,
    Xβ[i = 1:n_fit], β0[data_fit.Period[i]] 
                        + β[data_fit.Period[i],1] * data_fit.Admit_Ind[i]
                        + β[data_fit.Period[i],2] * data_fit.hs_gpa_std[i]
                        + β[data_fit.Period[i],3] * data_fit.Admit_Honor_Ind[i]
                        + β[data_fit.Period[i],4] * data_fit.Diff_Major_Ind[i]
                        + β[data_fit.Period[i],5] * data_fit.Gender_Ind[i]
                        + β[data_fit.Period[i],6] * data_fit.inst_grant_rate[i]
                        + β[data_fit.Period[i],7] * data_fit.student_loan_rate[i]
                        + β[data_fit.Period[i],8] * data_fit.fed_efc_rate[i]
                        + β[data_fit.Period[i],9] * data_fit.Pell_Ind[i]
                        + β[data_fit.Period[i],10] * data_fit.Eth_ASIAN_Ind[i]
                        + β[data_fit.Period[i],11] * data_fit.Eth_BLACK_Ind[i]
                        + β[data_fit.Period[i],12] * data_fit.Eth_HISPA_Ind[i]
                        + β[data_fit.Period[i],13] * data_fit.Eth_WHITE_Ind[i]
                        + β[data_fit.Period[i],14] * data_fit.Eth_Multi_Ind[i]
                        + β[data_fit.Period[i],15] * data_fit.major_finder_std[i]
                        + β[data_fit.Period[i],16] * data_fit.SFS_std[i]
                        + β[data_fit.Period[i],17] * data_fit.Blue_Gold_Ind[i]
                        + β[data_fit.Period[i],18] * data_fit.Early_Event_Ind[i]
                        + β[data_fit.Period[i],19] * data_fit.Delay_Review_Ind[i]

    )

#=       
@NLconstraint(
    surv_fit, [i = 1:n_fit], 
    Xβ[i] >= 1e-4
)

### λ[i] = Xβ[i]
@NLobjective(
    surv_fit,
        Max,
        sum( data_fit.Dpst_Ind[i]*log(Xβ[i]) - Xβ[i]*data_fit.Period_length[i]  for i = 1:n_fit)
    )
=#

### λ[i] = exp(Xβ[i])

#r = 0.00001
r = 0
@NLobjective(
    surv_fit,
        Min,
        -sum( data_fit.Dpst_Ind[i]*Xβ[i] - exp(Xβ[i])*data_fit.Period_length[i]  for i = 1:n_fit) / n_fit
        + r * sum( (β[j,k])^2  for j = 1:n_Period, k = 1:n_var)
        
             
    )

#print(surv_test)

optimize!(surv_fit)

termination_status(surv_fit)
objective_value(surv_fit)

value.(β0)
value.(β)

print(value.(β0))

para_df = DataFrame(
    β0 = value.(β0)
    , β_Admit = value.(β)[:,1]
    , β_hs_gpa = value.(β)[:,2]
    , β_Admit_Honor = value.(β)[:,3]
    , β_Diff_Major = value.(β)[:,4]
    , β_Gender = value.(β)[:,5]
    , β_inst_grant = value.(β)[:,6]
    , β_loan = value.(β)[:,7]
    , β_fed_efc = value.(β)[:,8]
    , β_Pell = value.(β)[:,9]
    , β_ASIAN = value.(β)[:,10]
    , β_BLACK = value.(β)[:,11]
    , β_HISPA = value.(β)[:,12]
    , β_WHITE = value.(β)[:,13]
    , β_Multi = value.(β)[:,14]
    , β_Major_Finder = value.(β)[:,15]
    , β_SFS = value.(β)[:,16]
    , β_Blue_Gold = value.(β)[:,17]
    , β_Early_Event = value.(β)[:,18]
    , β_Delay_Review = value.(β)[:,19]
)

print(para_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/results/para_2208_Adm_Bio_20220628.csv", para_df)
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_20220817.csv", para_df)

# para_df = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2208_Baseline_20220628.csv", DataFrame)

β0_fit = value.(β0)
β_fit = value.(β)

#=
    Xβ[i = 1:n_fit], β0[data_fit.Period[i]] 
                        + β[data_fit.Period[i],1] * data_fit.Admit_Ind[i]
                        + β[data_fit.Period[i],2] * data_fit.hs_gpa_std[i]
                        + β[data_fit.Period[i],3] * data_fit.Admit_Honor_Ind[i]
                        + β[data_fit.Period[i],4] * data_fit.Diff_Major_Ind[i]
                        + β[data_fit.Period[i],5] * data_fit.Gender_Ind[i]
                        + β[data_fit.Period[i],6] * data_fit.inst_grant_rate[i]
                        + β[data_fit.Period[i],7] * data_fit.student_loan_rate[i]
                        + β[data_fit.Period[i],8] * data_fit.fed_efc_rate[i]
                        + β[data_fit.Period[i],9] * data_fit.Pell_Ind[i]
                        + β[data_fit.Period[i],10] * data_fit.Eth_ASIAN_Ind[i]
                        + β[data_fit.Period[i],11] * data_fit.Eth_BLACK_Ind[i]
                        + β[data_fit.Period[i],12] * data_fit.Eth_HISPA_Ind[i]
                        + β[data_fit.Period[i],13] * data_fit.Eth_WHITE_Ind[i]
                        + β[data_fit.Period[i],14] * data_fit.Eth_Multi_Ind[i]
                        + β[data_fit.Period[i],15] * data_fit.major_finder_std[i]
                        + β[data_fit.Period[i],16] * data_fit.SFS_std[i]
                        + β[data_fit.Period[i],17] * data_fit.Blue_Gold_Ind[i]
                        + β[data_fit.Period[i],18] * data_fit.Early_Event_Ind[i]
                        + β[data_fit.Period[i],19] * data_fit.Delay_Review_Ind[i]
=#

Xβ_fit = (β0_fit[data_fit.Period] 
.+ β_fit[data_fit.Period,1] .* data_fit.Admit_Ind
.+ β_fit[data_fit.Period,2] .* data_fit.hs_gpa_std
.+ β_fit[data_fit.Period,3] .* data_fit.Admit_Honor_Ind
.+ β_fit[data_fit.Period,4] .* data_fit.Diff_Major_Ind
.+ β_fit[data_fit.Period,5] .* data_fit.Gender_Ind
.+ β_fit[data_fit.Period,6] .* data_fit.inst_grant_rate
.+ β_fit[data_fit.Period,7] .* data_fit.student_loan_rate
.+ β_fit[data_fit.Period,8] .* data_fit.fed_efc_rate
.+ β_fit[data_fit.Period,9] .* data_fit.Pell_Ind
.+ β_fit[data_fit.Period,10] .* data_fit.Eth_ASIAN_Ind
.+ β_fit[data_fit.Period,11] .* data_fit.Eth_BLACK_Ind
.+ β_fit[data_fit.Period,12] .* data_fit.Eth_HISPA_Ind
.+ β_fit[data_fit.Period,13] .* data_fit.Eth_WHITE_Ind
.+ β_fit[data_fit.Period,14] .* data_fit.Eth_Multi_Ind
.+ β_fit[data_fit.Period,15] .* data_fit.major_finder_std
.+ β_fit[data_fit.Period,16] .* data_fit.SFS_std
.+ β_fit[data_fit.Period,17] .* data_fit.Blue_Gold_Ind
.+ β_fit[data_fit.Period,18] .* data_fit.Early_Event_Ind
.+ β_fit[data_fit.Period,19] .* data_fit.Delay_Review_Ind

)

λ_fit = exp.(Xβ_fit)
S_fit = exp.(-λ_fit .* data_fit.Period_length)
histogram(S_fit)
minimum(S_fit)

data_fit.S_fit = exp.(-λ_fit .* data_fit.Period_length)
data_fit.θ_fit = 1 .- data_fit.S_fit

sum(data_fit.θ_fit)
sum(data_fit.Dpst_Ind)


combine(groupby(data_fit, :Period), :θ_fit => sum)
combine(groupby(data_fit, :Period), :Dpst_Ind => sum)

Actual_vs_Fit = DataFrame(
                    Period = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Period
                    , Actual = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Dpst_Ind_sum
                    , Fit = combine(groupby(data_fit, :Period), :θ_fit => sum).θ_fit_sum )

#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/Performance_2228_20220817.csv", Actual_vs_Fit)



#################################


```
Period 1 to Period 8
Only non-resident
Hierarchical
Dynamic σ
```

raw_MLE_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2218_20220817.csv", DataFrame)
raw_MLE_high_para = CSV.read("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/para_2218_high_20220817.csv", DataFrame)
print(raw_MLE_para)

n_var = 19
n_Period

struct SurvivalProblem{Ty, TX}
    y::Ty
    X::TX
end

β0_MLE = raw_MLE_para.β0

function (problem::SurvivalProblem)(θ)
  @unpack y, X = problem   # extract the data           
  @unpack β0, β, β_high, σ_h, σ_l  = θ      # extract the parameters 
  

  μ_high = zeros(n_var)
  σ_0 = 0.1
  λ_h = 0.5
  λ_l = 0.5
  #σ_h = 10.0
  #σ_l = 10.0

  Xβ = (
    β0[X] 
    .+ β[0*n_Period .+ X] .* data_fit.Admit_Ind 
    .+ β[1*n_Period .+ X] .* data_fit.hs_gpa_std
    .+ β[2*n_Period .+ X] .* data_fit.Admit_Honor_Ind 
    .+ β[3*n_Period .+ X] .* data_fit.Diff_Major_Ind 
    .+ β[4*n_Period .+ X] .* data_fit.Gender_Ind 
    .+ β[5*n_Period .+ X] .* data_fit.inst_grant_rate 
    .+ β[6*n_Period .+ X] .* data_fit.student_loan_rate 
    .+ β[7*n_Period .+ X] .* data_fit.fed_efc_rate 
    .+ β[8*n_Period .+ X] .* data_fit.Pell_Ind 
    .+ β[9*n_Period .+ X] .* data_fit.Eth_ASIAN_Ind 
    .+ β[10*n_Period .+ X] .* data_fit.Eth_BLACK_Ind 
    .+ β[11*n_Period .+ X] .* data_fit.Eth_HISPA_Ind 
    .+ β[12*n_Period .+ X] .* data_fit.Eth_WHITE_Ind 
    .+ β[13*n_Period .+ X] .* data_fit.Eth_Multi_Ind 
    .+ β[14*n_Period .+ X] .* data_fit.major_finder_std
    .+ β[15*n_Period .+ X] .* data_fit.SFS_std
    .+ β[16*n_Period .+ X] .* data_fit.Blue_Gold_Ind
    .+ β[17*n_Period .+ X] .* data_fit.Early_Event_Ind
    .+ β[18*n_Period .+ X] .* data_fit.Delay_Review_Ind
    )

  loglike = sum(y .* Xβ .- exp.(Xβ) .* data_fit.Period_length)
  logpri_β0 = sum(logpdf(MultivariateNormal(β0_MLE, σ_0), β0))
  logpri_high = sum(logpdf(MultivariateNormal(μ_high, σ_h), β_high))
  logpri_σ_h = sum(logpdf(Exponential(λ_h), σ_h))
  
  μ_β = []
  for i in 1:n_var
    μ_temp = β_high[i]*ones(n_Period)
    if i == 1
        μ_β = μ_temp
        else μ_β = vcat(μ_β, μ_temp)
    end
  end

  logpri_low = sum(logpdf(MultivariateNormal(μ_β, σ_l), β))
  logpri_σ_l = sum(logpdf(Exponential(λ_l), σ_l))

  loglike + logpri_β0 + logpri_high + logpri_low + logpri_σ_h + logpri_σ_l
end


p = SurvivalProblem(data_fit.Dpst_Ind, data_fit.Period)

p((β0 = zeros(n_Period), β = zeros(n_var*n_Period), β_high = zeros(n_var), σ_h = 1.0, σ_l = 1.0 ))
p((β0 = zeros(n_Period), β = zeros(n_var*n_Period), β_high = zeros(n_var), σ_h = 10.0, σ_l = 10.0 ))

p((β0 = β0_MLE, β = zeros(n_var*n_Period), β_high = zeros(n_var), σ_h = 1.0, σ_l = 1.0  ))
p((β0 = β0_MLE, β = zeros(n_var*n_Period), β_high = ones(n_var), σ_h = 1.0, σ_l = 1.0  ))

p((β0 = β0_MLE, β = zeros(n_var*n_Period), β_high = zeros(n_var), σ_h = 10.0, σ_l = 10.0  ))
p((β0 = β0_MLE, β = zeros(n_var*n_Period), β_high = ones(n_var), σ_h = 10.0, σ_l = 10.0  ))


β0_init = ones(n_Period)
β_init = zeros(n_var*n_Period)
β_h_init = zeros(n_var)

t = as((β0 = as(Array, length(β0_init)), β = as(Array, length(β_init)), β_high = as(Array, length(β_h_init)), σ_h = asℝ₊, σ_l = asℝ₊))
P = TransformedLogDensity(t, p)
#∇P = ADgradient(:ForwardDiff, P);
#∇P = ADgradient(:ReverseDiff, P);
∇P = ADgradient(:Zygote, P);



q₀ = vcat(raw_MLE_para.β0
            , raw_MLE_para.β_Admit
            , raw_MLE_para.β_hs_gpa
            , raw_MLE_para.β_Admit_Honor
            , raw_MLE_para.β_Diff_Major
            , raw_MLE_para.β_Gender
            , raw_MLE_para.β_inst_grant
            , raw_MLE_para.β_loan
            , raw_MLE_para.β_fed_efc
            , raw_MLE_para.β_Pell 
            , raw_MLE_para.β_ASIAN
            , raw_MLE_para.β_BLACK
            , raw_MLE_para.β_HISPA
            , raw_MLE_para.β_WHITE
            , raw_MLE_para.β_Multi
            , raw_MLE_para.β_Major_Finder
            , raw_MLE_para.β_SFS
            , raw_MLE_para.β_Blue_Gold
            , raw_MLE_para.β_Early_Event
            , raw_MLE_para.β_Delay_Review

            , raw_MLE_high_para.β_Admit[1]
            , raw_MLE_high_para.β_hs_gpa[1]
            , raw_MLE_high_para.β_Admit_Honor[1]
            , raw_MLE_high_para.β_Diff_Major[1]
            , raw_MLE_high_para.β_Gender[1]
            , raw_MLE_high_para.β_inst_grant[1]
            , raw_MLE_high_para.β_loan[1]
            , raw_MLE_high_para.β_fed_efc[1]
            , raw_MLE_high_para.β_Pell[1] 
            , raw_MLE_high_para.β_ASIAN[1]
            , raw_MLE_high_para.β_BLACK[1]
            , raw_MLE_high_para.β_HISPA[1]
            , raw_MLE_high_para.β_WHITE[1]
            , raw_MLE_high_para.β_Multi[1]
            , raw_MLE_high_para.β_Major_Finder[1]
            , raw_MLE_high_para.β_SFS[1]
            , raw_MLE_high_para.β_Blue_Gold[1]
            , raw_MLE_high_para.β_Early_Event[1]
            , raw_MLE_high_para.β_Delay_Review[1] 

            , log(1.0)
            , log(0.2) )

results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, 1000;
                            initialization = (q=q₀, ) )


# To get the posterior for ``α``, we need to use `get_position` and
# then transform

#posterior = transform.(t, results.chain);
posterior = results.chain

posterior_α = first.(posterior);

# check the mean

mean(posterior_α)

# check the effective sample size

ess_α = effective_sample_size(posterior_α)

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
#CSV.write("H:/My Drive/FSAN/5_Adm Yield Proj/Temp results/chain_2218_20220819.csv", chain_df)

print(mean(Matrix(chain_df), dims = 2))
μ_chain_df = mean(Matrix(chain_df), dims = 2)

n_β = (n_var+1)*n_Period
β0_fit = μ_chain_df[1:n_Period]
β_fit = reshape(μ_chain_df[(n_Period+1):n_β], (n_Period,n_var))
β_h_fit = μ_chain_df[(n_β+1):(n_β+n_var)]
σ_h_fit = exp(μ_chain_df[n_β+n_var+1])
σ_l_fit = exp(μ_chain_df[n_β+n_var+2])
# exp(σ_h_init)
# exp(σ_l_init)

#=
    β0[X] 
    .+ β[0*n_Period .+ X] .* data_fit.Admit_Ind 
    .+ β[1*n_Period .+ X] .* data_fit.hs_gpa_std
    .+ β[2*n_Period .+ X] .* data_fit.Admit_Honor_Ind 
    .+ β[3*n_Period .+ X] .* data_fit.Diff_Major_Ind 
    .+ β[4*n_Period .+ X] .* data_fit.Gender_Ind 
    .+ β[5*n_Period .+ X] .* data_fit.inst_grant_rate 
    .+ β[6*n_Period .+ X] .* data_fit.student_loan_rate 
    .+ β[7*n_Period .+ X] .* data_fit.fed_efc_rate 
    .+ β[8*n_Period .+ X] .* data_fit.Pell_Ind 
    .+ β[9*n_Period .+ X] .* data_fit.Eth_ASIAN_Ind 
    .+ β[10*n_Period .+ X] .* data_fit.Eth_BLACK_Ind 
    .+ β[11*n_Period .+ X] .* data_fit.Eth_HISPA_Ind 
    .+ β[12*n_Period .+ X] .* data_fit.Eth_WHITE_Ind 
    .+ β[13*n_Period .+ X] .* data_fit.Eth_Multi_Ind 
    .+ β[14*n_Period .+ X] .* data_fit.major_finder_std
    .+ β[15*n_Period .+ X] .* data_fit.SFS_std
    .+ β[16*n_Period .+ X] .* data_fit.Blue_Gold_Ind
    .+ β[17*n_Period .+ X] .* data_fit.Early_Event_Ind
    .+ β[18*n_Period .+ X] .* data_fit.Delay_Review_Ind
=#

Xβ_fit = (β0_fit[data_fit.Period] 
.+ β_fit[data_fit.Period,1] .* data_fit.Admit_Ind
.+ β_fit[data_fit.Period,2] .* data_fit.hs_gpa_std
.+ β_fit[data_fit.Period,3] .* data_fit.Admit_Honor_Ind
.+ β_fit[data_fit.Period,4] .* data_fit.Diff_Major_Ind
.+ β_fit[data_fit.Period,5] .* data_fit.Gender_Ind
.+ β_fit[data_fit.Period,6] .* data_fit.inst_grant_rate
.+ β_fit[data_fit.Period,7] .* data_fit.student_loan_rate
.+ β_fit[data_fit.Period,8] .* data_fit.fed_efc_rate
.+ β_fit[data_fit.Period,9] .* data_fit.Pell_Ind
.+ β_fit[data_fit.Period,10] .* data_fit.Eth_ASIAN_Ind
.+ β_fit[data_fit.Period,11] .* data_fit.Eth_BLACK_Ind
.+ β_fit[data_fit.Period,12] .* data_fit.Eth_HISPA_Ind
.+ β_fit[data_fit.Period,13] .* data_fit.Eth_WHITE_Ind
.+ β_fit[data_fit.Period,14] .* data_fit.Eth_Multi_Ind
.+ β_fit[data_fit.Period,15] .* data_fit.major_finder_std
.+ β_fit[data_fit.Period,16] .* data_fit.SFS_std
.+ β_fit[data_fit.Period,17] .* data_fit.Blue_Gold_Ind
.+ β_fit[data_fit.Period,18] .* data_fit.Early_Event_Ind
.+ β_fit[data_fit.Period,19] .* data_fit.Delay_Review_Ind


)


λ_fit = exp.(Xβ_fit)
S_fit = exp.(-λ_fit .* data_fit.Period_length)
histogram(S_fit)
minimum(S_fit)

data_fit.S_fit = exp.(-λ_fit .* data_fit.Period_length)
data_fit.θ_fit = 1 .- data_fit.S_fit

sum(data_fit.θ_fit)
sum(data_fit.Dpst_Ind)


combine(groupby(data_fit, :Period), :θ_fit => sum)
combine(groupby(data_fit, :Period), :Dpst_Ind => sum)

Actual_vs_Fit = DataFrame(
                    Period = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Period
                    , Actual = combine(groupby(data_fit, :Period), :Dpst_Ind => sum).Dpst_Ind_sum
                    , Fit = combine(groupby(data_fit, :Period), :θ_fit => sum).θ_fit_sum )



#################################