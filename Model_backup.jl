



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

using JuMP, Ipopt



```
Load raw data

```

###
#raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2208_20221027.csv", DataFrame)
#raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2218_20221027.csv", DataFrame)
raw = CSV.read("C:/Users/chuanc/University of Delaware - o365/Team-IRE-Staff Shares - chuanc - chuanc/Project/02_Analytical/20200102 ADM_Yield/Data/for_julia_2228_20221027.csv", DataFrame)

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

print(hcat(1:38, names(data_fit)))

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

n_var = 17
l_β = -10 .* ones(n_var)  
u_β = 10 .* ones(n_var)  
@variable(surv_fit, l_β[i] <= β[i = 1:n_Period, 1:n_var] <= u_β[i])

@NLexpression(
    surv_fit,
    Xβ[i = 1:n_fit], β0[data_fit.Period[i]] 
                        + β[data_fit.Period[i],1] * data_fit.FinAid_Rate[i]
                       # + β[data_fit.Period[i],2] * data_fit.student_loan_rate[i]
                        + β[data_fit.Period[i],2] * data_fit.fed_efc_rate[i]
                       # + β[data_fit.Period[i],3] * data_fit.Financing_Started[i]
                        + β[data_fit.Period[i],3] * data_fit.Pell_Ind[i]
                        + β[data_fit.Period[i],4] * data_fit.home_distance_std[i]
                        + β[data_fit.Period[i],5] * data_fit.Admit_Honor_Ind[i]
                        + β[data_fit.Period[i],6] * data_fit.Diff_Major_Ind[i]
                        + β[data_fit.Period[i],7] * data_fit.Gender_Ind[i]
                        + β[data_fit.Period[i],8] * data_fit.Eth_ASIAN_Ind[i]
                        + β[data_fit.Period[i],9] * data_fit.Eth_BLACK_Ind[i]
                        + β[data_fit.Period[i],10] * data_fit.Eth_HISPA_Ind[i]
                        + β[data_fit.Period[i],11] * data_fit.Eth_WHITE_Ind[i]
                        + β[data_fit.Period[i],12] * data_fit.Eth_Multi_Ind[i]
                        + β[data_fit.Period[i],13] * data_fit.Pros_Event_Ind[i]
                        + β[data_fit.Period[i],14] * data_fit.CampusTour_Ever_Ind[i]
                        + β[data_fit.Period[i],15] * data_fit.DecisionDay_Ever_Ind[i]
                        + β[data_fit.Period[i],16] * data_fit.Delay_Review_Ind[i]
                        #+ β[data_fit.Period[i],17] * data_fit.Loan_Ind[i]
                        

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

w = 1.05
#r = 0.00001
r = 0
@NLobjective(
    surv_fit,
        Min,
        -sum( w*data_fit.Dpst_Ind[i]*Xβ[i] - exp(Xβ[i])*data_fit.Period_length[i]  for i = 1:n_fit) / n_fit
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
    #, β_inst_grant = value.(β)[:,1]
    , β_FinAid = value.(β)[:,1]
    #, β_loan = value.(β)[:,2]
    , β_fed_efc = value.(β)[:,2]
    #, β_Financing = value.(β)[:,3]
    , β_Pell = value.(β)[:,3]
    , β_home = value.(β)[:,4]
    , β_Honor = value.(β)[:,5]
    , β_Diff_Major = value.(β)[:,6]
    , β_Gender = value.(β)[:,7]
    , β_ASIAN = value.(β)[:,8]
    , β_BLACK = value.(β)[:,9]
    , β_HISPA = value.(β)[:,10]
    , β_WHITE = value.(β)[:,11]
    , β_Multi = value.(β)[:,12]
    , β_Pros_Event = value.(β)[:,13]
    , β_CampusTour = value.(β)[:,14]
    , β_DecisionDay = value.(β)[:,15]
    , β_Delay_Review = value.(β)[:,16]
    #, β_loan = value.(β)[:,17]

)

print(para_df)


β0_fit = value.(β0)
β_fit = value.(β)

#=
Xβ[i = 1:n_fit], β0[data_fit.Period[i]] 
                        + β[data_fit.Period[i],1] * data_fit.FinAid_Rate[i]
                        + β[data_fit.Period[i],2] * data_fit.fed_efc_rate[i]
                        + β[data_fit.Period[i],3] * data_fit.Pell_Ind[i]
                        + β[data_fit.Period[i],4] * data_fit.home_distance_std[i]
                        + β[data_fit.Period[i],5] * data_fit.Admit_Honor_Ind[i]
                        + β[data_fit.Period[i],6] * data_fit.Diff_Major_Ind[i]
                        + β[data_fit.Period[i],7] * data_fit.Gender_Ind[i]
                        + β[data_fit.Period[i],8] * data_fit.Eth_ASIAN_Ind[i]
                        + β[data_fit.Period[i],9] * data_fit.Eth_BLACK_Ind[i]
                        + β[data_fit.Period[i],10] * data_fit.Eth_HISPA_Ind[i]
                        + β[data_fit.Period[i],11] * data_fit.Eth_WHITE_Ind[i]
                        + β[data_fit.Period[i],12] * data_fit.Eth_Multi_Ind[i]
                        + β[data_fit.Period[i],13] * data_fit.Pros_Event_Ind[i]
                        + β[data_fit.Period[i],14] * data_fit.CampusTour_Ever_Ind[i]
                        + β[data_fit.Period[i],15] * data_fit.DecisionDay_Ever_Ind[i]
                        + β[data_fit.Period[i],16] * data_fit.Delay_Review_Ind[i]
=#

Xβ_fit = (β0_fit[data_fit.Period] 
.+ β_fit[data_fit.Period,1] .* data_fit.FinAid_Rate
.+ β_fit[data_fit.Period,2] .* data_fit.fed_efc_rate
.+ β_fit[data_fit.Period,3] .* data_fit.Pell_Ind
.+ β_fit[data_fit.Period,4] .* data_fit.home_distance_std
.+ β_fit[data_fit.Period,5] .* data_fit.Admit_Honor_Ind
.+ β_fit[data_fit.Period,6] .* data_fit.Diff_Major_Ind
.+ β_fit[data_fit.Period,7] .* data_fit.Gender_Ind
.+ β_fit[data_fit.Period,8] .* data_fit.Eth_ASIAN_Ind
.+ β_fit[data_fit.Period,9] .* data_fit.Eth_BLACK_Ind
.+ β_fit[data_fit.Period,10] .* data_fit.Eth_HISPA_Ind
.+ β_fit[data_fit.Period,11] .* data_fit.Eth_WHITE_Ind
.+ β_fit[data_fit.Period,12] .* data_fit.Eth_Multi_Ind
.+ β_fit[data_fit.Period,13] .* data_fit.Pros_Event_Ind
.+ β_fit[data_fit.Period,14] .* data_fit.CampusTour_Ever_Ind
.+ β_fit[data_fit.Period,15] .* data_fit.DecisionDay_Ever_Ind
.+ β_fit[data_fit.Period,16] .* data_fit.Delay_Review_Ind

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

### Baseline
β_test = 1.6
θ_test = 0.5
1-exp(exp(β_test)log(1-θ_test))

### Pell
β_test = 0.56
θ_test = 0.15
1-exp(exp(β_test)log(1-θ_test))

β_test = 0.26
θ_test = 0.15
1-exp(exp(β_test)log(1-θ_test))

β_test = 0.5
θ_test = 0.1
1-exp(exp(β_test)log(1-θ_test))

β_test = 1
θ_test = 0.1
x_diff = 0.1
1-exp(exp(β_test*x_diff)log(1-θ_test))

θ_2= 0.159
log(-log(1-θ_2))
θ_1 = 0.1
log(-log(1-θ_1))
log(-log(1-θ_2)) - log(-log(1-θ_1))

θ_2= 0.11
log(-log(1-θ_2))
θ_1 = 0.1
log(-log(1-θ_1))
log(-log(1-θ_2)) - log(-log(1-θ_1))

θ_2= 0.2
log(-log(1-θ_2))
θ_1 = 0.1
log(-log(1-θ_1))
log(-log(1-θ_2)) - log(-log(1-θ_1))

θ_2= 0.105
log(-log(1-θ_2))
θ_1 = 0.1
log(-log(1-θ_1))
log(-log(1-θ_2)) - log(-log(1-θ_1))