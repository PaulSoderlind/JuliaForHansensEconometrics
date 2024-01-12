#####################################
### This file executes the
### empirical work reported
### in Chapter 8
#####################################

# Load the data & create variables

using DelimitedFiles, Printf, LinearAlgebra
include("jlFiles/printmat.jl")
include("jlFiles/ExtraFunctions.jl")      #some extra utility functions


#   Load the data & create subsamples
data = readdlm("Data/MRW1992.txt",skipstart=1)
replace!(data,"NA"=>NaN)             #using NaN instead of NA
data = data[:,2:end]                 #cut the first column

N          = convert.(Int,data[:,1])
(Y60,Y85,pop_growth,invest,school) = [convert.(Float64,data[:,i]) for i in [4,5,7,8,9]]

lndY  = log.(Y85)-log.(Y60)
lnY60 = log.(Y60)
lnI   = log.(invest/100)
lnG   = log.(pop_growth/100.0.+0.05)
lnS   = log.(school/100)

xx = [lnY60 lnI lnG lnS ones(size(lndY,1))]
x  = xx[N.==1,:]
y  = lndY[N.==1]
(n,k) = size(x)

# Unrestricted regression
invx     = inv(x'*x)
beta_ols = x\y
e_ols    = y - x*beta_ols
xe_ols   = x.*e_ols
V_ols    = (n/(n-k))*invx*(xe_ols'*xe_ols)*invx
se_ols   = sqrt.(diag(V_ols))
βse_ols  = MatrixInterleave(beta_ols,se_ols)


# Constrained regression
R        = [0 1 1 1 0]'
iR       = invx*R*inv(R'*invx*R)*R'
beta_cls = beta_ols - iR*beta_ols
e_cls    = y - x*beta_cls
xe_cls   = x.*e_cls
V_tilde  = (n/(n-k+1))*invx*(xe_cls'*xe_cls)*invx
V_cls    = V_tilde - iR*V_tilde - V_tilde*iR' + iR*V_tilde*iR'
se_cls   = sqrt.(diag(V_cls))
βse_cls  = MatrixInterleave(beta_cls,se_cls)


# (3) Efficient minimum distance
beta_emd = beta_ols - V_ols*R*inv(R'*V_ols*R)*R'*beta_ols
e_emd    = y - x*beta_emd
xe_emd   = x.*e_emd
V2       = (n/(n-k+1))*invx*(xe_emd'*xe_emd)*invx
V_emd    = V2 - V2*R*inv(R'*V2*R)*R'*V2
se_emd   = sqrt.(diag(V_emd))
βse_emd  = MatrixInterleave(beta_emd,se_emd)


xNames   = ["GDP60","I/GDP","n+g+δ","School","c"]
rowNames = MatrixInterleave(xNames,fill("",k))
println("Table 8.1")
printmat(βse_ols,βse_cls,βse_emd;rowNames,colNames=["ols","cls","emd"],prec=2)
