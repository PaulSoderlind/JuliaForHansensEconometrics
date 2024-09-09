#####################################
### This file executes the
### empirical work reported
### in Chapter 4
###
### Uses data files cps09mar.txt; DDK2011.dta
#####################################


using DelimitedFiles, Printf, LinearAlgebra, Statistics, Distributions, ReadStatTables
include("jlFiles/printmat.jl")


#diary("chapter4matlab.log")

#   Load the data & create subsamples
dat = readdlm("Data/cps09mar.txt")
experience = dat[:,1] - dat[:,4] .- 6
mbf = (dat[:,11].==2) .& (dat[:,12].<=2) .& (dat[:,2].==1) .& (experience.==12)
sam = (dat[:,11].==4) .& (dat[:,12].==7) .& (dat[:,2].==0)

dat1=dat[mbf,:]
dat2=dat[sam,:]

#   Table 4.1
y = log.(dat1[:,5]./(dat1[:,6].*dat1[:,7]))
x = [dat1[:,4] ones(size(dat1,1))]
beta = x\y

e         = y - x*beta
leverage = sum((x.*(x*inv(x'*x))),dims=2)

(n,k) = size(x)
a = n/(n-k)
sig2 = (e'*e)/(n-k)
u1 = x.*e
u2 = x.*(e./sqrt.(1.0.-leverage))
u3 = x.*(e./(1.0.-leverage))
xx = inv(x'*x)
v0 = xx*sig2
v1 = xx*(u1'*u1)*xx
v1a = a*xx*(u1'*u1)*xx
v2 = xx*(u2'*u2)*xx
v3 = xx*(u3'*u3)*xx
s0  = sqrt.(diag(v0))  # Homoskedastic formula
s1  = sqrt.(diag(v1))  # White formula
s1a = sqrt.(diag(v1a)) # HC1 formula
s2  = sqrt.(diag(v2))  # HC2 formula
s3  = sqrt.(diag(v3))  # HC3 formula

println(" Table 4.1 [Regression [3.12]] \n")
println(" Coefficient estimates \n")
printmat(beta')

println(" \n Regression [3.13] \n")
println(" Coefficient estimates \n")
printmat(beta')
colNames = ["edu","c"]
rowNames = ["Homoskedastic se","White se","HC1 se","HC2 se","HC3 se"]
printmat([s0';s1';s1a';s2';s3'];rowNames,colNames)


# Equation 3.13
y          = log.(dat2[:,5]./(dat2[:,6].*dat2[:,7]))
experience = dat2[:,1] - dat2[:,4] .- 6
exp2       = (experience.^2)/100
x          = [dat2[:,4] experience exp2 ones(size(dat2,1))]
beta = x\y

e = y - x*beta
leverage = sum((x.*(x*inv(x'*x))),dims=2)

(n,k) = size(x)
a = n/(n-k)
sig2 = (e'*e)/(n-k)
u1 = x.*e
u2 = x.*(e./sqrt.(1.0.-leverage))
u3 = x.*(e./(1.0.-leverage))
xx = inv(x'*x)
v0  = xx*sig2
v1  = xx*(u1'*u1)*xx
v1a = a*xx*(u1'*u1)*xx
v2  = xx*(u2'*u2)*xx
v3  = xx*(u3'*u3)*xx
s0  = sqrt.(diag(v0)) # Homoskedastic formula
s1  = sqrt.(diag(v1)) # White formula
s1a = sqrt.(diag(v1a)) # HC1 formula
s2  = sqrt.(diag(v2)) # HC2 formula
s3  = sqrt.(diag(v3)) # HC3 formula

println(" \n Regression [3.13] \n")
println(" Coefficient estimates \n")
printmat(beta')
colNames = ["x1","x2","x3","c"]
rowNames = ["Homoskedastic se","White se","HC1 se","HC2 se","HC3 se"]
printmat([s0';s1';s1a';s2';s3'];rowNames,colNames)


# Table 4.2.
edu12     = dat[:,4] .> 11
dat3      = dat[edu12,:]
marriedF  = (dat3[:,12].<=3) .& (dat3[:,2].==1)
marriedM  = (dat3[:,12].<=3) .& (dat3[:,2].==0)
unionF    = (dat3[:,8].==1) .& (dat3[:,2].==1)
unionM    = (dat3[:,8].==1) .& (dat3[:,2].==0)
fmarriedF = (dat3[:,12].<=6) .& (dat3[:,12].>3) .& (dat3[:,2].==1)
fmarriedM = (dat3[:,12].<=6) .& (dat3[:,12].>3) .& (dat3[:,2].==0)
black = dat3[:,11] .==2
american_indian = dat3[:,11].==3
asian = dat3[:,11].==4
mixed = dat3[:,11].>=6

y = log.(dat3[:,5]./(dat3[:,6].*dat3[:,7]))
experience = dat3[:,1] - dat3[:,4] .- 6
exp2 = (experience.^2)/100
x = hcat( dat3[:,4],experience,exp2,dat3[:,2],
          unionF,unionM,marriedF,marriedM,fmarriedF,fmarriedM,
          dat3[:,3],black,american_indian,asian,mixed,ones(size(dat3,1)) )
beta = x\y
e    = y - x*beta
leverage = sum((x.*(x*inv(x'*x))),dims=2)
(n,k) = size(x)
u2 = x.*(e./sqrt.(1.0.-leverage))
xx = inv(x'*x)
v2 = xx*(u2'*u2)*xx
s2 = sqrt.(diag(v2)) # HC2 formula

println("\n Table 4.2: Log Wage Regression with HC2 se\n")
printmat([beta s2];colNames=["Coef","HC2 se"],rowNames="x")

println("\n EXTRA Table 4.2 (extension): p-value discussion from Section 9.8\n")
tstat = beta./s2
pval  = 2*ccdf.(Normal(0,1),abs.(tstat))    # ccdf(x) = 1-cdf(x)
printmat([beta s2 tstat pval];colNames=["Coef","HC2 se","t-stat","p value"],rowNames="x")

println("\n EXTRA Table 4.2 (extension): Wald test discussion from Section 9.10\n")
R = zeros(2,k)
R[1,5] = 1            #pick out female union member
R[2,6] = 1            #pick out male union member
q      = size(R,1)
printmat(R*beta)
W    = (R*beta-zeros(2))'*inv(R*v2*R')*(R*beta-zeros(2))
F    = W/q
pvalW = ccdf(Chisq(q),W)
pvalF = ccdf(FDist(q,n-k),F)
printmat([W pvalW;F pvalF];colNames=["test stat","p-value"],rowNames=["Wald","F"])

# DDK [2011]
# Load the data & create variables
data = readstat("Data/DDK2011.dta")   #read data using the ReadStatTables.jl package
schoolid   = convert.(Int,data.schoolid)
tracking   = convert.(Int,data.tracking)
totalscore = convert.(Float64,data.totalscore)

y = (totalscore .- mean(totalscore))./std(totalscore)
x = [tracking ones(size(y,1))]
(n,k) = size(x)
xx   = x'*x
invx = inv(xx)
beta = x\y
e    = y - x*beta

# Clustered robust standard error()
schools = unique(schoolid)
G = length(schools)
cluster_sums = zeros(G,k)
for g = 1:G
    local vg
    vg = schoolid .== schools[g]
    cluster_sums[g,:] = sum(x[vg,:].*e[vg],dims=1)
end
omega = cluster_sums'*cluster_sums
scale = G/(G-1)*(n-1)/(n-k)
V_clustered = scale*invx*omega*invx
se_clustered = sqrt.(diag(V_clustered))
println("\n Regression [4.57] of test scores on the tracking dummy \n")
println(" Coefficient estimates \n")
printmat(beta')
println(" Clustered ste \n")
printmat(se_clustered')
#diary off
