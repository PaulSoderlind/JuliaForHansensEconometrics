#####################################
### This file executes the
### empirical work reported
### in Chapter 3
#####################################

#diary("chapter3matlab.log")

using DelimitedFiles, Printf
include("jlFiles/printmat.jl")


#   Load the data & create subsamples
dat = readdlm("Data/cps09mar.txt")

experience = dat[:,1] - dat[:,4] .- 6
mbf  = (dat[:,11].==2) .& (dat[:,12].<=2) .& (dat[:,2].==1) .& (experience.==12)
sam  = (dat[:,11].==4) .& (dat[:,12].==7) .& (dat[:,2].==0)
dat1 = dat[mbf,:]

#   First regression
y    = log.(dat1[:,5]./(dat1[:,6].*dat1[:,7]))
x    = [dat1[:,4] ones(size(dat1,1))]
beta = x\y
println("\n Regression [3.12] \n")
printmat(beta)

#   Second regression
dat2       = dat[sam,:]
y          = log.(dat2[:,5]./(dat2[:,6].*dat2[:,7]))
experience = dat2[:,1] - dat2[:,4] .- 6
exp2       = (experience.^2)/100
x          = [dat2[:,4] experience exp2 ones(size(dat2,1))]
beta = x\y
println("\n Regression [3.13] \n")
printmat(beta)

# Create leverage & influence
e        = y - x*beta
xx       = x'*x
xxi      = inv(xx)
leverage = sum((x.*(x*xxi)),dims=2)
d        = leverage.*e./(1.0.-leverage)
influence = maximum(abs.(d))
println("\n The influence of log wage regression [3.13], cf. [3.48] \n")
printmat(influence)

#  Regression with the restricted sample
index45 = experience .< 45
x_r     = x[index45,:]
y_r     = y[index45]
beta_r  = x_r\y_r
println("\n Regression [3.49] \n")
printmat(beta_r)

#diary off
