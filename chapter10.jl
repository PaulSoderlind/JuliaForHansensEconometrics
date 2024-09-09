#####################################
### This matlab file creates Tables 10.1; 10.2; & 10.3.
### It does the calculations reported in Section 10.8
### & generates the bootstrap draws for Figure 10.1
### (which are drawn in R in figure10_1.R)
###
### Uses data file cps09mar.xlsx
#####################################

#diary("chapter10matlab.log")

using DelimitedFiles, Printf, LinearAlgebra, Statistics, Random, Distributions
include("jlFiles/printmat.jl")
include("jlFiles/ExtraFunctions.jl")      #some extra utility functions


B = 10_000

#   Load the data & create subsamples
#dat = xlsread("cps09mar.xlsx")
dat = readdlm("Data/cps09mar.txt")

experience = dat[:,1] - dat[:,4] .- 6
mbf = (dat[:,11].==2) .& (dat[:,12].<=2) .& (dat[:,2].==1)

mbf12 = mbf .& (experience.==12)

dat1 = dat[mbf12,:]
dat2 = dat[mbf,:]

#   First regression
w   = dat1[:,5]./(dat1[:,6].*dat1[:,7])
y   = log.(w)
edu = dat1[:,4]

x     = [edu ones(size(dat1,1))]
(n,k) = size(x)
xx    = inv(x'*x)
beta  = x\y
e     = y - x*beta
sig2  = (e'*e)/n
betat = [beta;sig2]
g1    = [16,1]
g     = [g1;1/2]
mu    = exp(g1'*beta + sig2/2)
theta = [beta;sig2;mu]

#  Asymptotic s.e.
leverage = vec(sum((x.*(x*inv(x'*x)))',dims=1))
u   = x.*(e./sqrt.(1.0.-leverage))
v   = xx*(u'*u)*xx
s   = sqrt.(diag(v))
x0  = zeros(k)
xxt = [xx x0; x0' 1/n]
ut  = hcat(x.*e, abs2.(e).-sig2)
v   = xxt*(ut'*ut)*xxt
st  = mu*sqrt(g'*v*g)
e2  = abs2.(e)
sv  = sqrt(var(e2)/n)
stheta = [s' sv st]

# Jackknife
jack = zeros(n,4)
for i = 1:n
  local vi, yi, xi, betai, ei, sigi, mui
  vi    = deleteat!(collect(1:n),i)
  yi    = view(y,vi)
  xi    = view(x,vi,:)
  betai = xi\yi
  ei    = yi - xi*betai          #could speed up this, see comment below the loop
  sigi  = (ei'*ei)/(n-1)
  mui   = exp(g1'*betai+sigi/2)
  jack[i,:] = [betai;sigi;mui]
end
sn = ((n-1)^2)/n
seJ = sqrt.(sn*var(jack,dims=1))

# To speed up, do
# (ei,xb) = (fill(NaN,n),fill(NaN,n))    #before loop
# mul!(xb,x,b);  ei .= y .- xb           #instead of `ei = ` above

println("Table 10.1")
println("Jackknifed estimates")
printmat(jack)
println("Jackknife & asymptotic s.e.")
printmat([seJ;stheta])

# Bootstrap
Random.seed!(13)
u1 = rand(1:n,n)
  yb = view(y,u1)
  xb = view(x,u1,:)
  betahatb = xb\yb
  eb       = yb - xb*betahatb
  sigb     = (eb'*eb)/n
  mub      = exp(g1'*betahatb+sigb/2)
  theta1   = [betahatb;sigb;mub]

u2 = rand(1:n,n)
  yb = view(y,u2)
  xb = view(x,u2,:)
  betahatb = xb\yb
  eb       = yb - xb*betahatb
  sigb     = (eb'*eb)/n
  mub      = exp(g1'*betahatb+sigb/2)
  theta2   = [betahatb;sigb;mub]
  
println("Calculations in Section 10.8")
println("First two bootstrap sample indices")
printmat([u1 u2])
println("First two bootstrap estimate vectors")
printmat([theta1 theta2])

thetas = zeros(B,4)
tstats = zeros(B,4)
Random.seed!(13)
for bi = 1:B
  local u, yb, xb, xx, betahatb, eb, sigb, mub, leverage, v, s, x0, xxt, ut, v,
  st, e2, sv, sthetab
  u  = rand(1:n,n)
  yb = view(y,u)
  xb = view(x,u,:)
  xx = inv(xb'*xb)
  betahatb       = xb\yb
  thetas[bi,1:2] = betahatb
  eb   = yb - xb*betahatb
  sigb = (eb'*eb)/n
  mub  = exp(g1'*betahatb+sigb/2)
  thetas[bi,3] = sigb
  thetas[bi,4] = mub
  leverage = vec(sum((xb.*(xb*inv(xb'*xb)))',dims=1))
  u   = xb.*(eb./sqrt.(1.0.-leverage))
  v   = xx*(u'*u)*xx
  s   = sqrt.(diag(v))
  x0  = zeros(k)
  xxt = [xx x0; x0' 1/n]
  ut  = hcat(xb.*eb,abs2.(eb).-sigb)
  v   = xxt*(ut'*ut)*xxt
  st  = mub*sqrt(g'*v*g)
  e2  = abs2.(eb)
  sv  = sqrt(var(e2)/n)
  sthetab = [s; sv; st]
  tstats[bi,:] = (thetas[bi,:]-theta)./sthetab
end
# plenty of room for speeding up the loop by preallocating and filling

# Bootstrap Standard Errors
seboot = sqrt.(var(thetas,dims=1))

# Bootstrap Percentile Interval
q1 = quantileM(thetas,0.025)'
q2 = quantileM(thetas,0.975)'

# Bootstrap BC Percentile Interval
z0  = [quantile(Normal(0,1),mean(thetas[:,i] .<= theta[i])) for i=1:length(theta)]'
z1  = quantile(Normal(0,1),0.025)
z2  = quantile(Normal(0,1),0.975)
xa1 = cdf(Normal(0,1),z1.+2*z0)
xa2 = cdf(Normal(0,1),z2.+2*z0)
qa1 = quantileM(thetas,xa1)'
qa2 = quantileM(thetas,xa2)'

# Bootstrap BCa Percentile Interval
mj   = mean(jack,dims=1)
d    = sum(abs2,mj.-jack,dims=1).^(3/2)
a    = sum((mj.-jack).^3,dims=1)./d/6
xa1  = cdf(Normal(0,1),z0+(z1.+z0)./(1.0.-a.*(z1.+z0)))
xa2  = cdf(Normal(0,1),z0+(z2.+z0)./(1.0.-a.*(z2.+z0)))
qca1 = quantileM(thetas,xa1)'                #different from matlab file
qca2 = quantileM(thetas,xa2)'

qt1 = theta - quantileM(tstats,0.975)'.*stheta'
qt2 = theta - quantileM(tstats,0.025)'.*stheta'

# Bootstrap percentile-t
println("Table 10.2")
println("Estimates & s.e.")
coefNames = ["β₁","β₂","σ²","μ"]
printmat([theta';stheta;seJ;seboot];colNames=coefNames,
         rowNames=["estimate","asymp se","Jackknife se","bootstap se"])
println("95% Percentile Confidence Intervals")
printmat([q1 q2];colNames=["low","high"],rowNames=coefNames,prec=2)
println("Bias-Corrected Percentile 95% Confidence Intervals")
printmat([qa1 qa2];colNames=["low","high"],rowNames=coefNames,prec=2)
println("BCa Percentile 95% Confidence Intervals")
printmat([qca1 qca2];colNames=["low","high"],rowNames=coefNames,prec=2)
println("Percentile-t 95% Confidence Intervals")
printmat([qt1 qt2];colNames=["low","high"],rowNames=coefNames,prec=2)

# For Figure 10.1; save bootstrap replications
betas1 = thetas[:,1]
theta4 = thetas[:,4]
b1 = beta[1]
bootreps = [b1 mu; betas1 theta4]
fileID = open("Results/bootreps.txt","w")
    printmat(fileID,bootreps;width=12,prec=8)   #saves the bootreps matrix as is
close(fileID)
##------------------------------------------------------------

#   Second regression [Section 10.17]
y     = log.(dat2[:,5]./(dat2[:,6].*dat2[:,7]))
exp1  = dat2[:,1] - dat2[:,4] .- 6
exp2  = abs2.(exp1)/100
x     = [dat2[:,4] exp1 exp2 ones(size(dat2,1))]
(n,k) = size(x)
beta  = x\y
e     = y - x*beta
theta = -50*beta[2]/beta[3]

#  Asymptotic s.e.
leverage = vec(sum((x.*(x*inv(x'*x)))',dims=1))
u  = x.*(e./sqrt.(1.0.-leverage))
xx = inv(x'*x)
v  = xx*(u'*u)*xx
se = sqrt.(diag(v))
g  = [0; -50/beta[3]; -theta/beta[3]; 0]
s  = sqrt(g'*v*g)

# Jackknife
thetaJ = zeros(n)
for bi = 1:n
  local vi, yb, xb, betahatb
  vi = deleteat!(collect(1:n),bi)
  yb = view(y,vi)
  xb = view(x,vi,:)
  betahatb   = xb\yb
  thetaJ[bi] = -50*betahatb[2]/betahatb[3]
end
sn = ((n-1)^2)/n
sethetaJ = sqrt(sn*var(thetaJ))

# Bootstrap
thetas1 = zeros(B)
thetast = zeros(B)
tau = 25
for bi = 1:B
  local u, yb, xb, betahatb, thetab, thetat
  u  = rand(1:n,n)
  yb = view(y,u)
  xb = view(x,u,:)
  betahatb  = xb\yb
  thetab    = -50*betahatb[2]/betahatb[3]
  thetat    = thetab*(abs(thetab-theta) <= tau) + 
              (theta-tau)*(thetab < theta-tau) + (theta+tau)*(thetab > theta+tau)
  thetas1[bi] = thetab
  thetast[bi] = thetat
end
sethetaboot1 = sqrt(var(thetas1))
sethetaboott = sqrt(var(thetast))

# Repeat Bootstrap
thetas2 = zeros(B)
tau = 25
for bi = 1:B
  local u, yb, xb, betahatb, thetab
  u  = rand(1:n,n)
  yb = view(y,u)
  xb = view(x,u,:)
  betahatb = xb\yb
  thetab   = -50*betahatb[2]/betahatb[3]
  thetas2[bi] = thetab
end
sethetaboot2 = sqrt(var(thetas2))

println("Log[wage] equation")
printmat([beta se];rowNames=["edu","exp","exp2","c"],colNames=["coef","se"])

println("Tabled 10.3")
println("Experience with maximum wage")
println("Estimates & s.e.")
coefse = [theta, s, sethetaJ, sethetaboot1, sethetaboot2, sethetaboott]
rowNames=["est","asymp se","Jackknife se","bootstrap standard","bootstrap repeat","bootstrap trimmed"]
printmat(coefse;rowNames)

#diary off

