using Printf, LinearAlgebra, DelimitedFiles, Statistics

include("jlFiles/printmat.jl")
include("jlFiles/UtilityFunctions.jl")


(x,header) = readdlm("Data/DKK2011.csv",',',header=true)   #load data
replace!(x,"NA"=>NaN)
x      = x[:,2:end]
header = header[2:end]

XX = PutDataInNT(x,header)                           #put data in Named Tuple
X  = excise(hcat(XX.wordscore,XX.sentscore,XX.letterscore,XX.spellscore,XX.additions_score,XX.substractions_score,XX.multiplications_score))
##----------------------------------------------------------

"""
    PrinComp(X,corQ=true)

Principal component analysis.

### Input
- `X::Matrix`:           Txn matrix of data
- `corQ::Bool`:          true: pc based on correlation matrix, otherwise covariance matrix


"""
function PrinComp(X,corQ=true)
  C       = corQ ? cor(X) : cov(X)
  F       = svd(C)
  W       = F.U
  λ       = F.S
  relvar  = λ/sum(λ)
  XDemean = X .- mean(X,dims=1)
  pc      = XDemean*W
  return pc, relvar, W, λ
end
##----------------------------------------------------------

(pc, relvar, W, λ) = PrinComp(X)

rowNames = ["words","sentences","letters","spelling","additions","substractions","multiplications"]
println("\neigenvalues and relative proportion")
printmat(λ,relvar;colNames=["Eigenvalue","proportion"],prec=2,width=13)
println("eigenvectors")
printmat(W[:,1:2];colNames=["1st","2nd"],rowNames,prec=2)

println("Notice that the sign of an eigenvector can be flipped")
##----------------------------------------------------------
