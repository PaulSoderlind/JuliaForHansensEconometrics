using Printf, LinearAlgebra, DelimitedFiles, Statistics, ReadStatTables

include("jlFiles/printmat.jl")
include("jlFiles/UtilityFunctions.jl")

x  = readstat("Data/DDK2011.dta")   #read data using the ReadStatTables.jl package

X  = excise( hcat(x.wordscore,x.sentscore,x.letterscore,x.spellscore,x.additions_score,
                 x.substractions_score,x.multiplications_score) )
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
