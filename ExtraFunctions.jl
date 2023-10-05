
##------------------------------------------------------------------------------
"""
    MatrixInterleaveMPs(x...)

Interleave the rows of matrices

"""
function MatrixInterleave(x...)
  all([size(xi)==size(first(x)) for xi in x]) || error("the input do not have the same size")
  xNum  = length(x)
  (T,K) = (size(first(x),1),size(first(x),2))
  y = reshape(collect(Iterators.flatten(zip(x...))),xNum*T,K)
  (K == 1) && (y = vec(y))          #to vector
  return y
end
##------------------------------------------------------------------------------


##---------------------------------------------------------------
"""
    quantileM(x,p)

Compute quantiles of each column in `x` where `p` is either a
number of an array with the same number of elements as columns in x.
In the latter case, we calculate quantile `p[i]` of `x[:,i]`.
"""
function quantileM(x::VecOrMat,p::Number)
  q = permutedims( [quantile(col,p) for col in eachcol(x)] )
  return q
end

function quantileM(x::VecOrMat,p::VecOrMat)
  (size(x,2) != length(p)) && error("size(x,2) != length(p)")
  q = permutedims( [quantile(view(x,:,i),p[i]) for i=1:length(p)] )
  return q
end
##---------------------------------------------------------------
