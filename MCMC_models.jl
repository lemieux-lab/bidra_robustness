using Turing, MCMCChains, ReverseDiff
using Distributions

@model function modelNormal(x)
    μ ~ Normal()
    σ ~ LogNormal()

    x .~ Normal.(μ, σ)
end


@model function modelNormalBeta(x)
    α ~ LogNormal(0,1)
    β ~ LogNormal(0,1)
    σ ~ LogNormal(0.1,0.1)
    

    N = length(x)
    Bᵢ ~ filldist(Beta(α, β), N)

    ## Likelihood
    ŷ = Bᵢ .* 100
    x .~ Normal.(ŷ, σ)
end
