using Turing, MCMCChains, ReverseDiff
using Distributions
using FillArrays

include("utils.jl")

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

@model function BIDRA(xs, ys)    
    HDR_B ~ Beta(0.4, 0.7)
    HDR ~ SkewNormal(HDR_B, 10, 0)
    #λₖ = [.7, .3]
    #HDR ~ MixtureModel([SkewNormal(100*HDR_B, 5, 10), Normal(100., 10.)], λₖ)
    #HDR_SN ~ SkewNormal(0, 2, 10)
    #HDR_N ~ Normal(100, 10)
    
    LDR ~ Normal(100,10)
    ic50 ~ Normal(0,10)
    slope ~ LogNormal(0.5, 1)
    
    σ ~ LogNormal(1, 1)
    
    
    #k ~ Categorical(λₖ)
    #if k == 1 
    #    HDR = HDR_SN + (HDR_B * 100) 
    #else
    #    HDR = HDR_N
    #end

    for i in 1:length(xs)
        f = llogistic([LDR, HDR, ic50, slope])
        ys[i] ~ Normal(f(xs[i]), σ)
    end
end
