using LinearAlgebra, DSP, Plots

function StochasticSTF(n,r::Float64=1.0,d::Int64=2)
    n1=Int(floor(n*r/(1+r)))
    n2=n-n1+1
    B1=BesselBridge(n1+1,d)
    B2=BesselBridge(n2+1,d)
    STF=conv(B1[2:n1+1],B2[1:n2])
    return STF/sum(STF)
end

function BesselBridge(n,d)
    x=BrownBridges(n,d)
    return mapslices(norm,x,dims=2)
end

function BrownBridges(n,d)
    x=vcat(zeros(1,d), cumsum(randn(n-1,d),dims=1))
    return x-collect(0:n-1)*x[n,:]'/(n-1)
end

n=1000
r=1.0
STF=StochasticSTF(n,r)
plot(STF)
