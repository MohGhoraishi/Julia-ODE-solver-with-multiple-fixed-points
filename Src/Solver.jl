
using Plots

function EulerODE(secondDerX::Function, firstDerX0::Number, resolution::Int, xStart::Number, tStart::Number, tEnd::Number)

    h = 10^(-resolution)
    Steps = Int((tEnd - tStart) / h)
    X = Float64[xStart]
    D1X = Float64[firstDerX0]
    T = Float64[tStart]
    for i in 1:Steps
        push!(T, T[i] + h)
        push!(D1X, D1X[i] + secondDerX(X[i], D1X[i], T[i]) * h)
        push!(X, X[i] + D1X[i] * h)
    end


end


    