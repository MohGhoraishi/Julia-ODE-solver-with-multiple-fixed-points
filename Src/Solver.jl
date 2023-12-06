
using Plots

function euler_ode(secondDerX::Function, firstDerX0::Number, resolution::Int, xStart::Number, tStart::Number, tEnd::Number)
    h = float(10)^(-resolution)
    Steps = Int((tEnd - tStart) / h)
    X = Float64[xStart]
    D1X = Float64[firstDerX0]
    T = Float64[tStart]
    for i in 1:Steps
        push!(T, T[i] + h)
        push!(D1X, D1X[i] + secondDerX(X[i], D1X[i], T[i]) * h)
        push!(X, X[i] + D1X[i] * h)
    end
    return X, T
end

function euler_ode_end_only(secondDerX::Function, firstDerX0::Number, resolution::Int, xStart::Number, tStart::Number, tEnd::Number)
    h = float(10)^(-resolution)
    Steps = Int((tEnd - tStart) / h)
    X = xStart
    D1X = firstDerX0
    T = tStart
    for i in 1:Steps
        D1X2 = D1X + secondDerX(X, D1X, T) * h
        X2 = X + D1X * h
        T += h
        X = X2
        D1X = D1X2
    end
    return X
end

function fit_euler_ode(x1::Number, y1::Number, x2::Number, y2::Number, secondDerX::Function, resolution::Int)
    firstDer = 0
    shift = pi / 4
    recursive_euler_fitter(x1, y1, x2, y2, secondDerX, resolution, firstDer, shift)
end

function recursive_euler_fitter(x1::Number, y1::Number, x2::Number, y2::Number, secondDerX::Function, resolution::Int, currentFirstDer::Number, currentShift::Number)
    firstDerUp = currentFirstDer + tan(currentShift)
    firstDerDown = currentFirstDer - tan(currentShift)
    deltaUp = euler_endpoint_delta(x1, y1, x2, y2, secondDerX, firstDerUp, resolution)
    deltaDown = euler_endpoint_delta(x1, y1, x2, y2, secondDerX, firstDerDown, resolution)
    if deltaUp > deltaDown
        if abs(deltaDown) <= (float(10)^(-resolution))
            return firstDerDown
        else
            recursive_euler_fitter(x1, y1, x2, y2, secondDerX, resolution, firstDerDown, currentShift / 2)
        end
    elseif deltaUp < deltaDown
        if abs(deltaUp) <= (float(10)^(-resolution))
            return firstDerUp
        else
            recursive_euler_fitter(x1, y1, x2, y2, secondDerX, resolution, firstDerUp, currentShift / 2)
        end
    end
end

function euler_endpoint_delta(x1::Number, y1::Number, x2::Number, y2::Number, secondDerX::Function, firstDerX0::Number, resolution::Int)
    endPoint = euler_ode_end_only(secondDerX, firstDerX0, resolution, y1, x1, x2)
    return (endPoint - y2)
end

f(x, xdot, t) = - x

a = fit_euler_ode(0, 1, 1.5, 0.071, f, 4)