using Distributions, Symbolics, Polynomials, Expectations

function moment(d::Normal, n::Int)
    @variables σ μ x
    f = taylor_coeff(taylor(exp(μ*x + σ*x^2/2), x, 0, 1:n), x, n)
    return substitute(f, Dict(μ => mean(d), σ => var(d))) * factorial(n)
end

# computes expectated values of iid Normal of polynomials
function linearexpectation(f::AbstractUnivariatePolynomial, d::Normal)
    m = 0
    n = length(f) - 1
    @variables σ μ x
    t = taylor(exp(μ*x + σ*x^2/2), x, 0, 1:n)
    for i in eachindex(f)
        if i == 0
            m += f[i]
        else
            ex = substitute(taylor_coeff(t, x, i), Dict(μ => mean(d), σ => var(d))) * factorial(i)
            m += f[i]*ex
        end
    end
    return m
end

function linearexpectation(f::Num, d::Normal)
    @assert length(Symbolics.get_variables(f)) == 1 "Linear expression of IIDs should only have one variable"
    @variables x
    f = substitute(f, first(Symbolics.get_variables(f)) => x)
    symbols = Symbolics.arguments(Symbolics.value(f))
    m = 0
    for s in symbols
        f = build_function(s, x, expression = Val{false})
        temp = expectation(x -> f(x), d)
        (isinf(temp) || isnan(temp)) && @warn("$(s) is not finite")
        m += temp
    end
    return m
end

function linearexpectation(f::Num, d::Vector{<:ContinuousUnivariateDistribution})
    n = length(d)
    @assert length(Symbolics.get_variables(f)) == n "Number of variables ≠ number of distributions"
    symbols = Symbolics.arguments(Symbolics.value(f))
    m = 0
    for (i,s) in enumerate(symbols)
        f = build_function(s, first(Symbolics.get_variables(s)), expression = Val{false})
        temp = expectation(x -> f(x), d[i])
        (isinf(temp) || isnan(temp)) && @warn("E[$(s)] wrt $(d[i]) is not finite")
        m += temp
    end
    return m
end

