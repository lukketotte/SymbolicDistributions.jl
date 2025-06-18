using Distributions, Symbolics, Polynomials, TaylorSeries

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