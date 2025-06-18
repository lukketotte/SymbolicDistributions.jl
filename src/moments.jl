function moment(d::Normal, n::Int)
    @variables σ μ x
    f = taylor_coeff(taylor(exp(μ*x + σ*x^2/2), x, 0, 1:n), x, n)
    return substitute(f, Dict(μ => mean(d), σ => var(d))) * factorial(n)
end

moment(d::Chisq, n::Int) = 2^n * gamma(n+dof(d)/2) / gamma(dof(d)/2)