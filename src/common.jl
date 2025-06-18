const var_dists = Dict{Symbol, Distribution}()
# @variables_with_dist x ~ Normal(0, 1)
# var_dists[:x] gives Normal(0,1)
macro variables_with_dist(expr)
    # Expect something like :(x ~ Normal(0, 1))
    if !(expr isa Expr && expr.head == :call && expr.args[1] === :~)
        error("@variables_with_dist expects syntax like `x ~ Distribution(...)`")
    end
    
    var_sym = expr.args[2]
    dist_expr = expr.args[3]
    
    # Generate the variable as a Symbolics.Num
    var_def = :(Symbolics.@variables $var_sym)
    
    # Evaluate the distribution expression to get the Distribution object
    dist_val = eval(dist_expr)
    
    # Store mapping in var_dists dictionary
    push!(var_dists, var_sym => dist_val)
    
    # Return the variable symbol for use
    return var_def
end