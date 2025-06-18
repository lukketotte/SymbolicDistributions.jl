using Test, SymbolicDistributions, Distributions, Symbolics

@testset "moments" begin
    μ,σ = 1.2, 2
    # test first 8 moments 
    for i in 1:8
        m = moment(Normal(μ,σ), i)
        if i == 1
            @test m ≈ μ
        elseif i == 2
            @test m ≈ μ^2 + σ^2
        elseif i == 3
            @test m ≈ μ^3 + 3*μ*σ^2
        elseif i == 4
            @test m ≈ μ^4 + 6*μ^2*σ^2 +  3*σ^4
        elseif i == 5
            @test m ≈ μ^5 + 10*μ^3*σ^2 + 15*μ*σ^4
        elseif i == 6
            @test m ≈ μ^6 + 15*μ^4*σ^2 + 45*μ^2*σ^4 + 15*σ^6
        elseif i == 7
            @test m ≈ μ^7 + 21*μ^5*σ^2 + 105*μ^3*σ^4 + 105*μ*σ^6
        elseif i == 8
            @test m ≈ μ^8 + 28*μ^6*σ^2 + 210*μ^4*σ^4 + 420*μ^2*σ^6 + 105*σ^8
        end
    end
end

@testset "expectations" begin
    # odd functions 
    @variables x 
    linearexpectation(sinh(x), Normal(0,1)) ≈ 0
    linearexpectation(cos(x), Normal(0,1)) ≈ 0
    linearexpectation(sin(x), Normal(0,1)) ≈ 0
    
    

end