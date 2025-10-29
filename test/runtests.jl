using Test
using SymmetricLogarithmicNumbers

const SLN = SymmetricLogarithmicNumbers

realvalue(x::SLN.SymLogarithmic) = symexp(symlog(x))

@testset "SymmetricLogarithmicNumbers" begin
    @testset "symexp ↔ symlog" begin
        real_inputs = [-100.0, -10.0, -2.0, -1.0, -0.5, -0.1, -1e-6, 0.0, 1e-6, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0]
        for x in real_inputs
            y = symlog(x)
            @test isapprox(symexp(y), x; atol=1e-12, rtol=1e-12)
        end

        log_inputs = [-5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
        for ix in log_inputs
            z = symexp(ix)
            @test isapprox(symlog(z), ix; atol=1e-12, rtol=1e-12)
        end
    end

    @testset "SymLogarithmic construction" begin
        xs = [-3.5, -1.0, -0.25, 0.0, 0.25, 1.0, 3.5]
        for x in xs
            hx = SymLogarithmic(x)
            @test hx isa SLN.SymLogarithmic{Float64}
            @test isapprox(realvalue(hx), x; atol=eps(x)*10, rtol=0)
        end

        @test SymLogarithmic{Float32}(2.0) isa SLN.SymLogarithmic{Float32}
        @test isapprox(realvalue(SymLogarithmic{Float32}(2.0)), 2f0; atol=eps(Float32)*10, rtol=0)
    end

    @testset "Arithmetic" begin
        pairs = [(-3.0, 0.5), (-1.5, -0.25), (0.75, 1.5), (2.0, -4.0)]
        for (x, y) in pairs
            hx, hy = SymLogarithmic(x), SymLogarithmic(y)
            @test isapprox(realvalue(hx + hy), x + y; atol=1e-12)
            @test isapprox(realvalue(hx - hy), x - y; atol=1e-12)
            @test isapprox(realvalue(hx * hy), x * y; atol=1e-12)
            @test isapprox(realvalue(hx / hy), x / y; atol=1e-12)
        end

        neg = SymLogarithmic(-3.0)
        @test isapprox(realvalue(neg ^ 2), 9.0; atol=1e-12)
        @test isapprox(realvalue(neg ^ 3), -27.0; atol=1e-12)
        @test_throws DomainError neg ^ 0.5

        pos = SymLogarithmic(4.0)
        @test isapprox(realvalue(pos ^ 0.5), 2.0; atol=1e-12)
    end

    @testset "Logarithms" begin
        pos = SymLogarithmic(5.0)
        @test isapprox(log(realvalue(pos)), log(pos); atol=1e-12)
        @test isapprox(log(realvalue(pos)) / log(2), log2(pos); atol=1e-12)
        @test isapprox(log(realvalue(pos)) / log(10), log10(pos); atol=1e-12)
        @test_throws DomainError log(SymLogarithmic(-0.5))
    end

    @testset "Zero and identity" begin
        zero_val = zero(SymLogarithmic{Float64})
        one_val = one(SymLogarithmic{Float64})
        @test realvalue(zero_val) == 0.0
        @test realvalue(one_val) == 1.0
        @test iszero(zero_val)
        @test isone(one_val)
    end
end
