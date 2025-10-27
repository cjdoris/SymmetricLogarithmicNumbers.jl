module SymmetricLogarithmicNumbers

export SymLogarithmic, SymLogarithmicFloat16, SymLogarithmicFloat32, SymLogarithmicFloat64, symexp, symlog

issmall(x) = abs(x) < one(x)

function h(x)
    ax = abs(x)
    ix = inv(ax)
    r1 = ax - one(ax)
    r2 = one(ix) - ix
    r1, r2 = promote(r1, r2)
    return ifelse(ax < one(ax), r2, r1)
end

function invh(x)
    r1 = x + one(x)
    r2 = inv(one(x) - x)
    r1, r2 = promote(r1, r2)
    return ifelse(signbit(x), r2, r1)
end

"""
    symexp(x)

Compute the `symexp` function at `x`.

It is defined as:
- `exp(x - 1)` if `x ≥ 1`
- `exp(1 - 1/x)` if `x ≥ 0`
- `-exp(1 + 1/x)` if `x ≥ -1`
- `-exp(-x - 1)` otherwise
    
It has these nice arithmetic properties:
- It is an increasing, continuous, differentiable, invertible function on the real numbers.
- It grows exponentially for large `x` and shrinks exponentially for small `x`.
- `symexp(x) = x` for `x` in `-Inf`, `-1`, `0`, `1`, `Inf`.
- `-symexp(x) = symexp(-x)`
- `1/symexp(x) = symexp(1/x)`
- Its derivative is `1` at `±1`, and `0` at `0`.

The inverse function is [`symlog`](@ref).
"""
symexp(x) = copysign(exp(h(x)), x)

"""
    symlog(x)

The inverse function of [`symexp`](@ref).
"""
symlog(x) = copysign(invh(log(abs(x))), x)


### CONSTRUCTION / CONVERSION / PROMOTION

"""
    SymLogarithmic{[T]}(x)

Convert `x` to a `SymLogarithmic` number, which is represented by `symlog(x)`.

If you already know `ix = symlog(x)`, then you may call `symexp(SymLogarithmic, ix)` instead.

If you already know `lx = log(x)`, then you may call `exp(SymLogarithmic, ix)` instead.
"""
struct SymLogarithmic{T<:Real} <: Real
    symlog::T
    global symexp(::Type{SymLogarithmic{T}}, x::T) where {T<:Real} = new{T}(x)
end

const SymLogarithmicFloat16 = SymLogarithmic{Float16}
const SymLogarithmicFloat32 = SymLogarithmic{Float32}
const SymLogarithmicFloat64 = SymLogarithmic{Float64}

decon(x::SymLogarithmic) = (ix = symlog(x); (signbit(ix), abs(ix)))
recon(sb::Bool, ix::Real) = symexp(SymLogarithmic, ifelse(sb, -ix, +ix))

"""
    symexp(SymLogarithmic, x)

Compute `symexp(x)` but the result is stored exactly as a `SymLogarithmic`.
"""
symexp(::Type{SymLogarithmic{T}}, x::Real) where {T<:Real} = symexp(SymLogarithmic{T}, convert(T, x))
symexp(::Type{SymLogarithmic}, x::T) where {T<:Real} = symexp(SymLogarithmic{T}, x)

symlog(x::SymLogarithmic) = x.symlog
symexp(x::SymLogarithmic) = symexp(SymLogarithmic, symexp(symlog(x)))

Base.convert(::Type{SymLogarithmic{T}}, x::SymLogarithmic{T}) where {T} = x
Base.convert(::Type{SymLogarithmic{T}}, x::Real) where {T} = symexp(SymLogarithmic{T}, symlog(x))

Base.convert(::Type{SymLogarithmic}, x::SymLogarithmic) = x
Base.convert(::Type{SymLogarithmic}, x::Real) = symexp(SymLogarithmic, symlog(x))

SymLogarithmic(x::Real) = convert(SymLogarithmic, x)

SymLogarithmic{T}(x::Real) where {T} = convert(SymLogarithmic{T}, x)

Base.AbstractFloat(x::SymLogarithmic) = AbstractFloat(symexp(AbstractFloat(symlog(x))))

Base.BigFloat(x::SymLogarithmic) = BigFloat(symexp(BigFloat(symlog(x))))

Base.promote_rule(::Type{SymLogarithmic{T}}, ::Type{SymLogarithmic{S}}) where {T<:Real, S<:Real} = SymLogarithmic{promote_type(T, S)}
Base.promote_rule(::Type{SymLogarithmic{T}}, ::Type{S}) where {T<:Real, S<:Real} = promote_type(SymLogarithmic{T}, typeof(SymLogarithmic(zero(S))))


### COMPARISONS / PREDICATES

for op in [:signbit, :isfinite, :isinf, :isnan, :iszero, :isone]
    @eval Base.$op(x::SymLogarithmic) = $op(symlog(x))
end

for op in [:isequal, :isless, :cmp, :(==), :(!=), :(<), :(<=), :(>), :(>=)]
    @eval Base.$op(x::SymLogarithmic, y::SymLogarithmic) = $op(symlog(x), symlog(y))
end

for op in [:sign, :abs, :(+), :(-), :inv, :nextfloat, :prevfloat]
    @eval Base.$op(x::SymLogarithmic) = symexp(SymLogarithmic, $op(symlog(x)))
end

for op in [:zero, :one, :typemin, :typemax]
    @eval Base.$op(::Type{SymLogarithmic{T}}) where {T} = symexp(SymLogarithmic, $op(T))
end

for op in [:nextfloat, :prevfloat]
    @eval Base.$op(x::SymLogarithmic, n::Integer) = symexp(SymLogarithmic, $op(symlog(x), n))
end

Base.zero(::Type{SymLogarithmic}) = zero(SymLogarithmic{Int})

Base.one(::Type{SymLogarithmic}) = one(SymLogarithmic{Int})

# We hash the inner value, which for free means that hash(SymLogarithmic(x)) == hash(x) for x in
# [-Inf, -1, 0, 1, Inf], which are also the only values likely to overlap with other types.
Base.hash(x::SymLogarithmic, h::UInt) = hash(symlog(x), h)


### TYPES

Base.big(::Type{SymLogarithmic}) = SymLogarithmic{BigFloat}
Base.big(::Type{SymLogarithmic{T}}) where {T} = SymLogarithmic{big(T)}
Base.big(x::SymLogarithmic) = convert(big(typeof(x)), x)

Base.widen(::Type{SymLogarithmic{T}}) where {T} = SymLogarithmic{widen(T)}
Base.widen(x::SymLogarithmic) = convert(widen(typeof(x)), x)


### ARITHMETIC

function Base.:(*)(x::SymLogarithmic{T}, y::SymLogarithmic{T}) where {T}
    sx, ix = decon(x)
    sy, iy = decon(y)
    sr = xor(sx, sy)
    ir = invh(h(ix) + h(iy))
    return recon(sr, ir)
end

function Base.:(/)(x::SymLogarithmic{T}, y::SymLogarithmic{T}) where {T}
    sx, ix = decon(x)
    sy, iy = decon(y)
    sr = xor(sx, sy)
    ir = invh(h(ix) - h(iy))
    return recon(sr, ir)
end

function Base.:(^)(x::SymLogarithmic, p::Real) where {T}
    sx, ix = decon(x)
    if sx && !iszero(ix) && !isnan(ix)
        if isinteger(p)
            sr = isodd(p)
        else
            throw(DomainError(x, "Can only take integer powers of negative SymLogarithmic numbers."))
        end
    else
        sr = false
    end
    ir = invh(h(ix) * p)
    return recon(sr, ir)
end

function Base.:(^)(x::SymLogarithmic, p::Integer) where {T}
    sx, ix = decon(x)
    sr = sx && isodd(p)
    ir = invh(h(ix) * p)
    return recon(sr, ir)    
end

Base.:(+)(x::SymLogarithmic{T}, y::SymLogarithmic{T}) where {T} = _add(decon(x)..., decon(y)...)

function _add(sx, ix, sy, iy)
    se = xor(sx, sy)
    sr = se ? xor(ix > iy, sy) : sx
    a, b = minmax(ix, iy)
    ha = h(a)
    hb = h(b)
    c = ha - hb
    if se
        if c < oftype(c, -1)
            d = log(-expm1(c))
        else
            d = log1p(-exp(c))
        end
    else
        d = log1p(exp(c))
    end
    ir = invh(hb + d)
    return recon(sr, ir)    
end

Base.:(-)(x::SymLogarithmic{T}, y::SymLogarithmic{T}) where {T} = _sub(decon(x)..., decon(y)...)

_sub(sx, ix, sy, iy) = _add(sx, ix, !sy, iy)

function Base.log(x::SymLogarithmic)
    ix = symlog(x)
    if signbit(ix) && !iszero(ix) && !isnan(ix)
        throw(DomainError(x))
    elseif issmall(ix)
        return one(ix) - inv(ix)
    else
        return ix - one(ix)
    end
end

function Base.log2(x::SymLogarithmic)
    logx = log(x)
    return logx / log(oftype(logx, 2))
end

function Base.log10(x::SymLogarithmic)
    logx = log(x)
    return logx / log(oftype(logx, 10))
end

function Base.exp(x::SymLogarithmic)
    return symexp(SymLogarithmic, invh(symexp(symlog(x))))
end

function Base.exp(::Type{SymLogarithmic}, x::Real)
    if signbit(x)
        ir = inv(one(x) - x)
    else
        ir = x + inv(inv(one(x)))
    end
    return symexp(SymLogarithmic, ir)
end

function Base.exp(::Type{SymLogarithmic{T}}, x::Real) where {T}
    if signbit(x)
        ir = inv(one(T) - convert(T, x))
    else
        ir = convert(T, x) + one(T)
    end
    return symexp(SymLogarithmic{T}, ir)
end


## IO

function Base.show(io::IO, x::SymLogarithmic)
    print(io, "symexp(")
    if get(io, :typeinfo, Any) != typeof(x)
        show(io, typeof(x))
        print(io, ", ")
    end
    show(io, symlog(x))
    print(io, ")")
end

function Base.write(io::IO, x::SymLogarithmic)
    write(io, symlog(x))
end

function Base.read(io::IO, ::Type{SymLogarithmic{T}}) where {T}
    symexp(SymLogarithmic{T}, read(io, T))
end

end # module
