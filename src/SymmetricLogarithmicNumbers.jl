module SymmetricLogarithmicNumbers

export SymLog, SymLogF16, SymLogF32, SymLogF64, symexp, symlog

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
    SymLog{[T]}(x)

Convert `x` to a `SymLog` number, which is represented by `symlog(x)`.

If you already know `ix = symlog(x)`, then you may call `symexp(SymLog, ix)` instead.

If you already know `lx = log(x)`, then you may call `exp(SymLog, ix)` instead.
"""
struct SymLog{T<:Real} <: Real
    symlog::T
    global symexp(::Type{SymLog{T}}, x::T) where {T<:Real} = new{T}(x)
end

const SymLogF16 = SymLog{Float16}
const SymLogF32 = SymLog{Float32}
const SymLogF64 = SymLog{Float64}

decon(x::SymLog) = (ix = symlog(x); (signbit(ix), abs(ix)))
recon(sb::Bool, ix::Real) = symexp(SymLog, ifelse(sb, -ix, +ix))

"""
    symexp(SymLog, x)

Compute `symexp(x)` but the result is stored exactly as a `SymLog`.
"""
symexp(::Type{SymLog{T}}, x::Real) where {T<:Real} = symexp(SymLog{T}, convert(T, x))
symexp(::Type{SymLog}, x::T) where {T<:Real} = symexp(SymLog{T}, x)

symlog(x::SymLog) = x.symlog
symexp(x::SymLog) = symexp(SymLog, symexp(symlog(x)))

Base.convert(::Type{SymLog{T}}, x::SymLog{T}) where {T} = x
Base.convert(::Type{SymLog{T}}, x::Real) where {T} = symexp(SymLog{T}, symlog(x))

Base.convert(::Type{SymLog}, x::SymLog) = x
Base.convert(::Type{SymLog}, x::Real) = symexp(SymLog, symlog(x))

SymLog(x::Real) = convert(SymLog, x)

SymLog{T}(x::Real) where {T} = convert(SymLog{T}, x)

Base.AbstractFloat(x::SymLog) = AbstractFloat(symexp(AbstractFloat(symlog(x))))

Base.BigFloat(x::SymLog) = BigFloat(symexp(BigFloat(symlog(x))))

Base.promote_rule(::Type{SymLog{T}}, ::Type{SymLog{S}}) where {T<:Real, S<:Real} = SymLog{promote_type(T, S)}
Base.promote_rule(::Type{SymLog{T}}, ::Type{S}) where {T<:Real, S<:Real} = promote_type(SymLog{T}, typeof(SymLog(zero(S))))


### COMPARISONS / PREDICATES

for op in [:signbit, :isfinite, :isinf, :isnan, :iszero, :isone]
    @eval Base.$op(x::SymLog) = $op(symlog(x))
end

for op in [:isequal, :isless, :cmp, :(==), :(!=), :(<), :(<=), :(>), :(>=)]
    @eval Base.$op(x::SymLog, y::SymLog) = $op(symlog(x), symlog(y))
end

for op in [:sign, :abs, :(+), :(-), :inv, :nextfloat, :prevfloat]
    @eval Base.$op(x::SymLog) = symexp(SymLog, $op(symlog(x)))
end

for op in [:zero, :one, :typemin, :typemax]
    @eval Base.$op(::Type{SymLog{T}}) where {T} = symexp(SymLog, $op(T))
end

for op in [:nextfloat, :prevfloat]
    @eval Base.$op(x::SymLog, n::Integer) = symexp(SymLog, $op(symlog(x), n))
end

Base.zero(::Type{SymLog}) = zero(SymLog{Int})

Base.one(::Type{SymLog}) = one(SymLog{Int})

# We hash the inner value, which for free means that hash(SymLog(x)) == hash(x) for x in
# [-Inf, -1, 0, 1, Inf], which are also the only values likely to overlap with other types.
Base.hash(x::SymLog, h::UInt) = hash(symlog(x), h)


### TYPES

Base.big(::Type{SymLog}) = SymLog{BigFloat}
Base.big(::Type{SymLog{T}}) where {T} = SymLog{big(T)}
Base.big(x::SymLog) = convert(big(typeof(x)), x)

Base.widen(::Type{SymLog{T}}) where {T} = SymLog{widen(T)}
Base.widen(x::SymLog) = convert(widen(typeof(x)), x)


### ARITHMETIC

function Base.:(*)(x::SymLog{T}, y::SymLog{T}) where {T}
    sx, ix = decon(x)
    sy, iy = decon(y)
    sr = xor(sx, sy)
    ir = invh(h(ix) + h(iy))
    return recon(sr, ir)
end

function Base.:(/)(x::SymLog{T}, y::SymLog{T}) where {T}
    sx, ix = decon(x)
    sy, iy = decon(y)
    sr = xor(sx, sy)
    ir = invh(h(ix) - h(iy))
    return recon(sr, ir)
end

function Base.:(^)(x::SymLog, p::Real) where {T}
    sx, ix = decon(x)
    if sx && !iszero(ix) && !isnan(ix)
        if isinteger(p)
            sr = isodd(p)
        else
            throw(DomainError(x, "Can only take integer powers of negative SymLog numbers."))
        end
    else
        sr = false
    end
    ir = invh(h(ix) * p)
    return recon(sr, ir)
end

function Base.:(^)(x::SymLog, p::Integer) where {T}
    sx, ix = decon(x)
    sr = sx && isodd(p)
    ir = invh(h(ix) * p)
    return recon(sr, ir)    
end

Base.:(+)(x::SymLog{T}, y::SymLog{T}) where {T} = _add(decon(x)..., decon(y)...)

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

Base.:(-)(x::SymLog{T}, y::SymLog{T}) where {T} = _sub(decon(x)..., decon(y)...)

_sub(sx, ix, sy, iy) = _add(sx, ix, !sy, iy)

function Base.log(x::SymLog)
    ix = symlog(x)
    if signbit(ix) && !iszero(ix) && !isnan(ix)
        throw(DomainError(x))
    elseif issmall(ix)
        return one(ix) - inv(ix)
    else
        return ix - one(ix)
    end
end

function Base.log2(x::SymLog)
    logx = log(x)
    return logx / log(oftype(logx, 2))
end

function Base.log10(x::SymLog)
    logx = log(x)
    return logx / log(oftype(logx, 10))
end

function Base.exp(x::SymLog)
    return symexp(SymLog, invh(symexp(symlog(x))))
end

function Base.exp(::Type{SymLog}, x::Real)
    if signbit(x)
        ir = inv(one(x) - x)
    else
        ir = x + inv(inv(one(x)))
    end
    return symexp(SymLog, ir)
end

function Base.exp(::Type{SymLog{T}}, x::Real) where {T}
    if signbit(x)
        ir = inv(one(T) - convert(T, x))
    else
        ir = convert(T, x) + one(T)
    end
    return symexp(SymLog{T}, ir)
end


## IO

function Base.show(io::IO, x::SymLog)
    print(io, "symexp(")
    if get(io, :typeinfo, Any) != typeof(x)
        show(io, typeof(x))
        print(io, ", ")
    end
    show(io, symlog(x))
    print(io, ")")
end

function Base.write(io::IO, x::SymLog)
    write(io, symlog(x))
end

function Base.read(io::IO, ::Type{SymLog{T}}) where {T}
    symexp(SymLog{T}, read(io, T))
end

end # module
