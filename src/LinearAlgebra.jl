Base.@pure function intlog(b, x)
    i = 0
    b = oftype(x, b)
    while true
        b^i == x && return i
        b^i > x && return i - 1
        i += 1
    end
end

function norm(F::Type{<:AbstractGaloisField}, x::AbstractGaloisField)
    q = length(F)
    qn = length(typeof(x))
    n = intlog(q, qn)
    e = div(qn - 1, q - 1)
    F(x^e)
end

function tr(F::Type{<:AbstractGaloisField}, x::AbstractGaloisField)
    qn = length(typeof(x))
    q = oftype(qn, length(F))
    n = intlog(q, qn)
    res = zero(x)
    xqi = x  # loop invariant: xqi = x^(q^i)
    for i in 0 : n - 1
        res += xqi
        xqi = xqi^q
    end
    F(res)
end

norm(x::AbstractGaloisField) = norm(GaloisField(char(x)), x)
tr(x::AbstractGaloisField)   = tr(GaloisField(char(x)), x)
