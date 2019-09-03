module Util

function hilo_mul(u::Int128, v::Int128)
    local u0::UInt128, v0::UInt128, w0::UInt128
    local u1::Int128, v1::Int128, w1::UInt128, w2::Int128, t::UInt128

    u0 = u % UInt64; u1 = u >> 64
    v0 = v % UInt64; v1 = v >> 64
    w0 = u0 * v0
    t = reinterpret(UInt128, u1) * v0 + (w0 >>> 64)
    w2 = reinterpret(Int128, t) >> 64
    w1 = u0 * reinterpret(UInt128, v1) + (t % UInt64)
    hi = u1 * v1 + w2 + (reinterpret(Int128, w1) >> 64)
    lo = w0 % UInt64 + (w1 << 64)

    hi, lo
end

end
