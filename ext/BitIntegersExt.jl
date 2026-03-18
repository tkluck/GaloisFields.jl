module BitIntegersExt

import BitIntegers
import GaloisFields.widen_bits

widen_bits(x::Int128) = BitIntegers.Int256(x)
widen_bits(::Type{Int128}) = BitIntegers.Int256

end #module
