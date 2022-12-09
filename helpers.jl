include("curves.jl")
using Mods

########################
#OTHER HELPERS
#######################
function product(V::Vector, n::Integer)
    W = V
    for i=1:n-1
        W = [vcat(i,j) for i in V for j in W]
    end
    return W
end

function supersingular_gens(E::EllipticCurve)
    """
    Compute generators of E, assuming E is supersingular
    with smooth order (p+1)^2 with factors 2 and 3 only.
    This is faster than the PARI method.
    """
    # Find a random point of order (p+1) (probability 1/3)
    p = BigInt(characteristic(E.field))
    k2 = div(p+1,2)
    k3 = div(p+1,3)
    while true
        global P = random_point(E)
        if k2*P != E(0) && k3*P != E(0)
            break
        end
    end

    while true
        global Q = random_point(E)
        if k2*Q != 0 && k3*Q != 0
            # but is it linearly independent? (probability 1/3)
            w = weil_pairing(P,Q, p+1)
            #println(w)
            if w^k2 != 1 && w^k3 != 1
                return P, Q
            end
        end
    end
end

function fast_log3(x, _base)
    """
    Fast discrete log when elements are known to have order
    dividing 3^k
    """
    _one = one(parent(x))
    powers = [_base]
    b = _base
    log_order = nothing
    for i = 0:10_000
        b = b^3
        if b == _one
            log_order = i+1
            break
        end
        push!(powers,b)
    end
    if b != _one
        error("impossible")
    end
    _digits = []
    @assert x^(BigInt(3)^log_order) == 1
    @assert _base^(BigInt(3)^(log_order-1)) != 1
    for i=1:log_order
        for d=0:2
            if (x * powers[i]^d)^(BigInt(3)^(log_order-i)) == 1
                push!(_digits, mod(-d, 3))
                if d != 0 
                    x = div(x,powers[i]^(3-d))
                end
                break
            end
        end
        if x == 1
            break
        end
    end
    dlog = sum(d[2]*BigInt(3)^(d[1]-1) for d in enumerate(_digits))
    return dlog
end