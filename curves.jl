using Nemo
using Random
import Base.:+
import Base.:-
import Base.*


#########################
# ELLIPTIC CURVES
#########################
struct EllipticCurve
    ainvs::Vector
    field
    
    function EllipticCurve(ainvs::Vector,field) 
        if length(ainvs)== 2
            ainvs = [0,0,0,ainvs[1],ainvs[2]]
        elseif length(ainvs) != 5
            error("Invalid number of coefficients: two or five should be given")
        end
        
        ainvs = [field(a) for a in ainvs]
        
        # a invariants
        a1,a2,a3,a4,a6 = ainvs

        # b invariants
        b2 = a1^2+4*a2
        b4 = 2*a4+a1*a3
        b6 = a3^2+4*a6
        b8 = a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2
        
        # c invariants
        c4 = b2^2-24*b4
        c6 = -b2^3+36*b2*b4-216*b6
        
        # Discriminant
        Δ = -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6
        
        if Δ == 0
            error("The curve you are trying to define is singular")
        end
        new(ainvs,field)
    end
end

# a invariants
function a1(E::EllipticCurve)
    return E.ainvs[1]
end

function a2(E::EllipticCurve)
    return E.ainvs[2]
end

function a3(E::EllipticCurve)
    return E.ainvs[3]
end

function a4(E::EllipticCurve)
    return E.ainvs[4]
end

function a6(E::EllipticCurve)
    return E.ainvs[5]
end

# b invariants
function b2(E::EllipticCurve)
    return a1(E)^2+4*a2(E)
end

function b4(E::EllipticCurve)
    return 2*a4(E)+a1(E)*a3(E)
end

function b6(E::EllipticCurve)
    return a3(E)^2+4*a6(E)
end

function b8(E::EllipticCurve)
    return a1(E)^2*a6(E)+4*a2(E)*a6(E)-a1(E)*a3(E)*a4(E)+a2(E)*a3(E)^2-a4(E)^2
end

# c invariants
function c4(E::EllipticCurve)
    return b2(E)^2-24*b4(E)
end

function c6(E::EllipticCurve)
    return -b2(E)^3+36*b2(E)*b4(E)-216*b6(E)
end

# Discriminant
function discriminant(E::EllipticCurve)
    return -b2(E)^2*b8(E)-8*b4(E)^3-27*b6(E)^2+9*b2(E)*b4(E)*b6(E)
end

# j-invariant
function jinvariant(E::EllipticCurve)
    return c4(E)^3//discriminant(E)
end

# Curve equation
f(x,y,E::EllipticCurve) = y^2+a1(E)*x*y+a3(E)*y-x^3-a2(E)*x^2-a4(E)*x-a6(E)

# Get base field
function field(E::EllipticCurve)
    return E.field
end

# Define equality of curves
function Base.:(==)(E1::EllipticCurve, E2::EllipticCurve)
    return (E1.ainvs == E2.ainvs) && (E1.field == E2.field)
end

# Pretty printing
function Base.show(io::IO, E::EllipticCurve)
    a = E.ainvs
    K = E.field
    
    s = "Elliptic Curve defined by y^2"
    
    if a[1] == K(-1)
        s *= " - x*y"
    elseif a[1] == K(1)
        s *= " + x*y"
    elseif a[1] != E.field(0)
        s *= " + $(a[1])*x*y"
    end
    if a[3] == K(-1)
        s *= " - y"
    elseif a[3] == K(1)
        s *= " + y"
    elseif a[3] != E.field(0)
        s *= " + $(a[3])*y"
    end
    s *= " = x^3"
    if a[2] == K(-1)
        s *= " - x^2"
    elseif a[2] == K(1)
        s *= " + x^2"
    elseif a[2] != E.field(0)
        s *= " + $(a[2])*x^2"
    end
    if a[4] == K(-1)
        s *= " - x"
    elseif a[4] == K(1)
        s *= " + x"
    elseif a[4] != E.field(0)
        s *= " + $(a[4])*x"
    end
    if a[5] == K(-1)
        s *= " - 1"
    elseif a[5] == K(1)
        s *= " + 1"
    elseif a[5] != E.field(0)
        s *= " + $(a[5])"
    end
    
    s*= " over $(K)"
    print(io,replace(s,"+ -" => "- "))
end

#############################
# POINTS ON ELLIPTIC CURVES
#############################

struct EllipticCurvePoint
    coords::Tuple
    curve::EllipticCurve
    
    function EllipticCurvePoint(coords::Tuple,curve::EllipticCurve)
        if length(coords)==2
            coords = (coords[1],coords[2],1)
        elseif length(coords) != 3
            error("Two or three coordinates should be given")
        end
        
        coords = (curve.field(coords[1]),curve.field(coords[2]),curve.field(coords[3]))
        valid = true
        
        if coords[3] != curve.field(0)
            x = coords[1]//coords[3]
            y = coords[2]//coords[3]
            z = curve.field(1)
            
            if f(x,y,curve) != curve.field(0)
                valid = false
            end
            
        elseif coords[1] == curve.field(0) && coords[2] != curve.field(0)
            x = curve.field(0)
            y = curve.field(1)
            z = curve.field(0)
        else
            valid = false
        end
        
        if !valid
            error("The point does not belong to the elliptic curve")
        end
        
        coords = (x,y,z)
        new(coords,curve)       
    end
    
end

# Extract coordinates
x(P::EllipticCurvePoint) = P.coords[1]
y(P::EllipticCurvePoint) = P.coords[2]
z(P::EllipticCurvePoint) = P.coords[3]

# Create a point using E(x,y)
function (E::EllipticCurve)(x,y)
    return EllipticCurvePoint((x,y),E)
end

function (E::EllipticCurve)(x,y,z)
    return EllipticCurvePoint((x,y,z),E)
end

function (E::EllipticCurve)(x)
    if x == 0
        return EllipticCurvePoint((0,1,0),E)
    else
        error("Please specify two or three coordinates")
    end
end

# Define equality of curve points
function Base.:(==)(P1::EllipticCurvePoint, P2::EllipticCurvePoint)
    return (P1.curve == P2.curve) && (P1.coords == P2.coords)
end

# Pretty printing
Base.show(io::IO, P::EllipticCurvePoint) = print(io,"($(x(P)) : $(y(P)) : $(z(P)))")

# Test if a point is zero
function iszero(P::EllipticCurvePoint)
    K = P.curve.field
    if P == P.curve(0)
        return true
    else
        return false
    end
end

# Negate a point
function -(P::EllipticCurvePoint)
    if iszero(P)
        return P
    else
        return EllipticCurvePoint((x(P),-y(P) - a1(P.curve)*x(P) - a3(P.curve)),P.curve)
    end
end

# Point addition
function +(P1::EllipticCurvePoint,P2::EllipticCurvePoint)
    if P1.curve != P2.curve
        error("The points do not belong to the same curve")
    end
    
    if iszero(P1)
        return P2
    elseif iszero(P2)
        return P1
    end
    
    E = P1.curve
    a1, a2, a3, a4, a6 = E.ainvs
    x1, y1 = x(P1), y(P1)
    x2, y2 = x(P2), y(P2)
    
    if x1 == x2 && y1 == -y2 - a1*x2 - a3
        return E(0)  
    elseif x1 == x2 && y1 == y2
        m = (3*x1*x1 + 2*a2*x1 + a4 - a1*y1) // (2*y1 + a1*x1 + a3)
    else
        m = (y1-y2)//(x1-x2)
    end
    if a1!=0 println(a1) end
    x3 = -x1 - x2 - a2 + m*(m+a1)
    y3 = -y1 - a3 - a1*x3 + m*(x1-x3)
    
    return EllipticCurvePoint((x3,y3), E)
end

# Point substraction
function -(P1::EllipticCurvePoint,P2::EllipticCurvePoint)
    return P1 + (-P2)
end

# Scalar multiplication
function *(n::Integer, P::EllipticCurvePoint)
    E = P.curve
    Q = P
    R = E(0)
    while n>0
        if n&1 == 1
            R = R+Q
        end
        Q = Q+Q
        n = n>>1
    end
    return R
end

# Generate a random point (for curves over finite field only)
function random_point(E::EllipticCurve)
    F = E.field
    if(characteristic(F)==0 || a1(E) != 0 || a3(E) != 0)
        error("random_point generation is only supported for elliptic curves over finite field of the form y^2 = ...")
    end
    while true
        a = rand(F)
        y2 = a^3+a2(E)*a^2+a4(E)*a+a6(E)
        b = 2*bitrand(1)[1]-1
        if issquare(y2)
            return E(a,b*sqrt(y2))
        end
    end
end

####################################
# ISOGENIES
####################################

struct EllipticCurveIsogeny
    #Elliptic curve isogeny of degree n with cyclic kernel generated by Q
    kergen::EllipticCurvePoint
    degree::Integer
    domain::EllipticCurve
    codomain::EllipticCurve
    
    function EllipticCurveIsogeny(Q::EllipticCurvePoint, n::Integer)
        kergen = Q
        degree = n
        domain = Q.curve
        
        a1,a2,a3,a4,a6 = domain.ainvs
        v,w = vwSum(Q,n)
        
        A4 = a4 - 5*v
        A6 = a6 - (a1^2+4*a2)*v - 7*w
        
        codomain = EllipticCurve([a1,a2,a3,A4,A6],domain.field)
        new(kergen,degree,domain,codomain)
    end
end

# These functions compute u and v in Velu's formulas
function uvQ(Q::EllipticCurvePoint)
    E = Q.curve
    
    fQ = 3*x(Q)^2+2*a2(E)*x(Q)+a4(E)-a1(E)*y(Q)
    gQ = -2*y(Q)-a1(E)*x(Q)-a3(E)
    uQ = gQ^2
    
    if iszero(2*Q)
        vQ = fQ
    else
        vQ = 2*fQ - a1(E)*gQ
    end
    
    return (uQ,vQ)
end

# Intermediate step for Velu's formulas
function vwSum(Q::EllipticCurvePoint, n::Integer)
    K = Q.curve.field
    v = K(0)
    w = K(0)
    for i=1:(n>>1)
        uQ, vQ = uvQ(Q)
        v += vQ
        w += uQ+x(Q)*vQ
        Q += Q
    end
    return (v,w)
end

# Evaluate isogeny phi on a point P
function (phi::EllipticCurveIsogeny)(P::EllipticCurvePoint)
    E = phi.codomain
    
    if iszero(P)
        return E(0)
    end
    
    Q = phi.kergen
    n = phi.degree
    a1,a2,a3,a4,a6 = phi.domain.ainvs
    X,Y = x(P),y(P)
    
    for i=1:(n>>1)
        if x(P)==x(Q) # Point is in Kernel
            return E(0) 
        end
        uQ, vQ = uvQ(Q)
        fQ = 3*x(Q)^2+2*a2*x(Q)+a4-a1*y(Q)
        gQ = -2*y(Q)-a1*x(Q)-a3
        X += vQ//(x(P)-x(Q)) + uQ//(x(P)-x(Q))^2
        Y -= uQ*(2*y(P)+a1*x(P)+a3)//(x(P)-x(Q))^3 + vQ*(a1*(x(P)-x(Q))+y(P)-y(Q))//(x(P)-x(Q))^2 + (a1*uQ-fQ*gQ)//(x(P)-x(Q))^2
        Q += Q
    end
    return E(X,Y)
end

###################################
#  PAIRINGS
###################################

# Computes the value at P of line through P and R (translated from sage source code)
function line(P::EllipticCurvePoint, R::EllipticCurvePoint, Q::EllipticCurvePoint)
    if iszero(Q)
        error("Q must be nonzero.")
    end

    if iszero(P) || iszero(R)
        if P == R
            return P.curve.field(1)
        elseif iszero(P)
            return x(Q) - x(R)
        elseif iszero(R)
            return x(Q) - x(P)
        end
    elseif P != R
        if x(P) == x(R)
            return x(Q) - x(P)
        else
            l = (y(R) - y(P))//(x(R) - x(P))
            return y(Q) - y(P) - l * (x(Q) - x(P))
        end
    else
        a1, a2, a3, a4, a6 = P.curve.ainvs
        numerator = (3*x(P)^2 + 2*a2*x(P) + a4 - a1*y(P))
        denominator = (2*y(P) + a1*y(P) + a3)
        if denominator == 0
            return x(Q) - x(P)
        else
            l = numerator//denominator
            return y(Q) - y(P) - l * (x(Q) - x(P))
        end
    end
end

# Miller's algorithm to compute f_P(A_Q) (translated from sage source code
function miller(P::EllipticCurvePoint, Q::EllipticCurvePoint, l::Integer)
    if P.curve != Q.curve
        error("Points must be on the same curve")
    end
    E = P.curve
    F = E.field
    q = characteristic(F)
    if q == 0
        error("Pairing only implemented for elliptic curves over finite fields")
    end
    k = degree(F)
    V = P
    f1 = 1
    f2 = 1
    t = 1
    lbin = digits(l, base=2)
    i = length(lbin)-1
    while i > 0
        S = 2*V
        ell = line(V,V,Q)
        vee = line(S,-S,Q)
        t = t^2*(ell//vee)
        V = S
        if lbin[i] == 1
            S = V+P
            ell = line(V,P, Q)
            vee = line(S,-S, Q)
            t = t*(ell//vee)
            V = S
        end
        i = i-1
    end
    return t
end

# Computes the weil pairing of P and Q (translated from sage source code)
function weil_pairing(P::EllipticCurvePoint, Q::EllipticCurvePoint, n::Integer)
    E = P.curve

    if Q.curve != E
        error("points must both be on the same curve")
    end


    # Test if P, Q are both in E[n]
    if n*P != E(0) || n*Q != E(0)
        error("points must both be n-torsion")
    end

    one = E.field(1)

    # Case where P = Q
    if P == Q
        return one
    # Case where P = O or Q = O
    elseif iszero(P) || iszero(Q)
        return one
    end

    num = miller(P,Q,n)
    den = miller(Q,P,n)
    if den == 0
        return one
    else
        return (-1)^(n&1)*num // den
    end
end
    


####################################
# DIVISORS ON HYPERELLIPTIC CURVES
####################################

# Hyperelliptic curve y^2 = f(x)
struct HyperellipticCurve
    f::PolyElem
end

# Pretty printing of a hyperelliptic curve
function Base.show(io::IO, H::HyperellipticCurve)
    print("Hyperelliptic Curve defined by y^2 = $(H.f) over $(base_ring(H.f))")
end

# Divisor on a hyperelliptic curve stored in mumford coordinates [u,v]
struct Divisor
    u::PolyElem
    v::PolyElem
    H::HyperellipticCurve
    function Divisor(u::PolyElem,v::PolyElem,H::HyperellipticCurve)
        u = div(u,coefficients(u)[degree(u)]) # Make u monic
        v = rem(v,u) # Reduce v to a degree one polynomial
        new(u,v,H)
    end
end

# Addition of two divisors using Cantor composition. Follows https://en.wikipedia.org/wiki/Imaginary_hyperelliptic_curve#The_divisor_and_the_Jacobian
function +(D1::Divisor,D2::Divisor)
    u1, v1 = D1.u, D1.v
    u2, v2 = D2.u, D2.v
    H = D1.H
    f = H.f
    g = 2
    x = gen(parent(f))
    
    d1, e1, e2 = gcdx(u1,u2)
    d, c1, c2 = gcdx(d1,v1+v2)
    s1, s2, s3 = c1*e1, c1*e2, c2
    
    u = div(u1*u2,d^2)
    v = rem(div(s1*u1*v2+s2*u2*v1+s3*(v1*v2+f),d),u)
    up = x^(g+1)
    vp = 0
    while degree(up) > g
        up = div(f-v^2,u)
        vp = rem(-v,up)
    end
    up = div(up, coefficients(up)[degree(up)])
    
    return Divisor(up,vp,H)
end

# Pretty printing of a divisors    
function Base.show(io::IO, D::Divisor)
    print("[$(D.u), $(D.v)]")
end
    
    