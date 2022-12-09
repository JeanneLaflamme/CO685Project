using Nemo
using Mods
include("curves.jl")

function FromProdToJac(C::EllipticCurve, E::EllipticCurve, P_c::EllipticCurvePoint, Q_c::EllipticCurvePoint, P::EllipticCurvePoint, Q::EllipticCurvePoint, a::Integer)
    Fp2 = field(C)
    Rx,X = PolynomialRing(Fp2,"x")
    MS1 = MatrixSpace(Fp2, 3, 3)
    MS2 = MatrixSpace(Fp2, 3, 1)

    P_c2 = BigInt(2)^(a-1)*P_c
    Q_c2 = BigInt(2)^(a-1)*Q_c
    P2 = BigInt(2)^(a-1)*P
    Q2 = BigInt(2)^(a-1)*Q

    a1, a2, a3 = x(P_c2), x(Q_c2), x(P_c2 + Q_c2)
    b1, b2, b3 = x(P2), x(Q2), x(P2 + Q2)

    # Compute coefficients
    M = MS1([
        a1*b1 a1 b1 ;
        a2*b2 a2 b2 ;
        a3*b3 a3 b3])
    b = MS2([1;1;1])
    R, S, T = solve(M,b)
    RD = R * det(M)
    da = (a1 - a2)*(a2 - a3)*(a3 - a1)
    db = (b1 - b2)*(b2 - b3)*(b3 - b1)

    s1, t1 = - da // RD, db // RD
    s2, t2 = -T//R, -S//R

    a1_t = (a1 - s2) // s1
    a2_t = (a2 - s2) // s1
    a3_t = (a3 - s2) // s1
    h = s1 * (X^2 - a1_t) * (X^2 - a2_t) * (X^2 - a3_t)

    H = HyperellipticCurve(h)

    function isogeny(pair::Tuple)
        # Argument members may be None to indicate the zero point.

        # The projection maps are:
        # H->C: (xC = s1/x²+s2, yC = s1 y)
        # so we compute Mumford coordinates of the divisor f^-1(P_c): a(x), y-b(x)
        Pc, P = pair
        if Pc != nothing
            xPc, yPc = x(Pc), y(Pc)
            JPc = Divisor(s1 * X^2 + s2 - xPc, Rx(yPc //s1), H)
        end
        # Same for E
        # H->E: (xE = t1 x² + t2, yE = t1 y/x^3)
        if P != nothing
            xP, yP = x(P), y(P)
            JP = Divisor((xP - t2) * X^2 - t1, yP *div(X^3,t1), H)
        end
        if Pc != nothing && P != nothing
            return JPc + JP
        elseif Pc != nothing
            return JPc
        elseif P != nothing
            return JP
        end
    end

    imPcP = isogeny((P_c, P))
    imQcQ = isogeny((Q_c, Q))

    return H, imPcP,  imQcQ, isogeny
end

function RichelotCorr(G1, G2, H1, H2, hnew, D)
        "Computes the image of D"
        U, V = D.u, D.v
        U = div(U,coefficients(U)[2])
        V = rem(V,U)
        H11 = H1*H1
        H12 = H1*H2
        H22 = H2*H2
        Hnew = HyperellipticCurve(hnew)
        x = gen(parent(hnew))
        # Sum and product of (xa, xb)
        s, p = -coefficients(U)[1], coefficients(U)[0]
        # Compute X coordinates (non reduced, degree 4)
        g1red = G1 - U
        g2red = G2 - U
        g11, g10 = coefficients(g1red)[1], coefficients(g1red)[0]
        g21, g20 = coefficients(g2red)[1], coefficients(g2red)[0]
        # see above
        Px = (g11*g11*p + g11*g10*s + g10*g10) * H11+ 
           (2*g11*g21*p + (g11*g20+g21*g10)*s + 2*g10*g20) * H12+
           (g21*g21*p + g21*g20*s + g20*g20) * H22

        # Compute Y coordinates (non reduced, degree 3)
        v1, v0 = coefficients(V)[1], coefficients(V)[0]
        Py2 = v1*v1*p + v1*v0*s + v0*v0
        Py1 = (2*v1*g11*p + v1*g10*s + v0*g11*s + 2*v0*g10)*x - (v1*g11*s*p + 2*v1*g10*p + v0*g11*(s*s-2*p) + v0*g10*s)
        Py1 *= H1
        Py0 = H11 * U * (g11*g11*p + g11*g10*s + g10*g10)

        # Now reduce the divisor, and compute Cantor reduction.
        _, Py1inv, _ = gcdx(Py1,Px)
        Py = rem(- Py1inv * (Py2 * hnew + Py0), Px)

        Dx = div(hnew - Py^2, Px)
        Dy = rem(-Py, Dx)

        return Divisor(Dx, Dy, Hnew)
end

function FromJacToJac(H::HyperellipticCurve, D1::Divisor, D2::Divisor, a::Integer)
    h = H.f
    x = gen(parent(h))
    Fp2 = base_ring(h)
    MS = MatrixSpace(Fp2,3,3)

    G1, _ = jacobian_iter_double(h, D1.u, D1.v, a-1)
    G2, _ = jacobian_iter_double(h, D2.u, D2.v, a-1)
    G3, r3 = divrem(h, G1 * G2)

    delta = MS([coefficients(G)[i] for G in [G1,G2,G3], i = 0:2])
    delta = inv(delta)
    H1 = -delta[1,1]*x^2 + 2*delta[2,1]*x - delta[3,1]
    H2 = -delta[1,2]*x^2 + 2*delta[2,2]*x - delta[3,2]
    H3 = -delta[1,3]*x^2 + 2*delta[2,3]*x - delta[3,3]

    hnew = H1*H2*H3
    Hnew = HyperellipticCurve(hnew)
   
  
    # Now compute image points:
    isogeny(D) =  RichelotCorr(G1, G2, H1, H2, hnew, D)                        
    imD1 = isogeny(D1)
    imD2 = isogeny(D2)

    return Hnew, imD1, imD2, isogeny
end

function jacobian_double(h, u, v)
    """
    Computes the double of a jacobian point (u,v)
    given by Mumford coordinates: except that u is not required
    to be monic, to avoid redundant reduction during repeated doubling.
    See SAGE cantor_composition() and cantor_reduction
    """
    @assert degree(u) == 2
    # Replace u by u^2
    # Compute h3 the inverse of 2*v modulo u
    # Replace v by (v + h3 * (h - v^2)) % u
    R = parent(u)
    q, r = divrem(u,2*v)
    if coefficients(r)[0] == 0 # gcd(u, v) = v, very improbable
        a = q^2
        b = rem(div(v + (h - v^2), v), a)
        return a, b
    else # gcd(u, v) = 1
        h3 = 1 // (-coefficients(r)[0]) * q
        a = u*u
        b = rem(v + h3 * (h - v^2), a)
        # Cantor reduction
        Dx = div(h - b^2, a)
        Dy = rem(-b, Dx)

        return Dx, Dy
    end
end

function jacobian_iter_double(h, u, v, n)
    for i=1:n
        u, v = jacobian_double(h, u, v)
    end

    return div(u,coefficients(u)[degree(u)]), v
end

function FromJacToProd(G1::PolyElem, G2::PolyElem, G3::PolyElem)
    """
    Construct the "split" isogeny from Jac(y^2 = G1*G2*G3)
    to a product of elliptic curves.
    This computation is the same as Benjamin Smith
    see 8.3 in http://iml.univ-mrs.fr/~kohel/phd/thesis_smith.pdf
    """
    h = G1*G2*G3
    R = parent(h)
    Fp2 = base_ring(R)
    xx = gen(R)
    MS = MatrixSpace(Fp2,3,3)
    VS = MatrixSpace(Fp2,3,1)

    M = MS([coefficients(G)[i] for G in [G1,G2,G3], i = 0:2])
    # Find homography
    u,v,w = right_kernel(M)[2]
    v = v//u
    w = w//u
    u = u//u
    
    d = u//2
    ad,b = roots(xx^2 - v*xx + w*d//2)
    a = ad//d

    # Apply transform G(x) -> G((a*x+b)/(x+d))*(x+d)^2
    # The coefficients of x^2 are M * (1, a, a^2)
    # The coefficients of 1 are M * (d^2, b*d, b^2)
    H11, H21, H31 = M*VS([Fp2(1);a;a*a])
    H10, H20, H30 = M * VS([d*d; b*d; b*b])
    @assert G1((a*xx+b)//(xx+d))*(xx+d)^2 == H11*xx^2+H10

    h2 = (H11*xx^2+H10)*(H21*xx^2+H20)*(H31*xx^2+H30)
    H2 = HyperellipticCurve(h2)

    p1 = (H11*xx+H10)*(H21*xx+H20)*(H31*xx+H30)
    p2 = (H11+H10*xx)*(H21+H20*xx)*(H31+H30*xx)
    # We will need to map to actual elliptic curve
    p1norm = (xx + H10*H21*H31)*(xx + H20*H11*H31)*(xx + H30*H11*H21)
    p2norm = (xx + H11*H20*H30)*(xx + H21*H10*H30)*(xx + H31*H10*H20)
    p1norm = coefficients(p1norm)
    p2norm = coefficients(p2norm)
    E1 = EllipticCurve([0, p1norm[2], 0, p1norm[1], p1norm[0]],Fp2)
    E2 = EllipticCurve([0, p2norm[2], 0, p2norm[1], p2norm[0]],Fp2)

    morphE1(x, y)=(H11*H21*H31*x, H11*H21*H31*y)
    morphE2(x, y)=(H10*H20*H30*x, H10*H20*H30*y)
    # The morphisms are:
    # inverse homography:
    # H->H2: x, y => ((b-dx) / (x-a), y/(x-a)^3)
    # then H2->E1:(x,y) => (x^2,y)
    #   or H2->E2:(x,y) => (1/x^2,y/x^3)

    function isogeny(D::Divisor)
        # To map a divisor, perform the change of coordinates
        # on Mumford coordinates
        U, V = coefficients(D.u), coefficients(D.v)
        # apply homography
        # y = v1 x + v0 => 
        U_ = U[0] * (xx+d)^2 + U[1]*(a*xx+b)*(xx+d) + U[2]*(a*xx+b)^2
        V_ = V[0] * (xx+d)^3 + V[1]*(a*xx+b)*(xx+d)^2
        V_ = rem(V_,U_)
        U_ = coefficients(U_)
        V_ = coefficients(V_)
        v1, v0 = V_[1], V_[0]
        # prepare symmetric functions
        s = - U_[1] // U_[2]
        p = U_[0] // U_[2]
        # compute Mumford coordinates on E1
        # Points x1, x2 map to x1^2, x2^2
        U1 = xx^2 - (s*s - 2*p)*xx + p^2
        # y = v1 x + v0 becomes (y - v0)^2 = v1^2 x^2
        # so 2v0 y-v0^2 = p1 - v1^2 xH^2 = p1 - v1^2 xE1
        V1 = div(p1 - v1^2 * xx + v0^2,2*v0)
        # Reduce Mumford coordinates to get a E1 point
        V1 = rem(V1, U1)
        U1red = div(p1 - V1^2,U1)
        U1red = coefficients(U1red)
        xP1 = -U1red[0] // U1red[1]
        yP1 = V1(xP1)
        @assert yP1^2 == p1(xP1)
        # Same for E2
        # Points x1, x2 map to 1/x1^2, 1/x2^2
        U2 = xx^2 - (s*s-2*p)//p^2*xx + 1//p^2
        # yE = y1/x1^3, xE = 1/x1^2
        # means yE = y1 x1 xE^2
        # (yE - y1 x1 xE^2)(yE - y2 x2 xE^2) = 0
        # p2 - yE (x1 y1 + x2 y2) xE^2 + (x1 y1 x2 y2 xE^4) = 0
        V21 = xx^2 * (v1 * (s*s-2*p) + v0*s)
        V20 = p2 + xx^4 * (p*(v1^2*p + v1*v0*s + v0^2))
        # V21 * y = V20
        _, V21inv, _ = gcdx(V21,U2)
        V2 = rem(V21inv * V20, U2)
        #assert V2**2 % U2 == p2 % U2
        # Reduce coordinates
        U2red = div(p2 - V2^2,U2)
        U2red = coefficients(U2red)
        xP2 = -U2red[0] // U2red[1]
        yP2 = V2(xP2)
        
        xP1, yP1 = morphE1(xP1, yP1)
        xP2, yP2 = morphE2(xP2, yP2)

        return (E1(xP1,yP1), E2(xP2, yP2))
    end

    return isogeny, (E1, E2)
end

function Does22ChainSplit(C::EllipticCurve, E::EllipticCurve, P_c::EllipticCurvePoint, Q_c::EllipticCurvePoint,
        P::EllipticCurvePoint, Q::EllipticCurvePoint, a::Integer)

    chain = []
    # gluing step
    H, D1, D2, f = FromProdToJac(C, E, P_c, Q_c, P, Q, a)
    push!(chain,f)
    
    for i in range(1,a-2)
        H, D1, D2, f = FromJacToJac(H, D1, D2, a-i)
        push!(chain,f)
    end
        
    # now we are left with a quadratic splitting: is it singular?
    Fp2 = base_ring(H.f)
    MS = MatrixSpace(Fp2,3,3)
    
    G1 = D1.u
    G2 = D2.u
    G3, r3 = divrem(H.f, G1 * G2)
    @assert r3 == 0

    delta = MS([coefficients(G)[i] for G in [G1,G2,G3], i = 0:2])
    if det(delta) !=0
        return nothing
    end
    
    # Finish chain
    f, codomain = FromJacToProd(G1, G2, G3)
    push!(chain,f)
    return chain, codomain
end

function Pushing3Chain(E, P, i)
    """
    Compute chain of isogenies quotienting
    out a point P of order 3^i
    https://trac.sagemath.org/ticket/34239
    """
    function rec(Q, k)
        if k == 1  # base case
            return [EllipticCurveIsogeny(Q, 3)]
        end

        k1 = Int(floor(k * .8 + .5))
        k1 = max(1, min(k-1, k1))  # clamp to [1, k-1]

        Q1 = BigInt(3)^k1 * Q
        L = rec(Q1, k-k1)
    
        Q2 = Q
        for psi in L
            Q2 = psi(Q2)
        end
        R = rec(Q2, k1)

        return vcat(L,R)
    end

    chain = rec(P, i)
    return last(chain).codomain, chain
end

function AuxiliaryIsogeny(i::Int, u::BigInt, v::BigInt, E_start::EllipticCurve, P2::EllipticCurvePoint, Q2::EllipticCurvePoint, tauhatkernel::EllipticCurvePoint, two_i)
    """
    Compute the distored  kernel using precomputed u,v and the
    automorphism two_i.
    This is used to construct the curve C from E_start and we
    compute the image of the points P_c and Q_c
    """
    tauhatkernel_distort = u*tauhatkernel + v*two_i(tauhatkernel)
    C, tau_tilde = Pushing3Chain(E_start, tauhatkernel_distort, i)
    function chain(P)
        global Pc = u*P + v*two_i(P)
        for taut in tau_tilde
            Pc = taut(Pc)
        end
        return Pc
    end
    return C, chain(P2), chain(Q2), chain
end