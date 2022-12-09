include("richelot_aux.jl")
include("uvtable.jl")
include("curves.jl")

function CastryckDecruAttack(E_start::EllipticCurve, P2::EllipticCurvePoint, Q2::EllipticCurvePoint, EB::EllipticCurve, PB::EllipticCurvePoint, QB::EllipticCurvePoint, two_i)
    skB = [] # TERNARY DIGITS IN EXPANSION OF BOB'S SECRET KEY

    # gathering the alpha_i, u, v from table
	expdata = [[BigInt(0), BigInt(0), BigInt(0)] for i=1:b-3]
	for i=1:b-3
    	if (b-i)%2 == 1
        	index = Int((b - i + 1) // 2)
        	expo = uvtable[index][2]
        	if expo <= a
            	global u = uvtable[index][3]
            	global v = uvtable[index][4]
            	expdata[i] = [expo, u, v]
        	end
    	end
	end

	# gather digits until beta_1
	global bet1 = 1
	while expdata[bet1][1] == 0
    	global bet1 += 1
	end

	ai = expdata[bet1][1]
	u  = expdata[bet1][2]
	v  = expdata[bet1][3]

	println("Determination of first $(bet1) ternary digits. We are working with 2^$(ai)-torsion.")

	bi = b - bet1
	alp = a - ai

    function CheckGuess(first_digits)
        println("Testing digits: $(first_digits)")

        tauhatkernel = BigInt(3)^bi*P3
    	for k=1:bet1
       		tauhatkernel += (BigInt(3)^(k-1)*first_digits[k])*BigInt(3)^bi*Q3
    	end
        tauhatkernel_distort = (u*tauhatkernel) + (v*two_i(tauhatkernel))

        C, P_c, Q_c, chainC = AuxiliaryIsogeny(bet1, u, v, E_start, P2, Q2, tauhatkernel, two_i)
        
        # We have a diagram
        #  C <- Eguess <- E_start
        #  |    |
        #  v    v
        #  CB-> EB
        split = Does22ChainSplit(C, EB, 2^alp*P_c, 2^alp*Q_c, 2^alp*PB, 2^alp*QB, ai)
        if split != nothing
            Eguess, _ = Pushing3Chain(E_start, tauhatkernel, bet1)

            chain, (E1, E2) = split
            # Compute the 3^b torsion in C
            P3c = chainC(P3)
            Q3c = chainC(Q3)
            
          
            # Map it through the (2,2)-isogeny chain
            if jinvariant(E2) == jinvariant(Eguess)
                CB, index = E1, 1
            else
                CB, index = E2, 2
            end
            function apply_chain(c, X::EllipticCurvePoint)
                X = (X, nothing) # map point to C x {O_EB}
                for f in c
                    X = f(X)
                end
                return X[index]
            end
            println("Computing image of 3-adic torsion in split factor CB")
            P3c_CB = apply_chain(chain, P3c)
            Q3c_CB = apply_chain(chain, Q3c)
            

            m = BigInt(3)^b
            # Determine kernel of the 3^b isogeny.
            # The projection to CB must have 3-adic rank 1.
            # To compute the kernel we choose a symplectic basis of the
            # 3-torsion at the destination, and compute Weil pairings.
            P_CB, Q_CB = supersingular_gens(CB)
            K = div(BigInt(p+1), BigInt(3)^b)
            P3_CB =  K*P_CB
            Q3_CB = K*Q_CB
            
            w = weil_pairing(P3_CB,Q3_CB, BigInt(3)^b)
            # Compute kernel
            for G in (P3_CB, Q3_CB)
                xP = fast_log3(weil_pairing(P3c_CB,G, BigInt(3)^b), w)
                xQ = fast_log3(weil_pairing(Q3c_CB,G, BigInt(3)^b), w)
                if xQ % 3 != 0
                    sk = mod(-xP*invmod(xQ,m), m)
                    return sk
                end
            end
            return true
        end
    end

    guesses = product([0,1,2], bet1)

    for guess in guesses
    	result = CheckGuess(guess)
        sk = result
        if sk != nothing
            println("Glue-and-split! These are most likely the secret digits.")
            global bobskey = sk
            break
        end
    end

    # Sanity check
    bobscurve, _ = Pushing3Chain(E_start, P3 + bobskey*Q3, b)
    found = jinvariant(bobscurve) == jinvariant(EB)

    if found
        println("Bob's secret key revealed as: $(bobskey)")
        println("In ternary, this is: $(digits(Integer(bobskey),base=3))")
        return bobskey
    else
        println("Something went wrong.")
        return nothing
    end
end