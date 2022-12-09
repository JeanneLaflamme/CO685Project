include("richelot_aux.jl")
include("uvtable.jl")
include("curves.jl")
include("helpers.jl")
include("attack.jl")
include("public_values_aux.jl")
using BenchmarkTools

SIKE_parameters = Dict([
    ("SIKEp434" , (216, 137)),
    ("SIKEp503" , (250, 159)),
    ("SIKEp610" , (305, 192)),
    ("SIKEp751" , (372, 239)),
    ("SIKEp964" , (486, 301)), # removed after NIST round 1
])

# Change me to attack different parameter sets
NIST_submission = "SIKEp434"
a, b = SIKE_parameters[NIST_submission]

println("Running the attack against $(NIST_submission) parameters, which has a prime: 2^$(a)*3^$(b) - 1")

println("Generating public data for the attack...")
# Set the prime, finite fields and starting curve
# with known endomorphism
p = ZZ(2)^a*ZZ(3)^b - 1
R,xx = PolynomialRing(GF(p),"x")
Fp2, ii = FiniteField(xx^2+1, "i")

E_start = EllipticCurve([0,6,0,1,0],Fp2)

# Generation of the endomorphism 2i
two_i(P) = generate_distortion_map(NIST_submission, E_start,P)

# Generate public torsion points, for SIKE implementations
# these are fixed but to save loading in constants we can
# just generate them on the fly
P2, Q2, P3, Q3 = generate_torsion_points(NIST_submission,E_start)

# Generate Bob's key pair
bob_private_key, EB, PB, QB = gen_bob_keypair(E_start, b, P2, Q2, P3, Q3)
solution = digits(Integer(bob_private_key),base=3)

println("If all goes well then the following digits should be found: $(solution)")

# ===================================
# =====  ATTACK  ====================
# ===================================
@time begin
	CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i)
end


