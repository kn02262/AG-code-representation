
"""
Hermitian_example.sage

Constructs AG-code on Hermitian curve
Computes a triple (Y,Q,F) with MMP[1] attack
Computes an equivalent triple (Y',Q,F') using the shortened code

Global parameters:
F - the finite field;
p = char(F)
v,y,R_S2c - lists of variables y[0]...y[k] for the MMP[1] routine
(!) Important: Variables y[0]...y[k] must be initialized as soon as we know dimensions of Generator matrix of the code
"""

F.<a> = GF(3^2)
p = F.characteristic()
# Variables y[0]...y[k] for the MMP[1] routine
v=[]
R_S2c = PolynomialRing(F, len(v), v, order="degrevlex")
y=[]


load("include/MMP.sage")
load("include/Code_shortening.sage")


def Hermitian_genus(q): # Genus of the Hermitian curve over F_{p^2}
	return q*(q-1)/2


def Curve_construct(): 
	
	"""
	Constructs Hermitian curve defined over finite field F
	
	Output:
	curve - algebraic curve
	g - genus
	pts - rational points
	p_inf - point at infinity
	
	"""
	
	global F,p
	
	P.<X,Y,Z> = ProjectiveSpace(F, 2)
	curve = Curve([Y^p*Z+Y*Z^p-X^(p+1)], P)
	g=Hermitian_genus(p)
	RatPts = curve.rational_points()
	for Q in RatPts: # Remove points at infinity
		if Q[2] == 0:
			RatPts.remove(Q)
	#print(curve)
	#print(g)
	#print(RatPts)
	#print(Q)

	return curve, g, RatPts, Q
	
	
	
	
def AG_Gen_matrix(crv, pts, E_div):
	
	"""
	Constructs generator matrix of AG-code
	
	Input: 
	crv - algebraic curve
	pts - rational points on curve
	E_div - one-point divisor E=m*p_inf
	
	"""
	
	global F,p
	
	assert len(E_div)==1, "E must be a one-point divisor!"
	n = len(pts)
	m=E_div[0][0]
	k = m+1-g
	
	# Construct a generator matrix
	G=matrix(F, k, n)
	# Precompute indeces (i,j) that form basis elements x^i*y^j of RR-Space
	ind=[]
	cnt=0
	for i in [0..9999999]:
		for j in [0..p-1]:
			if(i*p+j*(p+1) <= m):
				ind.append([i,j])
				cnt+=1
			if(cnt >= k):
				break
		if(cnt >= k):
			break
	assert len(ind)==k, "Riemann-Roch basis brecomputation failed!"
	#print(ind)	

	for i in [0..n-1]:
		# Values of coordinate functions
		x=pts[i][0] / pts[i][2]
		y=pts[i][1] / pts[i][2]
		#print(x,y)
		for j in [0..k-1]:
			G[j,i] = x^ind[j][0]*y^ind[j][1]
	return G
	
	
	

# Cunstruct the Hermitian curve over F
E, g, pts, p_inf = Curve_construct()

# Construct AG-code on Hermitian curve
m=2*g+3
G=AG_Gen_matrix(E, pts, E.divisor([ (m, p_inf) ]) )
#G=AG_Gen_equiv(G)
n=G.ncols()
k=G.nrows()
init_vars(k)
print(f"Constructed AG-code of length n={n} and dimension k={k}")
I_gens = MMP_GetCrvGens(G)
I=Ideal(I_gens)
# List places Q
Q=List_places(G)
print("+++++++++ (Y,Q,F) triple +++++++++");
print("Ideal I(Y) constructed following MMP:")
print(I)
print(f"dimI(Y)={I.dimension()}")
print(f"Q={Q}")
print("Divisor F:")
F_div=Ideal(I+[y[0]]) # Intersection with hyperplane that is disjoint of Q
print(F_div)
print(f"dimF={F_div.dimension()}")
print("++++++++++++++++++++++++++++++++++")

# Try to shorten the code
t = Heuristic_t(G)
assert t>=1, "Shortening heuristic failed, could not shorten"
G_sh=G
for i in [1..t]:
	G_sh=G_sh.delete_rows([G_sh.nrows()-1]).delete_columns([G_sh.ncols()-1])
#print(n,k, G.ncols(), G.nrows(), G_sh.ncols(), G_sh.nrows())
print(f"Shortening on t={t} positions")
print(f"Constructed shortened AG-code of length n={n-t} and dimension k={k-t}")
I_sh_gens = MMP_GetCrvGens(G_sh)
print("Generators of I(Y') constructed following MMP:")
print(Ideal(I_sh_gens))
I_sh_gens_affine = []
for f in I_sh_gens:
	I_sh_gens_affine.append(f.subs(y0=1))
print("Generators of I(Y') affine part with y0=1:")
print(Ideal(I_sh_gens_affine))
print(f"dimI(Y')={Ideal(I_sh_gens_affine).dimension()}")
print(f"dimI(Y'_affine)={Ideal(I_sh_gens_affine).dimension()}")
print("Polynomials lifting I(Y'):")
for i in [k-t..k-1]:
	f = lifting_polynomials(G, i)
	print(f)
	I_sh_gens_affine.append(f)
	print(f"dimI(Y')={Ideal(I_sh_gens_affine).dimension()}")
I_sh_affine=Ideal(I_sh_gens_affine)

# Homogenizing ideal I_sh
GB = I_sh_affine.groebner_basis()
I_sh_gens_hom = []
for f in GB:
	I_sh_gens_hom.append(f.homogenize(var=y[0]))
I_sh_hom=Ideal(I_sh_gens_hom)

print("+++++++++ (Y',Q,F') triple +++++++++")
print("Ideal I(Y') (affine) constructed following Prop.15:")
print(I_sh_affine)
print(f"dimI(Y')={I_sh_affine.dimension()}")
print(f"is_prime(I(Y'))={I_sh_affine.is_prime()}")

print("Ideal I(Y') (projective curve) constructed following Prop.15:")
print(I_sh_hom)
print(f"dimI(Y')={I_sh_hom.dimension()}")
print(f"is_prime(I(Y'))={I_sh_hom.is_prime()}")

print(f"Q={Q}")
print("Divisor F:")
F_prime_div=Ideal(I_sh_gens_hom+[y[0]]) # Intersection with hyperplane that is disjoint of Q
print(F_prime_div)
print(f"dimF={F_prime_div.dimension()}")
print("++++++++++++++++++++++++++++++++++")
