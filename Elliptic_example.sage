
"""
Elliptic_example.sage

Constructs AG-code on elliptic curve
Computes a triple (Y,Q,F) with MMP[1] attack
Computes an equivalent triple (Y',Q,F') using the shortened code

Global parameters:
F - the finite field
p = char(F)
R - a polynomial ring in two variables X,Y (used to construct the AG code)
v,y,R_S2c - lists of variables y[0]...y[k] for the MMP[1] routine
(!) Important: Variables y[0]...y[k] must be initialized as soon as we know dimensions of Generator matrix of the code
"""


F.<a> = GF(2^4)
p = F.characteristic()
# Variables to construct the AG code on Elliptic curve
R, (X, Y) = PolynomialRing(F, 2, 'XY', order='lex').objgens()
# Variables y[0]...y[k] for the MMP[1] routine
v=[]
R_S2c = PolynomialRing(F, len(v), v, order="degrevlex")
y=[]


load("include/MMP.sage")
load("include/Code_shortening.sage")


def Curve_construct(eq): 
    
    """
    Constructs elliptic curve defined over finite field F
    
    Input:
    eq - coefficients (EC) or the curve equation
    
    Output:
    curve - algebraic curve
    g - genus
    pts - rational points
    p_inf - point at infinity
    
    """
    
    global F,p,R,X,Y;
    
    curve = EllipticCurve(F, eq)
    pts = list(curve.points())
    p_inf = curve.point([0 , 1 , 0])
    pts.remove(p_inf)
    g = curve.genus()
    return curve, g, pts, p_inf
    
    
    

def RR_basis(m):
    
    """
    
    Constructs Rieman-Roch basis in case G = m*P_infty

    """
    
    global F,p,R,X,Y
    
    basis = [1,X,Y,X^2]
    flag = False
    for i in range(m//2 + 1):
        for j in range(2):
            if 2*i + 3*j <= m and X^i*Y^j not in basis:
                basis.append(X^i*Y^j)
                if len(basis) == m:
                    flag = True
                    break
        if flag == True:
            break
    return basis
    
   
    

def AG_Gen_matrix(crv, pts, E_div):
    
    """
    Constructs generator matrix of AG-code
    
    Input: 
    crv - algebraic curve
    pts - rational points on curve
    E_div - one-point divisor E=m*p_inf
    
    """
    
    global F,p,R,X,Y;
    
    assert len(E_div)==1, "E must be a one-point divisor!"
    m=E_div[0][0]
    rr_basis = RR_basis(m)
    n = len(pts)
    k = len(rr_basis)
    G = Matrix(F, k, n)
    for i in range(0, k):
        for j in range(n):
            if rr_basis[i] == 1:
                G[i, j] = 1
            else:
                G[i, j] = (rr_basis[i](pts[j][0], pts[j][1]))
    return G




# Select an elliptic curve over F
coeffs = [0, 0, a^7, a^5, a^2]
E, g, pts, p_inf = Curve_construct(coeffs)

# Construct AG-code on elliptic curve
m=7
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
print("+++++++++ (Y,Q,F) triple +++++++++")
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
#print(f"is_prime(I(Y'))={I_sh_hom.is_prime()}")

print(f"Q={Q}")
print("Divisor F:")
F_prime_div=Ideal(I_sh_gens_hom+[y[0]]) # Intersection with hyperplane that is disjoint of Q
print(F_prime_div)
print(f"dimF={F_prime_div.dimension()}")
print("++++++++++++++++++++++++++++++++++")
