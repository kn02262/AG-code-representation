
"""
MMP.sage

Contains routines for the attack described in 

[1] Márquez-Corbella, I., Martínez-Moro, E., Pellikaan, R., Ruano, D.: Computational
aspects of retrieving a representation of an algebraic geometry code. Journal of Symbolic
Computation 64, 67–87 (2014) https://doi.org/10.1016/j.jsc.2013.12.007 . Mathematical
and computer algebra techniques in cryptology

Global parameters:
F - the finite field;
p = char(F)
y[] - a list of variables y[0], ... , y[k]
R_S2c - a polynomial ring over F in variables y[0], ... , y[k]

"""

def init_vars(k):
	"""
	
	init_vars(k):
	Initializes global variables:
	
	y[] - a list of variables y[0], ... , y[k]
	R_S2c - a polynomial ring over F in variables y[0], ... , y[k]
	
	"""
	
	global F,p,v,R_S2c,y
	
	v = [SR.var('y{}'.format(i)) for i in range(k)]
	R_S2c = PolynomialRing(F, len(v), v, order="degrevlex")
	y = [R_S2c('y{}'.format(i)) for i in range(k)]
	
	return
	
	
	



def AG_test_VSAG(m,n,g):
	"""
	
	AG_test_VSAG(m,n,g):
	Checks whether the WAG code C_L(X,P,E) with parameters m=degE, n=[P], g=g(X) is VSAG or not.
	
	"""
	return (m in [2*g+2..floor_mod(n/2)])






def AG_test_VSAG_dual(m,n,g):
	"""
	
	AG_test_VSAG_dual(m,n,g):
	Checks whether the WAG code C_L(X,P,E)^dual with parameters m=degE, n=[P], g=g(X) is VSAG or not.
	
	"""
	return (m in [ceil_mod(n/2+2*g-2)..n-4])
	
	
	
	
	
	
def AG_Gen_equiv(G):
    """
    
    AG_Gen_equiv(G):
    For given generator matrix G (of size k*n),
    Outputs M*G*P, where M is a square k*k non-singular matrix, P is a permutation matrix
    
    """
    global F,p
    
    k=G.nrows()
    n=G.ncols()
    while True:
        M = random_matrix(F, k)
        if not M.is_singular():
            break
    P = Permutations(n).random_element()
    P = matrix(P)
    return M*G*P






def MMP_GetCrvGens(G):
    """
    
    For given generator matrix G (of dimension k*n),
    computes a generating set for ideal I(Y) with curve Y isomorphic to X.
    
    """
    
    global F,p,v,R_S2c,y
    
    k=G.nrows()
    n=G.ncols()
    assert len(y)>=k, "Global variables y[0]...y[k] are not initialized properly, forgot to call init_vars(k)?"
    
    # Compute G2 that represents schur product C*C
    G2m = binomial(k+1,2)
    G2n = n
    G2 = matrix(F, G2m, G2n)
    ln = 0
    for i in range(k):
        for j in range(i, k):
            for z in range(n):
                G2[ln, z] = G[i,z]*G[j,z]
            ln=ln+1
    
    # Compute left_kernel(G2):
    N2 = G2.left_kernel_matrix()
    N2n = N2.nrows()
    N2m = binomial(k+1,2)
    
    # Compute S2c_basis
    S2c_basis = []
    for i in range(k):
        for j in range(i,k):
            S2c_basis.append(y[i]*y[j])
    
    # Construct ideal I(Y):
    I2Gens = []
    for i in [0..N2n-1]:
        IGen = 0
        for j in [0..N2m-1]:
            IGen += N2[i][j]*S2c_basis[j]
        I2Gens.append(IGen)
    
    return I2Gens
