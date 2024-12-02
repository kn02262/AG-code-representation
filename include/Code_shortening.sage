
"""
Code_shortening.sage

Contains routines for shortening the code and heuristic to determine the maximum shortening value

Global parameters:
F - the finite field
p = char(F)
y[] - a list of variables y[0], ... , y[k]
R_S2c - a polynomial ring over F in variables y[0], ... , y[k]

"""

def floor_mod(a):
	if frac(a) == 0:
		return a-1
	else:
		return floor(a)

def ceil_mod(a):
	if frac(a) == 0:
		return a+1
	else:
		return ceil(a)






def List_places(G):
	
	"""
	List_places(G):
	
	Returns list of places Q (projective coordinates are columns of generator matrix G)
	
	"""
	
	global F,p
	QPlaces = []
	for s in [0..len(G.columns())-1]:
		QPlaces.append(list(G.column(s)))
	return QPlaces






def Heuristic_t(G):

	"""
	
	Heuristic_t(G,q)
	Heuristically computes maximal value t, such that G_t is a VSAG code
	Input:
		G is a Generator matrix G of size (k * n) defined over F_{q^2}
	Output: 
		t>=1 or 'false'

	"""
	global F,p
	n=G.ncols()
	k=G.nrows()
	t=k-1
	pre_check = false
	cur_check = false
	pre_g=0
	cur_g=0
	pre_m=0
	cur_m=0
	while t>=1:
		G_1=G.echelon_form()
		# Remove t columns and t rows from G_1:
		for i in [1..t]:
			s=len(G_1.rows())
			G_1=G_1.delete_rows([1])
			s=len(G_1.columns())
			G_1=G_1.delete_columns([1])
		G_2=[]
		rows = G_1.rows()
		for i in [0..len(rows)-1]:
			for j in [i..len(rows)-1]:
				G_2.append(rows[i].pairwise_product(rows[j]))
		#print(len(G_2), binomial(len(G_1.rows())+1,2))
		G_2=Matrix(G_2)
		# Parameters
		k1_t=G_1.rank()
		k2_t=G_2.rank()
		mt = k + k2_t - 2*k1_t
		gt = k2_t-2*k1_t+1
		#print(t, k1_t, k2_t, gt, mt, [mt-t >= 2*gt+2, mt-t < (n-t)/2])
		if t==k-1:
			#pre_check = ((mt-t) in [2*gt+2..floor_mod(n-t/2)]) or ((mt-t) in [ceil_mod((n-t)/2+2*gt-2)..n-t-4])
			pre_check = ((mt-t>=2*gt+2) and (mt-t<(n-t)/2))
			pre_g = gt
			pre_m = mt
		else:
			#cur_check = ((mt-t) in [2*gt+2..floor_mod(n-t/2)]) or ((mt-t) in [ceil_mod((n-t)/2+2*gt-2)..n-t-4])
			cur_check = ((mt-t>=2*gt+2) and (mt-t<(n-t)/2))
			cur_g = gt
			cur_m = mt
			if (cur_check and not pre_check) and (cur_g==pre_g) and (cur_m == pre_m):
				return t
				break
			else:
				pre_check = cur_check
				pre_g = cur_g
				pre_m = cur_m
		t -= 1
	if t==0:
		return false






def lifting_polynomials(G_orig, i):

	"""
	
	Computes polynomials required to embed a curve Y' into appropriate projective space.
	Input:
		G_orig - a generator matrix of the code (before shortening)
		i - number of a polynomial f_i to be generated
	Output:
		polynomial y[i]-f(y0..y[i-1]), s.t.
		f(Q)=0 for all points Q which projective coordinates are columns of matrix G_orig

	"""
	global F,p,v,R_S2c,y
	QPlaces = List_places(G_orig)
	#print(QPlaces[i])
	
	k=G_orig.nrows();
	
	Vars = Combinations([0..i-1],2).list()
	for LM_j in [0..i-1]:
		Vars.append([LM_j, LM_j])
	#print("%%%")
	#print(len(Vars))
	#print(binomial(i+1,2))
	#print(Vars)
	#print(len(QPlaces))
	#print(QPlaces)
	#print("%%%")
	LM = matrix(F, len(QPlaces), binomial(i+1,2)+i)
	LM_i=0
	for Q in QPlaces:
		LM_j = 0
		for Var in Vars:
			LM[LM_i, LM_j] = Q[Var[0]]*Q[Var[1]]
			LM_j = LM_j + 1
		LM_i = LM_i + 1
	LM_V = []
	for Q in QPlaces:
		LM_V.append(-Q[i])
		
	# Have to solve LM*X=LM_V
	#print(f"Number of equations={LM.nrows()}")
	#print(f"Number of vars={LM.ncols()}")
	SOL = LM.solve_right(vector(LM_V))
	#assert SOL.nrows()>=1, "No lifting polynomials of degree 2 found"
	f=0
	for LM_i in [0..len(Vars)-1]:
		f=f+SOL[LM_i]*y[Vars[LM_i][0]]*y[Vars[LM_i][1]]
	lifting_poly = y[i]-f
	#print("***")
	#print(SOL)
	#print(f)
	#Check that f(Q)=0 for all places Q
	#for Q in QPlaces:
	#	print(lifting_poly.subs(y0=Q[0], y1=Q[1], y2=Q[2], y3=Q[3], y4=Q[4], y5=Q[5], y6=Q[6]))
	#print("***")
	return lifting_poly
