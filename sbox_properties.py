def absolute_indicator(matrix,M,N):
	m_shift = 1 << M
	n_shift = 1 << N
	#input_array = [0 for i in range(len(matrix)*len(matrix[0]))]
	#input_array = [0 for i in range(M*(2**N))]
	k=0
	'''for i in range(len(matrix)):
		for j in range(len(matrix[i])):
			input_array[k] = matrix[i][j]
			k+=1'''
	input_array = matrix
	truth_table = [[0 for i in range(n_shift-1)] for j in range(m_shift)]

	for i in range(m_shift):
		for j in range(N):
			truth_table[i][j] = input_array[i] >> (N - j - 1) & 0x01

	ll = [[0 for i in range(n_shift-1)] for j in range(m_shift)]

	if N==1:
		for i in range(m_shift):
			ll[i][0] = truth_table[i][0]
	else:		
		for i in range(1,n_shift):
			for j in range(N):
				if (i >> j & 0x01):
					for k in range(m_shift):
						ll[k][i-1] = ll[k][i-1] ^ truth_table[k][j]

	wt = [[0 for i in range(n_shift-1)] for j in range(m_shift)]			

	for z in range(n_shift-1):
		for i in range(m_shift):
			wt[i][z] = 1 if ll[i][z]==0 else -1

		for i in range(1,M+1):
			m = 1<<i
			halfm = int(m/2)
			for r in range(0,m_shift,m):
				t1 = r
				t2 = r+halfm
				for j in range(halfm):
					a = wt[t1][z]
					b = wt[t2][z]
					wt[t1][z] = a + b
					wt[t2][z] = a - b
					t1+=1
					t2+=1
	
	ac = [[0 for i in range(n_shift-1)] for j in range(m_shift)]

	for z in range(n_shift-1):
		for i in range(m_shift):
			ac[i][z] = -1* wt[i][z]*wt[i][z]
		for i in range(1,M+1):
			m  = 1 << i
			halfm = int(m/2)
			for r in range(0,m_shift,m):
				t1 = r
				t2 = r+halfm
				for j in range(halfm):
					a = ac[t1][z]
					b = ac[t2][z]
					ac[t1][z] = a + b
					ac[t2][z] = a - b
					t1+=1
					t2+=1
		for i in range(m_shift):
			ac[i][z] /= (1 << M) * (-1)

	#print(truth_table,"\n",ll,"\n",wt,"\n",ac,"\n")		

	maxm=0
	for j in range(n_shift-1):
		temp = abs(ac[1][j])
		for i in range(2,m_shift):
			temp2 = abs(ac[i][j])
			if temp2>temp:
				temp=temp2
		if(temp>maxm):
			maxm=temp
	
	return maxm
def differentialUniformity(S,k,m):
    c = 0
    delta = 0
    n = 2**k
    m = 2**m
    for alpha in range(1, n):
        for beta in range(m):
            c = 0
            for z in range(n):
                if ((S[(z ^ alpha)] ^ S[z]) == beta):
                    c = c + 1
            if (c > delta):
                    delta = c
            
    return delta

# program returns array specifying nonlinearity of each component function
# takes input as a 1-D array as S-Box.
###########################################################################
#finds Walsh Transform of truth table(f)
def fwt(f):  # f is a Boolean function represented as a TT(0/1) of length 2^n
    import math
    wf = []
    for x in f:
        if x == 0:
            wf.append(1)
        else:  
            wf.append(-1)
    order = len(f)  # order = 2^n
    n = int(math.log(order, 2))
    size = int(math.floor(order / 2))
    while size >= 1:
        left = 0
        while left < order - 1:
            for p in range(int(size)):
                right = left + int(size)
                a = wf[left]
                b = wf[right]
                wf[left] = a + b
                wf[right] = a - b
                left = left + 1
            left = right + 1
        size = int(math.floor(size / 2))
    # print"\tWalsh transform of function's truth table is",
    # print wf
    return wf

############################################################################
#finds non-linearity of an n-variable boolean function 'f'
def bf_nonlinearity(f, n):
    import math
    fw = fwt(f)
    #find modulus of each element in Walsh transform
    for i in range(len(fw)):
        fw[i] = abs(fw[i])
    # nonlinearity from the Walsh transform
    x = ((2 ** (n - 1)) - (max(fw) / 2))
    # print"\tNL of function is",
    # print x
    return x

##############################################################################
#converts num to binary form (no of bits in binary representation = length) 
def binary(num, length):
    binary_string_list = list(format(num, '0{}b'.format(length)))
    #print("num,binary String List",num,binary_string_list)
    return [int(digit) for digit in binary_string_list]

##############################################################################
#returns array of NL of each component function for sbox S



#Optimized NL rewritten
def non_linearity(S,n,m):  
    order = 2**n
    #print("Order is ",order)
    #print("n is ",n)
    nl_array = []  # nl_array[] stores calculated NL for each function yi
    binary_value_array= []
    for index in range(order):
        binary_value_array.append(binary(S[index], m))
    columns = zip(*binary_value_array) # Transpose the binary matrix
    columns = list(columns) # change zip object into list
    for bitno in range(m):
        bfnl = bf_nonlinearity(columns[bitno], n)
        nl_array.append(bfnl)
        
    return sum(nl_array)/len(nl_array)

