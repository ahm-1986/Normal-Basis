#!/local/bin/python




'''
This code is to generate the change basis matrix for inverse sbox such AES , Clefia.
The current used matrix is for CLEFIA as PoC. The 4 output arrays are used in sbox_NB_NW.c 

'''
#from IRlist import IR
import random


sage.repl.load.load("const.sage", globals())


def sort_lists(A,B):
    min_B=min(B,default="EMPTY")
   

    if(min_B!="EMPTY" and len(B)==8 and len(A)==8):
        new_A=[A[B.index(min_B)]]
        new_B=[min_B]

        for i in range(0,len(B)-1):
            min_B=min_B*2%255
            new_A.append(A[B.index(min_B)])
            new_B.append(min_B)


        return new_A,new_B
    else:
        return A,B



def evaluate_poly(E,i):
    sum=0
    for W in E:
        sum=sum+(i**W)
    return sum


def find_IR_expression(poly):
    sum=0;
    
    for a in poly:
        sum=sum+b^a 
    return sum

'''
The function below generate the primitive elements of each IR polynomial in refrence to IR0
'''
def find_primitive_elements(poly_num,IR_power,base_IR):
    PE_list=[]
    for i in range(0,256):
        D.<b> = GF(2^8, modulus=IR[base_IR])
        #b=b+1
        e=D.fetch_int(i)
        if(evaluate_poly(IR_power[poly_num],e)==0):
            PE_list.append(i)

    return PE_list



def find_IR_primitve(poly_num,IR_power,base_IR):
    #print(IR[poly_num])
    
    D.<b> = GF(2^8, modulus=IR[base_IR])   
    PE=find_primitive_elements(poly_num,IR_power,base_IR)

    PList=[]

    #ref=IR_expression(poly_primitive[base_IR])
    #print("ref:",ref,"  ")
    for P in  PE:
       
        T=1

        ref=D.fetch_int(poly_primitive[base_IR])
        #ref=b+1 # this is the primitive root and should be changed.
        #ref=b # this is the primitive root and should be changed.
        #ref=b^3+1 # this is the primitive root and should be changed. 
        for A in range(0,255):
           
            if(T==D.fetch_int(P)):
                PList.append(A)

            T=T*ref
              
    PE,PList=sort_lists(PE,PList)
   
    return PE,PList



o = open('stuff.txt','w')
#log table generation
def generate_logtable(base_IR):
    logB=[0 for i in range(0,256)]
    for i in range(0,len(IR)):
        print("----------")
        A,B=find_IR_primitve(i,IR_power,base_IR)
        print(A)
        print(IR[i],"|",B)
        
        assert(len(A)==8)
        assert(len(B)==8)
        for j in range(0,len(A)):
            logB[A[j]]=B[j]

    for i in range(0,len(IR4_power)):
        print("----------")
        A,B=find_IR_primitve(i,IR4_power,base_IR)
        print(A)
        print(B)
        assert(len(A)==4)
        assert(len(B)==4)
        for j in range(0,len(A)):
            logB[A[j]]=B[j]

    
    for i in range(0,len(IR2_power)):
        print("----------")
        A,B=find_IR_primitve(i,IR2_power,base_IR)
        print(A)
        print(B)
        assert(len(A)==2)
        assert(len(B)==2)
        for j in range(0,len(A)):
            logB[A[j]]=B[j]

    return logB


#------------------------

def concatenate_vertical(T):
	bin2hex=[]
	for i in range(0,8):
		b2x=""
		for j in range(0,8):
			
			b2x=b2x+str(T[j][i])
		
		bin2hex.append('0x%02x' % int(b2x,2))
	

	return bin2hex

def concatenate_horizontal(T):
	bin2hex=[]
	for i in range(0,8):
		b2x=""
		for j in range(0,8):
			
			b2x=b2x+str(T[i][j])
		
		bin2hex.append('0x%02x' % int(b2x,2))
	

	return bin2hex


def hex_to_bin(M):
	L=[[0 for i in range(8)] for i in range(0,8)]
	
	for i in range(0,8):
		
		#L.append([0*8])
		for j in range(0,8):
			#print(j,"P")
			#print(M[i],(M[i]>>j)&0x1)
			L[i][7-j]= (M[i]>>j)&0x1
	#print(L)
	return L

def change_matrix(Y_1, Y_16, Z_1, Z_4, w_1, w_2):
    L = [
        int((w_2 + Z_4 + Y_16) % 255),
        int((w_1 + Z_4 + Y_16) % 255),
        int((w_2 + Z_1 + Y_16) % 255),
        int((w_1 + Z_1 + Y_16) % 255),
        int((w_2 + Z_4 + Y_1) % 255),
        int((w_1 + Z_4 + Y_1) % 255),
        int((w_2 + Z_1 + Y_1) % 255),
        int((w_1 + Z_1 + Y_1) % 255)
    ]

    return L


def generate_antilogtable(T):
    ALT = [i for i in range(len(T)-1)]
    for idx, value in enumerate(T):
       
        ALT[value] = idx
    #ALT.insert(0,0)
    return ALT

def wzy_N_V(Y_1, Y_16, Z_1, Z_4, w_1, w_2, LT, ALT):
    V=(Y_1+Y_16)%255
    N=(Z_1+Z_4)%255
    W=(w_1+w_2)%255
    print("v:",V)
    print("N:",N)
    print("W:",W)

    print(LT)
    print("----------")
    print(ALT)
    if (w_1==1):
         w_1=w_2
    
    if(N==w_1):
        N=w_1
        print("N=Ω")
    else:
        N=2*w_1
        print("N=Ω^2")
        
    #print(N_w,N_w2)
    if((N+Z_4)%255 == V):
        print("C:N    D:0")
        print("#define OPTION 0")
    elif((2*N+Z_4)%255 == V):
        print("C:N^2    D:0")
        print("#define OPTION 2")
    
    elif((N+Z_1)%255 == V):
        print("C:0    D:N")
        print("#define OPTION 1")
    elif((2*N+Z_1)%255 == V):
        print("C:0    D:N^2")
        print("#define OPTION 3")


    elif(LT[ALT[(2*N+Z_4)%255]^^ALT[(Z_1)%255]] == V):
        print("C:N^2    D:1")
        print("#define OPTION 6")
    
    elif(LT[ALT[(N+Z_4)%255]^^ALT[(Z_1)%255]] == V):
        print("C:N    D:1")
        print("#define OPTION 4")

    elif(LT[ALT[(Z_4)%255]^^ALT[(2*N+Z_1)%255]] == V):
        print("C:1    D:N^2")
        print("#define OPTION 7")
    
    elif(LT[ALT[(Z_4)%255]^^ALT[(N+Z_1)%255]] == V):
        print("C:1    D:N")
        print("#define OPTION 5")


    N=w_1
    S=(2*N+Z_1) %255
    print(N,S,V)



def printout(S,L):
    print("static int   ",S,"[8] = {",L[0],",",L[1],",",L[2],",",L[3],",",L[4],",",L[5],",",L[6],",",L[7],"};")


def matrix_generation(X):
    print("------------")
    
    #print("X2A")
    XB=hex_to_bin(X)
    #print(XB)
    X2A=concatenate_horizontal(XB)
    printout("X2A",X2A)
    #print(concatenate_horizontal(XB))
    M=matrix(GF(2),XB)
    X_1=M.inverse()

    A2X=concatenate_horizontal(X_1)
    printout("A2X",A2X)

    CM=matrix(GF(2),AF_CLEFIA) # AES sbox selected  AF_AES, to select CLEFIA: AF_CLEFIA
    MX=M*CM.transpose()

    

    X2S=concatenate_horizontal(MX)
    printout("X2S",X2S)
    
    
    MX=matrix(GF(2),MX)

    MX_1=MX.inverse()
    S2X=concatenate_horizontal(MX_1)
    printout("S2X",S2X)
    print("#####")

    #N=w_1
    #S=(2*N+Z_1) %255
    #print(N,S,V)





#Root_IR=2
def calculate_basis(Root_IR,Y_1,Y_16,Z_1,Z_4,W_1,W_2):
    #F.<a> = GF(2^8, modulus=x^8+x^4+x^3+x+1)
    F.<a> = GF(2^8, modulus=IR[Root_IR])
    #IR_primitve(2,IR4_power,Root_IR) # base_IR : x^8 + x^4 + x^3 + x + 1
    #IR_primitve(2,IR2_power,Root_IR)


    LT= generate_logtable(Root_IR)
    ALT=generate_antilogtable(LT)


    CM=change_matrix(Y_1,Y_16,Z_1,Z_4,W_1,W_2)
    CMH=[]

    
    print("----------")
    for a in CM:
        #print(hex(LT.index(a)))
        CMH.append(int(LT.index(a)))

    #print(CMH)
    matrix_generation(CMH)
    wzy_N_V(Y_1,Y_16,Z_1,Z_4,W_1,W_2,LT,ALT)

'''

for i in range(10,20):
    print("########")
    calculate_basis(i)
'''

#Example 1 7 112 34 136 85 170 
def main():
    i = int(sys.argv[1])
    Y1 = int(sys.argv[2])
    Y16 = int(sys.argv[3])
    Z1 = int(sys.argv[4])
    Z4 = int(sys.argv[5])
    W1 = int(sys.argv[6])
    W2 = int(sys.argv[7])
    print(i)
    calculate_basis(i,Y1,Y16,Z1,Z4,W1,W2)

if __name__ == "__main__":
    main()



