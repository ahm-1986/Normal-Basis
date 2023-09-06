#!/local/bin/python
var('b');
var('x');

IR=[x^8 + x^4 + x^3 + x + 1,
x^8 + x^4 + x^3 + x^2 + 1,
x^8 + x^5 + x^3 + x + 1,
x^8 + x^5 + x^3 + x^2 + 1,
x^8 + x^5 + x^4 + x^3 + 1,
x^8 + x^5 + x^4 + x^3 + x^2 + x + 1,
x^8 + x^6 + x^3 + x^2 + 1,
x^8 + x^6 + x^4 + x^3 + x^2 + x + 1,
x^8 + x^6 + x^5 + x + 1,
x^8 + x^6 + x^5 + x^2 + 1,
x^8 + x^6 + x^5 + x^3 + 1,
x^8 + x^6 + x^5 + x^4 + 1,
x^8 + x^6 + x^5 + x^4 + x^2 + x + 1,
x^8 + x^6 + x^5 + x^4 + x^3 + x + 1,
x^8 + x^7 + x^2 + x + 1,
x^8 + x^7 + x^3 + x + 1,
x^8 + x^7 + x^3 + x^2 + 1,
x^8 + x^7 + x^4 + x^3 + x^2 + x + 1,
x^8 + x^7 + x^5 + x + 1,
x^8 + x^7 + x^5 + x^3 + 1,
x^8 + x^7 + x^5 + x^4 + 1,
x^8 + x^7 + x^5 + x^4 + x^3 + x^2 + 1,
x^8 + x^7 + x^6 + x + 1,
x^8 + x^7 + x^6 + x^3 + x^2 + x + 1,
x^8 + x^7 + x^6 + x^4 + x^2 + x + 1,
x^8 + x^7 + x^6 + x^4 + x^3 + x^2 + 1,
x^8 + x^7 + x^6 + x^5 + x^2 + x + 1,
x^8 + x^7 + x^6 + x^5 + x^4 + x + 1,
x^8 + x^7 + x^6 + x^5 + x^4 + x^2 + 1,
x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + 1]

IR4=[
    x^4+x+1,x^4+x^3+1,x^4+x^3+x^2+x+1
]

IR4_power=[
    [4,1,0],
    [4,3,0],
    [4,3,2,1,0]]

IR_power=[
[8 , 4 , 3 , 1 , 0],
[8 , 4 , 3 , 2 , 0],
[8 , 5 , 3 , 1 , 0],
[8 , 5 , 3 , 2 , 0],
[8 , 5 , 4 , 3 , 0],
[8 , 5 , 4 , 3 , 2 , 1 , 0],
[8 , 6 , 3 , 2 , 0],
[8 , 6 , 4 , 3 , 2 , 1 , 0],
[8 , 6 , 5 , 1 , 0],
[8 , 6 , 5 , 2 , 0],
[8 , 6 , 5 , 3 , 0],
[8 , 6 , 5 , 4 , 0],
[8 , 6 , 5 , 4 , 2 , 1 , 0],
[8 , 6 , 5 , 4 , 3 , 1 , 0],
[8 , 7 , 2 , 1 , 0],
[8 , 7 , 3 , 1 , 0],
[8 , 7 , 3 , 2 , 0],
[8 , 7 , 4 , 3 , 2 , 1 , 0],
[8 , 7 , 5 , 1 , 0],
[8 , 7 , 5 , 3 , 0],
[8 , 7 , 5 , 4 , 0],
[8 , 7 , 5 , 4 , 3 , 2 , 0],
[8 , 7 , 6 , 1 , 0],
[8 , 7 , 6 , 3 , 2 , 1 , 0],
[8 , 7 , 6 , 4 , 2 , 1 , 0],
[8 , 7 , 6 , 4 , 3 , 2 , 0],
[8 , 7 , 6 , 5 , 2 , 1 , 0],
[8 , 7 , 6 , 5 , 4 , 1 , 0],
[8 , 7 , 6 , 5 , 4 , 2 , 0],
[8 , 7 , 6 , 5 , 4 , 3 , 0]]

IR2_power=[[2,1,0]]


AF_AES=[[1,1,1,1,1,0,0,0],
		[0,1,1,1,1,1,0,0],
		[0,0,1,1,1,1,1,0],
		[0,0,0,1,1,1,1,1],
		[1,0,0,0,1,1,1,1],
		[1,1,0,0,0,1,1,1],
		[1,1,1,0,0,0,1,1],
		[1,1,1,1,0,0,0,1]]


# AF_CLEFIA=[[0,0,0,1,1,0,0,0],
#            [0,1,0,1,0,0,0,1],
#            [0,0,0,0,0,0,0,1],
#            [0,0,0,0,0,1,1,0],
#            [0,1,1,0,0,1,0,1],
#            [0,1,0,1,1,1,0,0],
#            [0,1,1,0,0,0,0,0],
#            [1,0,0,0,0,0,0,1]]


AF_CLEFIA=[[0,0,0,0,1,0,1,0],
           [0,1,0,0,0,0,0,1],
           [0,1,0,1,1,0,0,0],
           [0,0,1,0,0,0,0,0],
           [0,0,1,1,0,0,0,0],
           [0,0,0,0,0,0,1,0],
           [1,0,0,1,0,0,0,0],
           [0,1,0,0,0,1,0,0]] 





poly_primitive=[3,
              2,
              2,
              2,
              3,
              3,
              2,
              2,
              2,
              2,
              2,
              2,
              3,
              9,
              2,
              6,
              2,
              3,
              3,
              2,
              6,
              7,
              2,
              2,
              7,
              6,
              2,
              6,
              2,
              3]
