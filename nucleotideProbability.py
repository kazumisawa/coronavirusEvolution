import math
import numpy as np

def Rmatrix(a, h, t):
    # transpose on June 7th
    result = np.array( [a] * 16 ).reshape( (4,4) )
    result[0,1] = h
    result[0,0] = -2*a -h
    result[1,1] = result[2,2] = result[3,3] = -3*a
    return result


def Tmatrix(a, h, t):
    # transpose on June 7th
    # unchanged by transpose because it is a diagonal matrix
    result = np.zeros( (4,4) )
    b = 3*a + h
    result[1,1] = math.exp( -b* t )
    result[2,2] = result[3,3] = math.exp( -4*a*t )
    return result


def Qmatrix(a, h, t):
    # transpose on June 7th
    result = np.array( [
		[4*a,	2*(h+a),	h+3*a,	h+3*a],
		[1,	-1,	0,	0],
		[0,	1,	0,	-1],
		[0,	0,	1,	-1]    
    ])
    return result

def Qprimematrix(a, h, t):
    # transpose on June 7th
    result = np.array( [
		[1,	4*h+8*a,	2*h+6*a,	-h-3*a],
		[1,	-4*a,	2*h+6*a,	-h-3*a],
		[1,	-4*a,	-2*h-6*a,	3*h+9*a],
		[1,	-4*a,	-2*h-6*a,	-h-3*a]
    ])
    return result

def Pmatrix(a, h, t):
    # transpose on June 7th
    b = 3*a + h
    tmp1 = np.dot( Qprimematrix(a, h ,t), Tmatrix(a, h ,t) )
    tmp2 = np.dot( tmp1, Qmatrix(a, h, t) )
    result = tmp2 / b / 4
    return result

def diffTa(a, h, t):
    # unchanged by transpose because it is a diagonal matrix
    result = np.zeros( (4,4) )
    b = 3*a + h
    result[1,1] = -3 * t * math.exp( -b* t )
    result[2,2] = result[3,3] = -4 * t * math.exp( -4*a*t )
    return result

def diffQa(a, h, t):
    # transpose on June 7th
    result = np.array( [
		[4,	2,	3,	3],
		[0,	0,	0,	0],
		[0,	0,	0,	0],
		[0,	0,	0,	0]
        ])
    return result


def diffQprimea(a, h, t):
    # transpose on June 7th
    result = np.array( [
		[0,	8,	6,	-3],
		[0,	-4,	6,	-3],
		[0,	-4,	-6,	9],
		[0,	-4,	-6,	-3]
        ])
    return result

def diffTh(a, h, t):
    result = np.zeros( (4,4) )
    b = 3*a + h
    result[1,1] = - t * math.exp( -b* t )
    return result

def diffQh(a, h, t):
    # transpose on June 7th
    result = np.array( [
		[0,	2,	1,	1],
		[0,	0,	0,	0],
		[0,	0,	0,	0],
		[0,	0,	0,	0]
        ])
    return result

def diffQprimeh(a, h, t):
    # transpose on June 7th
    result = np.array( [
		[0,	4,	2,	-1],
		[0,	0,	2,	-1],
		[0,	0,	-2,	3],
		[0,	0,	-2,	-1]
        ])
    return result

def D(a, h, t):
    # D: 1/(4*b), a scalar
    b = 3*a + h
    return 1/(4*b)

def diffDa(a, h, t):
    b = 3*a + h
    return -3/(4*b*b)

def diffDh(a, h, t):
    b = 3*a + h
    return -1/(4*b*b)

