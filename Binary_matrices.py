### Classical problems modelled with binary matrices ###
# _ Island search

# Author : Mathis Antonetti
# Date : 2023/02

### PACKAGES ###
import numpy as np
import matplotlib.pyplot as plt

### TOOLS ###
def connected(matrix, i, j, di, dj, k, L, n, m):
    # Parameters :
# (matrix, n, m) : a matrix (numpy array expected) with n rows and m columns
# (i,j) : index of the node (matrix element)
# (i+di, j+dj) : index of the neighboor
# k : current layer
# L : list of indices
    # returns True if the neighboor is in L
    
    l = len(L)
    tc = min(l, 2*(n-2*k)+2*(m-2*(k+1)))
    for u in range(tc):
        if((i+di == L[l-1-u][0]) and (j+dj == L[l-1-u][1])):
            return True
    return False

def searchBorder(matrix, n , m):
    # Parameters :
# (matrix, n, m) : a matrix (numpy array expected) with n rows and m columns
    # returns the list of indices of the 1s in the border
    L = []
    for i in range(n):
        if(matrix[i,0] == 1):
            L.append([i, 0])
        if(matrix[i, m-1] == 1):
            L.append([i, m-1])
    for j in range(1,m-1):
        if(matrix[0,j] == 1):
            L.append([0,j])
        if(matrix[n-1,j] == 1):
            L.append([n-1,j])
    return L

def searchBorderIslands(matrix, n , m, indices, k, isLast, dir):
    # Parameters
# (matrix, n, m) : a matrix (numpy array expected) with n rows and m columns
# indices : list of indices of 1s that are connected to a 1 in the border
# k : current layer
    # returns (mod_list, mod_count) where : 
# _ mod_list is the current list of indices of 1s that are connected to a 1 in the border
# _ mod_count is the count of modifications needed for the layer k
    mod_list = indices.copy() 
    mod_count = 0

    current_mod = [] # current list of 1s
    status = False # True if an element of current_mod is connected to the border
    for side in range(4):

        # search path coordinator
        u, i0, j0, N, switch = 1, 0, 0, 0, 0
        if(side == 0): # left
            u, i0, j0, N, switch = 1, k, k, n-2*k, 0
        elif(side == 1): # bottom
            u, i0, j0, N, switch = 1, n-k-1, k, m-2*k, 1
        elif(side == 2): # right
            u, i0, j0, N, switch = -1, n-k-1, m-k-1, n-2*k, 0
        elif(side == 3): # top
            u, i0, j0, N, switch = -1, k, m-k-1, m-2*k, 1

        # search modifications along the path
        for addit in range(N):
            i, j, di, dj = i0+(1-switch)*u*addit, j0+switch*u*addit, switch*u, (1 - switch)*(-u)
            if(matrix[i,j] == 1):
                current_mod.append([i,j])
                if(connected(matrix, i, j, dir*di, dir*dj, k, mod_list, n, m)):
                    status = True
                if(isLast and (1 <= addit) and (addit <= N-2) and (matrix[i-di,j-dj] == 1)): # recuperate missed points for interior nodes
                    current_mod.append([i-di,j-dj])
            else:
                if(status == True):
                    mod_count += len(current_mod)
                    mod_list = np.concatenate((mod_list, current_mod))
                    status = False
                current_mod = []

    if(status == True):
        mod_count += len(current_mod)
        mod_list = np.concatenate((mod_list, current_mod))

    return mod_list, mod_count

# return binary matrix without islands cut at index K
def removeIslandsCut(matrix, K):
    # Parameters
# matrix : a binary matrix (as numpy array)
    # Returns the modified matrix without islands
    m, n = len(matrix[0]), len(matrix)
    matrix2 = np.zeros((n, m)) # initialization
    indices = searchBorder(matrix, n, m) # indices of non-island 1s
    k = 1 # incrementation index (also the layer)
    c = len(indices) # adds count for the current k

    # search the non-island 1s
    while(c >= 1 and k < min(K+1, min(n//2 + (n%2), m//2 + (m%2)))):
        indices, c = searchBorderIslands(matrix, n, m, indices, k, k == min(K+1, min(n//2 + (n%2), m//2 + (m%2)))-1, 1)
        k += 1

    k -= 1
    while(k >= 1):
        indices, c = searchBorderIslands(matrix, n, m, indices, k, False, -1)
        k -= 1

    # put non-island 1s in the new matrix
    for i in range(len(indices)):
        matrix2[indices[i][0], indices[i][1]] = 1

    return matrix2

# return binary matrix without islands
def removeIslands(matrix):
    # Parameters
# matrix : a binary matrix (as numpy array)
    # Returns the modified matrix without islands
    return removeIslandsCut(matrix, max(len(matrix[:,0]), len(matrix[0,:])))

### TESTS ###

# mid tools test (k-th iteration step)
def run_test(num_test, K):
    print("Running the tools searchBorder and searchBorderIslands on a ")
    if(num_test == 0):
        print("tiny test with the ")
        Mat = np.array([[1, 0, 1, 0], [0, 1, 1, 1], [0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 1, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [1, 0, 0, 0]])
    elif(num_test == 1):
        print("small test with the ")
        Mat = np.array([[0, 0, 0, 0, 0, 1, 0, 0, 1], [0, 0, 1, 1, 0, 1, 0, 0, 0], [0, 1, 1, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1, 1], [0, 0, 1, 1, 1, 0, 0, 0, 0], [0, 1, 0, 0, 1, 0, 1, 0, 0], [1, 1, 0, 0, 0, 0, 1, 0, 1], [0, 0, 1, 1, 1, 0, 0, 0, 1]])
    print("matrix : \n")
    print(Mat)
    print("\n Matrix without islands : \n")
    print(removeIslandsCut(Mat, K))
    print("\n")


# functionality test case
def run_globaltestcase(num_test):
    M = np.zeros((2, 2))
    print("Running the tool removeIslands on a ")
    if(num_test == 0):
        print("tiny test with the ")
        M = np.array([[0, 0, 0, 0, 0], [0, 0, 1, 1, 0], [0, 1, 1, 0, 0], [0, 0, 0, 0, 0], [0, 0, 1, 1, 1], [0, 1, 0, 0, 1], [1, 1, 0, 0, 0], [0, 0, 1, 1, 1]])
    elif(num_test == 1):
        print("simple and small test with the ")
        M = np.array([[0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 1, 1, 0, 1, 0, 0], [0, 1, 1, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 1, 1, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 1, 0], [1, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 1, 1, 0, 0, 0]])
    elif(num_test == 2):
        print("complex and small test with the ")
        M = np.array([[0, 0, 0, 0, 0, 1, 0, 0, 1], [0, 0, 1, 1, 0, 1, 0, 0, 0], [0, 1, 1, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 1, 1], [0, 0, 1, 1, 1, 0, 1, 1, 0], [0, 1, 0, 0, 1, 0, 1, 0, 0], [1, 1, 0, 0, 0, 1, 1, 0, 1], [0, 0, 1, 1, 1, 0, 0, 0, 1]])
    print("matrix tested : \n")
    print(M)
    print("\n Matrix without islands : \n")
    print(removeIslands(M))
    print("\n")

#run_test(1, 2)

run_globaltestcase(2)
run_globaltestcase(1)
run_globaltestcase(0)