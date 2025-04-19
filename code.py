import random

def set_n():
    while(n<8 and n>12):
        n = int(input())
        return n

def generate_ordered_pairs(n):
    ordered_pairs = []
    for i in n:
        ordered_pairs.append((random.random(),random.random()))
    return ordered_pairs
    
def get_h(ordered_pairs):
    h = []
    for i in range(ordered_pairs)-1:
        h.append(ordered_pairs[i+1][0]-ordered_pairs[i][0])
    return h

def get_left_side_matrix(h, n):
    left_side_matrix = [[0]]
    for i in range(n-3):
        left_side_matrix[i][i]=(h[i]+h[i+1])>>1 #check if it's *2 or /2
        left_side_matrix[i+1][i]=h[i+1]
        left_side_matrix[i][i+1]=h[i+1]
    left_side_matrix[n-2][n-2]=(h[n-2]+h[n-1])>>1
    return left_side_matrix


def get_right_side_matrix(h,ordered_pairs):
    right_side_matrix = []
    for i in range(ordered_pairs-2):
        right_side_matrix = 6*(ordered_pairs[i+1][1]/h[i+1] 
        - ordered_pairs[i][1]/h[i])
    return right_side_matrix


def SEL_solver():
    #Lineal ecuation system solver

def cubic_splin_interpolation(ordered_pairs):
    #include every function in here





