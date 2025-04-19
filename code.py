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
    for i in len(ordered_pairs)-1:
        h.append(ordered_pairs[i+1][0]-ordered_pairs[i][0])
    return h

def get_left_side_matrix(h, n):
    left_side_matrix = [[0]]
    for i in n-3:
        left_side_matrix[i][i]=(h[i]+h[i+1])>>1 #check if it's *2 or /2
        left_side_matrix[i+1][i]=h[i+1]
        left_side_matrix[i][i+1]=h[i+1]
    left_side_matrix[n-2][n-2]=(h[n-2]+h[n-1])>>1
    return left_side_matrix


def get_right_side_matrix(h,ordered_pairs, n):
    right_side_matrix = []
    for i in n-2:
        right_side_matrix = 6*(ordered_pairs[i+1][1]/h[i+1] 
        - ordered_pairs[i][1]/h[i])
    return right_side_matrix


def get_determinant_n(matrix, n):
    sum, substraction = 0
    for i in n:
        for j in n:
            #verificar si len(matrix) es n o n-1
            sum+=matrix[(i+j)%n][(i+j)%n]
            res+=matrix[n-1-j][n-1-j]
    return (sum-substraction)

def putting_answer_in_varible(right_matrix,left_matrix, pos_variable, n):
    new_matrix = left_matrix
    for i in n:
        #asumiendo que el primer indice refiere a fila y segundo columna
        new_matrix[i][pos_variable]=right_matrix[i]

def crammer(right_matrix, left_matrix, n):
        determinant = []
        determinant[0] = get_determinant_n(right_matrix)
        for i in range(1,n):
            determinant[i] = get_determinant_n(putting_answer_in_varible(right_matrix,left_matrix,i-1))
        return [determinant[i]/2 for determinant[i] in range(1,len(determinant))]



def cubic_splin_interpolation(ordered_pairs, n):
    h = get_h(ordered_pairs)
    left_matrix = get_left_side_matrix(h,n)
    right_matrix = get_right_side_matrix(h,ordered_pairs,n)
    values = crammer(right_matrix,left_matrix,n)


#include a function for function graphics

def main():
    n = set_n()
    ordered_pairs = generate_ordered_pairs(n)
    cubic_splin_interpolation(ordered_pairs, n)




