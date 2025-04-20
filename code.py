import random

def set_n():
    n = 0
    while(n<8 or n>12):
        n = int(input())
        return n

def generate_ordered_pairs(n):
    ordered_pairs = []
    for i in range(n):
        ordered_pairs.append((random.random(),random.random()))
    return ordered_pairs
    
def get_h(ordered_pairs):
    h = []
    for i in range(len(ordered_pairs)-1):
        h.append(ordered_pairs[i+1][0]-ordered_pairs[i][0])
    return h

def get_left_side_matrix(h, n):
    left_side_matrix = [[0 for _ in range(n-2)] for _ in range(n-2)]
    for i in range(n-3):
        left_side_matrix[i][i]=(h[i]+h[i+1])*2  #>>1 check if it's *2 or /2
        left_side_matrix[i+1][i]=h[i+1]
        left_side_matrix[i][i+1]=h[i+1]
    left_side_matrix[i+1][i+1]=(h[i+1]+h[i+2])*2
    return left_side_matrix


def get_right_side_matrix(h,ordered_pairs, n):
    right_side_matrix = []
    for i in range(1,n-1):
        a1 = ordered_pairs[i+1][1]
        a2 = ordered_pairs[i][1]
        a3 = ordered_pairs[i-1][1]
        right_side_matrix.append(6*((a1-a2)/h[i] -
                                    (a2-a3)/h[i-1]))
    return right_side_matrix


def get_determinant_n(matrix, n):
    sum=0
    substraction = 0
    for i in range(n-2):
        for j in range(n-2):
            #verificar si len(matrix) es n o n-1
            index = (i+j)%n-2
            sum+=matrix[index][index]
            substraction+=matrix[n-3-j][n-3-j]
    return (sum-substraction)

def putting_answer_in_varible(right_matrix,left_matrix, pos_variable, n):
    new_matrix = left_matrix
    for i in range(n-2):
        #asumiendo que el primer indice refiere a fila y segundo columna
        new_matrix[i][pos_variable]=right_matrix[i]
    return new_matrix

def crammer(right_matrix, left_matrix, n):
        determinant = []
        determinant.append(get_determinant_n(left_matrix,n))
        for i in range(1,n-2):
            determinant.append(get_determinant_n(putting_answer_in_varible(right_matrix,left_matrix,i-1,n),n))
        return [determinant[i]/2 for determinant[i] in range(1,n-2)]



def cubic_splin_interpolation(ordered_pairs, n):
    h = get_h(ordered_pairs)
    left_matrix = get_left_side_matrix(h,n)
    right_matrix = get_right_side_matrix(h,ordered_pairs,n)
    values = crammer(right_matrix,left_matrix,n)
    return values

#include a function for function graphics

def main():
    n = set_n()
    #ordered_pairs = generate_ordered_pairs(n)
    ordered_pairs = [(0.1,10),(0.2,5),(0.5,2),(1,1),(2,0.5),(5,0.2),(10,0.1),]
    points = cubic_splin_interpolation(ordered_pairs, n)
    print(points)


if __name__ == '__main__':
    main()
