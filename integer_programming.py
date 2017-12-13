#!/usr/bin/python
#-*- coding: utf-8 -*-
# encoding: utf-8

from __future__ import print_function
from copy import deepcopy
import sys
import math

PRECISION = 4 #float precision

max_int_obj_v = float('-inf')
max_int_opt_x = None

def identity(n):
    return [[1 if i==j else 0 for j in range(n)] for i in range(n)]

def read_input(filename):
    with open(filename) as file:
        method = int(file.readline().strip()) #0 = cutting planes - 1 = branch and bound
        m = int(file.readline().strip())
        n = int(file.readline().strip())
        tableau = eval(file.readline().strip())

    tableau = [map(float, x) for x in tableau]
    tableau[0] = [-1.0*x for x in tableau[0]]
    return method, m, n, tableau

def printf(l):
    for s in l:
        print(*s)

class PL:
    
    def __init__ (self, tableau, m, n, ext_tableau=None):
        self.tableau = tableau
        
        if ext_tableau:
            self.ext_tableau = ext_tableau
        else:   
            self.ext_tableau = identity(m-1)
            self.ext_tableau.insert(0, [0 for i in range(m-1)]) #create extended tableau
        
        self.m = m #n of rows
        self.n = n #n of columns
        self.original_m = m
        self.original_n = n
        self.basis = []
        self.enforced = False

    def transform_fpi(self):
        curr_n = self.n
        for j in range(self.m - 1):
            for i in range(self.m):
                if j == i-1:
                    self.tableau[i].insert(curr_n-1,1.0)
                else:
                    self.tableau[i].insert(curr_n-1,0.0)
            curr_n = curr_n + 1
        self.basis = range(self.n-1,curr_n-1)
        self.n = curr_n

    def get_aux_PL(self):
        aux_tableau = deepcopy(self.tableau)
        aux_pl = PL(aux_tableau, self.m, self.n, deepcopy(self.ext_tableau))
        aux_pl.transform_fpi()

        aux_tableau[0] = [0.0 for x in range(aux_pl.n)]
        for i in range(1, self.m): #adjusting the value of c
            aux_tableau[0][-(i+1)] = 1.0

        return aux_pl

    def enforce_b(self):
        for i in range(1, self.m):
            if self.tableau[i][-1] < 0:
                self.tableau[i] = [-1*x for x in self.tableau[i]]
                self.ext_tableau[i] = [-1*x for x in self.ext_tableau[i]]
                self.enforced = True

    def del_last_columns(self, column_c):
        for i in range(self.m):
            for j in range(column_c):
                self.tableau[i].pop(-2)
        self.n = self.n - column_c

    def is_cannonic(self):
        for item in self.basis:
            for i in range(self.m):
                if i == (self.basis.index(item) + 1): #if it is the pivot line
                    if not self.tableau[i][item] == 1:
                        return False
                else:
                    if not self.tableau[i][item] == 0: #all but the pivot line has to be 0, including c
                        return False
        return True

    def cannonize(self, pivot):
        pivot_i = pivot[0]
        pivot_j = pivot[1]
        pivot_v = self.tableau[pivot_i][pivot_j]

        self.tableau[pivot_i] = [x/pivot_v for x in self.tableau[pivot_i]]
        self.ext_tableau[pivot_i] = [x/pivot_v for x in self.ext_tableau[pivot_i]]

        for i in range(self.m):
            if i == pivot_i:
                continue
            alpha = -1*self.tableau[i][pivot_j]
            self.tableau[i] = [alpha*x + y for x,y in zip(self.tableau[pivot_i], self.tableau[i])]
            self.ext_tableau[i] = [alpha*x + y for x,y in zip(self.ext_tableau[pivot_i], self.ext_tableau[i])]

    def primal_simplex(self, verbose=True):
        if not self.is_cannonic():
            for i in range(self.m - 1):
                pivot = (i+1,self.basis[i])
                self.cannonize(pivot)

        while True:
            if verbose:
                print(self.tableau)

            for i in range(self.n - 1): #tab[0] = c, -1 to exclude b
                if round(self.tableau[0][i], PRECISION) < 0:
                    k = i               #k is the variable thats gonna enter the basis
                    optimal = False
                    break
                optimal = True

            if optimal:
                opt_x = [0.0 for x in range(self.n - 1)]
                y = [self.ext_tableau[0][j] for j in range(self.m - 1)] #dual solution

                for i, item in enumerate(self.basis):
                    opt_x[item] = round(self.tableau[i+1][-1], PRECISION)
                
                obj_value = round(self.tableau[0][-1], PRECISION)
                return obj_value, opt_x, y

            aux = []
            for j in range(1, self.m):
                if not round(self.tableau[j][k], PRECISION) > 0: #Ajk <= 0
                    aux.append(float('inf')) #we are looking for the minimum value, so this wont be one
                else:
                    aux.append(self.tableau[j][-1]/self.tableau[j][k])

            t = min(aux)
            if t == float('inf'): #Aj <= 0 -> PL is unbounded
                unb_cert = [0 for x in range(column_max-1)]
                unb_cert[k] = 1
                
                for i, item in enumerate(self.basis):
                    unb_cert[item] = -1*self.tableau[i+1][k]

                obj_value = float('inf') #if the PL is unbounded we return inf as the obj_v
                return obj_value, unb_cert, None

            r = aux.index(t) #list.index already returns the smallest index contaning the value, so theres no need to explictly implement Bland's rule
            #line 0 is c
            pivot = (r+1, k) #the pivot is the index of the min value of aux
            self.basis[r] = k

            if not self.is_cannonic():
                self.cannonize(pivot)
    
    def dual_simplex (self, verbose=True):
        if not self.is_cannonic():
            for i in range(self.m - 1):
                pivot = (i+1,self.basis[i])
                self.cannonize(pivot)

        while True:
            if verbose:
                print(self.tableau)

            for i in range(1, self.m):
                if round(self.tableau[i][-1], PRECISION) < 0:
                    k = i              #basis[k] is the column thats gonna LEAVE the basis
                    optimal = False
                    break
                optimal = True

            if optimal:
                opt_x = [0.0 for x in range(self.n - 1)]
                y = [self.ext_tableau[0][j] for j in range(self.m - 1)] #dual solution

                for i, item in enumerate(self.basis):
                    opt_x[item] = round(self.tableau[i+1][-1], PRECISION)
                
                obj_value = round(self.tableau[0][-1], PRECISION)
                return obj_value, opt_x, y

            aux = []
            for j in range(self.n - 1):
                if round(self.tableau[k][j], PRECISION) >= 0: #Ajk <= 0
                    aux.append(float('inf')) #we are looking for the minimum value, so this wont be one
                else:
                    aux.append(round(-1.0*(self.tableau[0][j]/self.tableau[k][j]), PRECISION))

            t = min(aux)
            if t == float('inf'): #infisible
                return None, None, None
            
            r = aux.index(t) #list.index already returns the smallest index contaning the value, so theres no need to explictly implement Bland's rule
            
            pivot = (k, r) #the pivot is the index of the min value of aux
            self.basis[k-1] = r #-1 cause the first line is c, so it doesnt count

            if not self.is_cannonic():
                self.cannonize(pivot)

    def add_constraint(self, constraint_row, leq=True):
        if leq:
            constraint_row.insert(self.n-1, 1.0)
        else:
            constraint_row.insert(self.n-1, -1.0)
            constraint_row = ([-1.0*x for x in constraint_row])
        
        for i in range(self.m):
            self.tableau[i].insert(self.n-1, 0.0)
            self.ext_tableau[i].insert(self.n-1, 0.0)

        self.n = self.n + 1
        self.m = self.m + 1
        self.tableau.append(constraint_row)
        self.ext_tableau.append([0.0 for x in range(self.m)])
        self.basis.append(self.n-2)
            
    def cutting_planes(self):
        done = False
        selected = []

        while not done:
            for i, var in enumerate(self.basis, start=1):
                if var < self.original_n - 1: #if it is a actual variable of the original problem
                    selected_row = [math.floor(x) for x in self.tableau[i]] #dont add same restriction multiple times
                    if not round(self.tableau[i][-1], PRECISION).is_integer() and selected_row not in selected:
                        selected.append(selected_row)
                        done = False
                        break
                    else:
                        done = True

            if done == True:
                return obj_v, opt_x, y

            self.add_constraint(selected_row)
            obj_v, opt_x, y = self.dual_simplex()
            if obj_v is None: #infisible
                return None, None, None
    
    def branch_bound(self):
        global max_int_obj_v
        global max_int_opt_x

        for i, item in enumerate(self.basis):
            if item < self.original_n - 1: #if it is a actual variable of the original problem
                if not round(self.tableau[i+1][-1], PRECISION).is_integer():
                    var = item
                    b_value = self.tableau[i+1][-1]
                    break

        const = [0.0 for i in range(self.n)]
        const[var] = 1.0

        l_pl = deepcopy(self)
        const[-1] = math.floor(b_value)
        l_pl.add_constraint(deepcopy(const))
        
        r_pl = deepcopy(self)
        const[-1] = math.ceil(b_value)
        r_pl.add_constraint(deepcopy(const), leq=False)

        l_obj_v, l_opt_x, __ = l_pl.dual_simplex()
        if l_obj_v is not None:
            l_all_integer = True
            for x in l_opt_x[:self.original_n - 1]:
                if not round(x, PRECISION).is_integer():
                    l_all_integer = False
                    break
            
            if l_obj_v > max_int_obj_v:
                if l_all_integer:
                    max_int_obj_v = l_obj_v
                    max_int_opt_x = deepcopy(l_opt_x)
                else:
                    l_pl.branch_bound()
   
        r_obj_v, r_opt_x, __ = r_pl.dual_simplex()
        if r_obj_v is not None: #if pl is feasible
            r_all_integer = True
            for x in r_opt_x[:self.original_n - 1]:
                if not round(x, PRECISION).is_integer():
                    r_all_integer = False
                    break

            if r_obj_v > max_int_obj_v:
                if r_all_integer:
                    max_int_obj_v = r_obj_v
                    max_int_opt_x = deepcopy(r_opt_x)
                else:
                    r_pl.branch_bound()
        return 


def main():
    method, m, n, tableau = read_input(sys.argv[1])
   
    pl = PL(tableau, m+1, n+1)
    pl.transform_fpi()
    pl.enforce_b()
    
    aux_pl = pl.get_aux_PL()
    obj_v, opt_x, y = aux_pl.primal_simplex(verbose=False)
    
    if round(obj_v, PRECISION) < 0:
        print(u'PI inviável')
        return

    if pl.enforced:
        aux_pl.del_last_columns(m)
        aux_pl.tableau[0] = deepcopy(pl.tableau[0])
        pl = aux_pl

    obj_v, opt_x, y = pl.primal_simplex()
    if obj_v == float('inf'):
        print('PL ilimitada, aqui está um certificado ' + str(opt_x[:n]))
        return

    all_integer = True
    for i in range(m+1):
        if not round(opt_x[i], PRECISION).is_integer():
            all_integer = False
            break

    if not all_integer:
        if method == 0:
            v, x, __= pl.cutting_planes()
        else:
            pl.branch_bound()
            v = max_int_obj_v
            x = max_int_opt_x

        if v is not None:
            print('\nSolução ótima xi =', x[:n], 'com valor objetivo', v, \
            '. A solução ótima da relaxão linear é x =', opt_x[:n], ', com valor objetivo', obj_v, 'e certificado', y)
        else:
            print(u'PI inviável')
    else:
        print('\nSolução ótima xi =', opt_x[:n], 'com valor objetivo', obj_v, \
            '. A solução ótima da relaxão linear é x =', opt_x[:n], ', com valor objetivo', obj_v, 'e certificado', y)


if __name__ == '__main__':
    main()