# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 12:02:29 2023

@author: Michael Ferguson
"""

import sys
from sympy import Symbol, symbols, symarray, oo, nonlinsolve
from sympy.combinatorics import Permutation
from sympy import init_printing
init_printing(perm_cyclic=False, pretty_print=False)

"""
Resources from custom packages:

"""
from chamber_package import element, chamber
from solutions_package import system_solution, frame_vector

"""
The symbol s represents the chosen uniformizing parameter for the field.

"""
s = Symbol('s')

"""
This defines the change of basis matrix using lists.

Notes:
======
1) This can be changed by the user, but it must be a square matrix consisting
   of elements with non-zero determinant. Requires more tweaking for matrices
   of other dimensions, but still fine to do.

2) This particular matrix was chosen at random from the family described in
   Theorem 5.2.1 from my thesis.

"""
Row_1 = [s**2+s+2,s,s+1]
Row_2 = [s,1+s,s]
Row_3 = [4+s,s,s+2]

CoB_Mat = [Row_1,Row_2,Row_3]
mat_dim = len(CoB_Mat)

"""
Elements representing the functions used by Algorithm A (see thesis). 

"""

x = symarray('x',mat_dim)

M = [element(sum([x[j]*CoB_Mat[i][j] for j in range(mat_dim)])) for i in range(mat_dim)]

"""
Symbols for our bases

"""

u_base = symarray('u',mat_dim)
v_base = symarray('v',mat_dim)

def printsw(foo,bar):
    """
    A convenience function.

    """
    if bar:
        print(f'{foo}\n')

def get_candidate(u_cham,v_cham,u_num,v_num,kl=oo,timeout=oo,show_work = False):
    """
    Returns a candidate vector in the form of a frame_vector from chamber_package
    
    Notes:
    ======
    get_candidate is step r for a given w_i in Algorithm A (see paper)
    
    u_num and v_num are used to keep track of lattice indices

    kl is the k value obtained from the previous call of this this function
    within get_w
    
    """
    sw = show_work
    cham_dim = v_cham.dim
    v_vert = v_cham.vertices[v_num]
    u_vert = u_cham.vertices[u_num]
    F_num = u_cham.perm((u_num-1)%cham_dim)

    DVR_symbols = symarray('f',cham_dim)
    generic_lattice_members = [DVR_symbols[i]*s**v_vert[i] for i in range(cham_dim)]
    Fsub = [F.subs(list(zip(x,generic_lattice_members))) for F in M]
    min_val = Fsub[F_num].LDeg
    k0 = min_val-u_vert[F_num]
    
    if sw:
        input_tuple = tuple(generic_lattice_members)
        printsw(f"Our input is of the form: {input_tuple}",sw)
        printsw("Our F's are then: ",sw)
        for i in range(len(M)):
            print(f"F_{i}{input_tuple} = {Fsub[i]}")

        printsw(f"\nWe are on step {v_num} for w_{u_num}. "+\
                f"We therefore require that v(F_{F_num}{input_tuple})"+\
                f" = {u_vert[F_num]}+k",sw)
        printsw(f"From above, the minimum possible valuation for F_{F_num}"+\
                f"{input_tuple} is {min_val}. Hence our initial k is {k0}",sw)
        
        printsw("Our system is thus:",sw)
        for i in range(len(M)):
            if i == F_num:
                print(f"v(F_{F_num}{input_tuple}) = {v_vert[i]}")
            else:
                print(f"v(F_{F_num}{input_tuple}) >= {v_vert[i]}")    
        printsw("\n"+"We now proceed with the algorithm.",sw)
    
    loop_count=0
    while loop_count>=0:
        k = k0 + loop_count
        
        if sw:
            print("....................................................................")
            print(f"Iteration {loop_count} with k = {k}")
            print("....................................................................\n")
        
        if k>=kl:
            printsw(f"k_{v_num} is greater than or equal to the k value of {kl}, "+\
                    "so we may skip this step.",sw)
            return "skip"
        
        """
        Explanation:
        ============
        ele_nums is the minimum degree needed for all relevant symbols to
        appear when substituting each DVR_symbol with a generic element
        
        any extra terms for each equation individually are removed when
        defining reduced_coeff_list
        
        """
        ele_nums=[]
        for i in range(cham_dim):
            dim_diff=[]
            for j in range(cham_dim):
                if j != F_num:
                    dim_diff+=[u_vert[j]-v_vert[i]+k-1]
                else:
                    dim_diff+=[u_vert[j]-v_vert[i]+k]
            ele_nums+=[max(dim_diff)]
        
        gen_eles = [element.generic_element(0,ele_nums[i],DVR_symbols[i])\
                    for i in range(cham_dim)]
        gen_eles_exprs = [e.expr for e in gen_eles]
        
        Feval = [F.subs(list(zip(DVR_symbols,gen_eles_exprs))) for F in Fsub]
        coeff_list = [F.coeffs_for_solver() for F in Feval]
        
        if sw:
            printsw("Choosing generic elements of appropriate degree, we have:",sw)
            for i in range(cham_dim):
                print(f"{DVR_symbols[i]} = {gen_eles[i]}")
            printsw("\nNote: If any of the polynomials are set to 0, then that means they "+\
                    "are irrelevant in determining a solution set.",sw)
            printsw("This then gives us:",sw)
            for i in range(cham_dim):
                print(f"F_{i}{tuple(DVR_symbols)} = {Feval[i]}")
            print()
        
        reduced_coeff_list = coeff_list[:]
        for i in range(cham_dim):
            if i!=F_num and len(coeff_list[i])>u_vert[i]-v_vert[i]+k:
                reduced_coeff_list[i]= \
                coeff_list[i][len(coeff_list[i])-u_vert[i]+v_vert[i]-k:len(coeff_list[i])]
        del reduced_coeff_list[F_num]
        
        system = []
        for c in reduced_coeff_list:
            system+=c

        if sw:
            start_deg = min(v_vert)
            longest = max([len(c) for c in reduced_coeff_list])
            if longest != 0:
                printsw("Our reduced system is then described as follows:",sw)
                for i in range(longest):
                    deg=i+start_deg
                    printsw(f"Degree {deg}:",sw)
                    for c in reduced_coeff_list:
                        if i < len(c):
                            print(f"{c[len(c)-i-1]} = 0")
                    print()
            else:
                printsw("Our reduced system is then empty.",sw)
            
            print()
        
        svars = []
        for i in range(len(DVR_symbols)):
            if gen_eles_exprs[i] == 0:
                svars+=[DVR_symbols[i]]
            else:
                svars+=gen_eles[i].coeffs_for_solver()
        system[:] = [val for val in system if val!=0]
        
        if system == []:
            result=[tuple(svars)]
        else:
            result=list(nonlinsolve(system,svars))

        potential_solutions = [list(zip(svars,list(result[i]))) for i in range(len(result))]
        
        Fs_from_solver = [Feval[F_num].subs(p) for p in potential_solutions]
        valid_solutions = [potential_solutions[i] for i in range(len(potential_solutions))
                           if Fs_from_solver[i].LDeg == k+u_vert[F_num]]
        
        if sw:
            printsw("The solutions to our reduced system are given by:",sw)
            for i in range(len(potential_solutions)):
                for (var,res) in potential_solutions[i]:
                    print(str(var)+" = "+str(res))
                printsw("\n"+f"Substituting into F_{F_num} for this solution gives us:",sw)
                printsw(f"F_{F_num}{tuple(generic_lattice_members)}"+\
                        f" = {Fs_from_solver[i]}",sw)
                if Fs_from_solver[i].LDeg == k+u_vert[F_num]:
                    printsw(f"Which can have valuation {k+u_vert[F_num]} and thus is valid.",sw)
                else:
                    printsw(f"Which cannot have valuation {k+u_vert[F_num]} and thus is not valid.",sw)
        
        if valid_solutions != []:
            sol_set = []
            for sol in valid_solutions:
                F = Feval[F_num].subs(sol)
                term = F.LTerm()
                free_syms = list(term.free_symbols)
                zeroes = nonlinsolve([term],free_syms)
                excl = [list(zip(free_syms,zeroes.args[i])) for i in range(len(zeroes))]
                adj_sols = [element(s**(-k+v_vert[i])*gen_eles[i].subs(sol))
                                      for i in range(len(gen_eles))]
                """
                Explanation:
                ============
                Here p_i represents a generic element that has non-negative valuation
                and N represents a generic integer.
                
                """
                p = symarray('p',len(gen_eles))
                N = Symbol('N')
                w_solution = []
                for i in range(len(gen_eles)):
                    if adj_sols[i].expr==0:
                        w_solution.append((x[i],p[i]*s**(N)))
                    else:
                        w_solution.append((x[i],adj_sols[i]+s**(adj_sols[i].UDeg+1)*p[i]))
                
                sol_set.append(system_solution(w_solution,excl))
            
            candidate = frame_vector(sol_set,k,v_num)
            if sw:
                printsw("Accounting for cases in which the lowest term of "+\
                        f"F_{F_num} is zero and adjusting according to "+\
                        "our k value, our candidate(s) are described by the following:",sw)
                printsw(candidate,sw)
            
            return candidate
        
        printsw("There are no valid solutions for this value of k, so we repeat "+\
                "this step with k = k+1.",sw)
        
        loop_count+=1
        if loop_count>=timeout:
            sys.exit(f"Timed Out after {timeout} iterations.")
            return("Timed Out")

def get_w(u_cham,v_cham,u_num,r_k_dict,timeout=oo,show_work = False):
    """
    Returns the vector of index u_num for the frame of the shared apartment.

    """
    sw = show_work
    r_vals = list(r_k_dict.keys())
    k = oo
    w = None
    for v_num in range(v_cham.dim):
        if not v_num in r_vals:
            if sw:
                print("--------------------------------------------------------------------")
                print(f"Performing step {v_num} for w_{u_num}")
                print("--------------------------------------------------------------------")
            
            candidate = get_candidate(u_cham,v_cham,u_num,v_num,kl=k,timeout=timeout,show_work=sw)
            if not candidate == "skip":
                w = candidate
                k = w.k
        else:
            if sw:
                print("--------------------------------------------------------------------")
                print(f"A previous w was determined at step {v_num}"+\
                      f", so we may skip this step for w_{u_num}")
                print("--------------------------------------------------------------------")
    return w

def get_frame(u_cham,v_cham,timeout=oo,show_work=False):
    """
    Returns a list [frame,r_k_dict]. The variable frame is the list of frame_vectors
    that will form the frame for the shared apartment. The variable r_k_dict
    is a dictionary that keeps track of the final k value obtained on step r.
    
    """
    sw=show_work
    r_k_dict = dict()
    frame = []
    for u_num in range(u_cham.dim):
        if sw:
            print("====================================================================")
            print(f"Finding w_{u_num}")
            print("====================================================================")
        result = get_w(u_cham,v_cham,u_num,r_k_dict,timeout=timeout,show_work=sw)
        frame.append(result)
        r_k_dict[result.r]=result.k
    return [frame,r_k_dict]
    
def get_w_cham(r_k_dict,show_work=False):
    """
    Returns the coordinate representation for the second chamber in the
    shared apartment as a chamber object.
    
    """
    perm = ~Permutation(list(r_k_dict.keys()))
    base = list(r_k_dict.values())
    w_cham = chamber(base,perm)
    w_cham = w_cham.change_base(1-w_cham.dim)
    return w_cham

def get_full(u_cham,v_cham,timeout=oo,show_work=False):
    result = get_frame(u_cham,v_cham,timeout=timeout,show_work=show_work)
    frame = result[0]
    r_k_dict = result[1]
    w_cham = get_w_cham(r_k_dict,show_work=show_work)
    for i in range(len(frame)):
        print(frame[i])
    print("\n"+f"The input chamber in this frame is:\n{w_cham}")

"""
========
Example:
========

"""

cham1 = chamber([0,0,0],[2,1,0])
cham2 = chamber([0,0,-1],[0,1,2])

"""
set show_work = False for just the results and not the explanation

"""
get_full(cham1,cham2,timeout=10,show_work=False)





