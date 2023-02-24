# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 12:28:09 2023

@author: micha
"""
from sympy import Symbol

"""
Custom package:
    
"""
from chamber_package import element


class system_solution:
    
    def __init__(self,equalities,exclusions):
        """
        equalities should be a list of tuples and exclusions should be a
        list of lists of tuples
        
        Notes:
        ======
        1) An element (sym,sol) in equalities represents the relation sym = sol
        
        2) A pair (sym,excl) in a list in exclusions represents the relation
           sym != sol
        
        """
        self.equalities = equalities
        self.exclusions = exclusions
    
    def __str__(self):
        """
        Returns a string representation of self
        
        Example:
        ========
        
        >>> x = symarray('x',3)
        >>> y, z, s = symbols('y z s')
        
        >>> equalities = [(x[0],-y/(2*s**2)+(-y/4-z/2)/s),\
                          (x[1],-y/(2*s)),(x[2],y/s**2+z/s)]
        >>> exclusions = [[y != 0]]
        
        >>> print(system_solution(equalities,exclusions))
        Solution:
        x_0 = -y/(2*s**2)+(-y/4-z/2)/s)
        x_1 = -y/(2*s)
        x_2 = y/s**2+z/s

        Exclusions:
        Set 1:
        y != 0
        
        """
        eq_str = "Solution:\n"
        for (var,sol) in self.equalities:
            eq_str += f"{var} = {sol}\n"
        
        ex_str = "\n"+"Exclusions:\n"
        for i in range(len(self.exclusions)):
            ex_str+=f"Set {i+1}:\n"
            for (var,exc) in self.exclusions[i]:
                ex_str+=f"{var} != {exc}\n"
            if i < len(self.exclusions)-1:
                ex_str+="\n"
        
        return eq_str+ex_str
    
    def __repr__(self):
        return str(self)

class frame_vector:
    """
    This class represnts a vector or candidate vector for the frame of the shared apartment.
    solutions should be a list of system_solution instances.
    
    """
    def __init__(self,solutions,k,r,sym=Symbol('w')):
        self.solutions = solutions
        self.k = k
        self.r = r
        self.sym = sym
    
    def __str__(self):
        out = f"The vector {self.sym}_{self.r} is determined as follows:\n\n"
        for i in range(len(self.solutions)):
            out+=f"System Solution {i+1}:\n\n" + f"{self.solutions[i]}"
        return out
    
    def __repr__(self):
        return str(self)

   
