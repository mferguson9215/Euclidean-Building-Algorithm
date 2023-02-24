# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 11:36:42 2023

@author: micha
"""

from sympy import *
from sympy.combinatorics import Permutation
init_printing(perm_cyclic=False, pretty_print=False)

def l_add(L1,L2):
    if len(L1)>=len(L2):
        return [L1[i]+L2[i] for i in range(len(L2))]+L1[len(L2):len(L1)]
    else:
        return [L1[i]+L2[i] for i in range(len(L1))]+L2[len(L1):len(L2)]

def l_subtract(L1,L2):
    if len(L1)>=len(L2):
        return [L1[i]-L2[i] for i in range(len(L2))]+L1[len(L2):len(L1)]
    else:
        return [L1[i]-L2[i] for i in range(len(L1))]+[-val for val in L2[len(L1):len(L2)]]

class element:
    """
    Represents an element in the underlying field.
    
    """
    def __init__(self,expr,s=Symbol('s')):
        """
        Should be passed a sympy expression 'expr' which can be written as a linear combination
        over the powers of s.
        
        self.coeff_dict is a dictionary using the powers of s as keys for the
        corresponding coefficients in the base field.
        
        self.degs is a list of the degrees of each key
        
        self.LDeg and self.UDeg are the minimum and maximum degrees of the terms
        
        """
        self.s=s
        if isinstance(expr,int) or isinstance(expr,float):
            self._expr=collect(expand(sympify(expr)),self.s)
        else:
            self._expr=collect(expand(expr),self.s)
        self.coeff_dict=self._expr.as_coefficients_dict(self.s)
        self.degs = [term.as_coeff_exponent(self.s)[1] for term in self.coeff_dict.keys()]
        if self.expr != 0:
            self.LDeg = min(self.degs)
            self.UDeg = max(self.degs)
        else:
            self.LDeg = -oo
            self.UDeg = -oo
    
    @property
    def expr(self):
        return self._expr
    
    @expr.setter
    def expr(self,new):
        """
        Equivalent to calling self.__init__(new)
        
        """
        
        if isinstance(expr,int) or isinstance(expr,float):
            self._expr=collect(expand(sympify(expr)),self.s)
        else:
            self._expr=collect(expand(expr),self.s)
        self.coeff_dict=self._expr.as_coefficients_dict(self.s)
        self.degs = [term.as_coeff_exponent(self.s)[1] for term in self.coeff_dict.keys()]
        if self.expr != 0:
            self.LDeg = min(self.degs)
            self.UDeg = max(self.degs)
        else:
            self.LDeg = -oo
            self.UDeg = -oo
    
    """
    Arithmetical operations:
    
    """
    def __add__(self,other):
        if isinstance(other,element):
            return element(self.expr+other.expr)
        else:
            return element(self.expr+sympify(other))

    def __radd__(self,other):
        if isinstance(other,element):
            return element(other.expr+self.expr)
        else:
            return element(sympify(other)+self.expr)
        
    def __sub__(self,other):
        if isinstance(other,element):
            return element(self.expr-other.expr)
        else:
            return element(self.expr-sympify(other))

    def __rsub__(self,other):
        if isinstance(other,element):
            return element(other.expr-self.expr)
        else:
            return element(sympify(other)-self.expr)
    
    def __neg__(self):
        return element(-self.expr)
    
    def __mul__(self,other):
        if isinstance(other,element):
            return element(self.expr*other.expr)
        else:
            return element(self.expr*sympify(other))

    def __rmul__(self,other):
        if isinstance(other,element):
            return element(other.expr*self.expr)
        else:
            return element(sympify(other)*self.expr)
    
    
    def _sympy_(self):
        return self.expr
    
    def __str__(self):
        return str(self.expr)
    
    def __repr__(self):
        return str(self)
    
    def LTerm(self):
        """
        Returns the coefficient of the lowest degree term.
        
        """
        return self.coeff_dict[self.s**self.LDeg]
    
    def UTerm(self):
        """
        Returns the coefficient of the highest degree term.
        
        """
        return self.coeff_dict[self.s**self.UDeg]
    
    def nTerm(self,n):
        """
        Returns the coefficient of the nth degree term.
        
        """
        return self.coeff_dict[self.s**n]
    
    def all_terms(self):
        """
        Returns output identical to that of the all_terms function for polynomials
        from sympy, but allows negative powers.
        
        Examples:
        ========
        1)
            >>> element(x**3 + 2*x - 1, Symbol('x')).all_terms()
            [((3,),1),((2,),0),((1,),2),((0,),-1)]
        
            >>> Poly(x**3 + 2*x - 1, x).all_terms()
            [((3,),1),((2,),0),((1,),2),((0,),-1)]
        
        2)
            >>> element(x**3 + 2*x - 1 + 1/x**2, Symbol('x')).all_terms()
            [((3,),1),((2,),0),((1,),2),((0,),-1),((-1,),0),((-2,),1)]
        
            >>> Poly(x**3 + 2*x - 1 + 1/x**2, x).all_terms()
            PolynomialError: 1/x contains an element of the set of generators.
        
        """
        if self.expr == 0:
            return []
        else:
            return [((self.UDeg-i,),self.nTerm(self.UDeg-i)) for i in range(self.UDeg-self.LDeg+1)]
    
    def coeffs_for_solver(self):
        """
        A more convenient version of the above when using nonlinsolve from sympy.
        
        """
        if self.expr == 0:
            return []
        else:
            return [self.nTerm(self.UDeg-i) for i in range(self.UDeg-self.LDeg+1)]
    
    def subs(self,*args,**kwargs):
        e = element(self.expr.subs(*args,**kwargs))
        return e
    
    @staticmethod
    def generic_element(min_deg,max_deg,coeff_sym=Symbol('q'),s=Symbol('s')):
        """
        Returns an element with symbols as the coefficients for each term.
        
        Example:
        ========
        >>> element.generic_element(-2,3)
        (q_-2)/s**2+(q_-1)/s+q_0+q_1*s+q_2*s**2+q_3*s**3
        
        """
        expr = 0
        if min_deg <= max_deg:
            q = dict([(i,Symbol(f'{coeff_sym}_{i}')) for i in range(min_deg,max_deg+1)])
            for i in range(min_deg,max_deg+1):
                expr = expr + q[i]*s**i
        return element(expr,s)

class chamber:
    """
    This class represents the chamber inside a generic apartment with respect to
    the given base vertex and permutation in the link of that vertex.
    
    Example:
    ========
    >>> perm = Permutation([0,2,1])
    >>> base = [2,0,-1]
    >>> cham = chamber(perm,base)
    
    >>> cham.dim
    3
    
    >>> cham.fund_chamber_coords
    [[0,0,0],[-1,0,0],[-1,0,-1]]
    
    >>> cham.vertices
    [[2,0,-1],[1,0,-1],[1,0,-2]]

    """
    
    def __init__(self,base,perm):
        self._base = base
        if not isinstance(perm, Permutation):
            self._perm = Permutation(perm)
        else:
            self._perm = perm
        self._dim = len(base)
        self._fund_chamber_coords = [[-1]*i+[0]*(self._dim-i) \
                                     for i in range(self._dim)]
        self._vertices = [l_add(self._base,self._perm(self._fund_chamber_coords[i])) \
                          for i in range(self._dim)]
    
    @property
    def base(self):
        return self._base
    
    @property
    def fund_chamber_coords(self):
        return self._fund_chamber_coords
    
    @property
    def perm(self):
        return self._perm
    
    @property
    def vertices(self):
        return self._vertices
    
    @property
    def dim(self):
        return self._dim
    
    def __str__(self):
        return f"Base Point: {self.base}\n"+ \
               f"Permutation (as list): {self.perm.list()}\n"+ \
               f"Coordinates: {self.vertices}"
    
    def __repr__(self):
        return str(self)
    
    @staticmethod
    def lattice_containment(lat1,lat2):
        """
        Returns an integer representing whether or not lat1 is contained in lat2
        if viewed as lattices with the given coordinates in the same apartment.
        
        Notes:
        ======
        1) lat1 and lat2 are lists of integers. A list of integers [a_0,..,a_n]
           represents the lattice s**(a_0)Ru_0+...+s**(a_n)Ru_n, where R is
           the DVR and u_0,...,u_n is an unspecified basis.
        
        2) An integer is returned instead of a boolean so that the number of lattices
           a lattice is contained in can be counted. See the function 'arrange_vertices'.
         
        """
        return int(all(lat1[i] <= lat2[i] for i in range(len(lat1))))
    
    @staticmethod
    def arrange_vertices(verts):
        """
        Returns verts sorted by lattice containment. verts should be a list of
        lists of integers of the same dimension
        
        Example:
        ========
        >>> chamber.arrange_vertices([[-4,0,2],[-3,0,2],[-4,0,1]])
        [[-3,0,2],[-4,0,2],[-4,0,1]]
        
        """
        num_contained_in = [sum(map(lambda v2 : chamber.lattice_containment(v1,v2),verts))-1 \
                                for v1 in verts]
        return Permutation(num_contained_in)(verts)
        
    @classmethod
    def from_verts(cls,vertex_list):
        """
        from_verts returns the chamber with vertices vertex_list.
        
        Example:
        ========
        >>> cham = chamber.from_verts([[-4,0,2],[-3,0,2],[-4,0,1]])
        
        >>> cham.vertices
        [[-3,0,2],[-4,0,2],[-4,0,1]]
        
        >>> cham.base
        [-3,0,2]
        
        >>> cham.perm
        [0,2,1]

        """
        vertices = chamber.arrange_vertices(vertex_list)
        base_point = vertices[0]
        shifted = [l_subtract(v,base_point) for v in vertices]
        p = [l_subtract(shifted[i],shifted[i-1]).index(-1) for i in range(1,len(shifted))]
        Perm=~Permutation(p+[next(iter(set(range(len(vertex_list)))-set(p)))])
        return chamber(base_point,Perm)
    
    def change_base(self,index):
        """
        Returns the chamber object whose base, permuation, and coordinates
        correspond to change of the base vertex to the lattice starting at
        index 'index' in the ascending chain of lattices.
        
        Example:
        ========
        >>> cham = chamber.from_verts([[-1,1,0],[-2,1,0],[-2,0,0]])
        
        >>> print(cham)
        Base Point: [-1, 1, 0]
        Permutation (as list): [0, 1, 2]
        Coordinates: [[-1, 1, 0], [-2, 1, 0], [-2, 0, 0]]
        
        >>> print(cham.change_base(-2))
        Base Point: [-1, 2, 1]
        Permutation (as list): [2, 0, 1]
        Coordinates: [[-1, 2, 1], [-1, 1, 1], [-1, 1, 0]]
        
        """
        mod = index % self.dim
        mag = int((index-mod)/self.dim)
        scaled_lats = [l_subtract(v,[mag]*self.dim) for v in self.vertices]
        new_verts = scaled_lats[mod:self.dim]+[l_subtract(L,[1]*self.dim) \
                                                 for L in scaled_lats[0:mod]]
        return chamber.from_verts(new_verts)








