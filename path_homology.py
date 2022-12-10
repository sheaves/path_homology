# -*- coding: utf-8 -*-
"""
Updated on Oct 6, 2022.

@author: LiangZe

A script to compute the path homology for a digraph.

Assumptions made:
    0) No self-loops or multiple edges (they will be ignored)
    1) Degeneracies are killed
    2) k = Real Numbers
"""

from sympy import Matrix
import numpy as np
import pandas as pd
import networkx as nx

def A_n(G, cutoff = 5):
    """
    Generate paths in G, up to length = cutoff.
    G : networkx DiGraph
    
    Paths are represented as tuples of vertices
    e.g. ('a', 'b', 'c', 'a')
    
    Returns a list of lists:
    [
     list of vertices, 
     list of edges, 
     list of length 2 paths, 
     ..., 
     list of length cutoff paths
    ]
    """
    # Initialize list of paths
    A = []

    # A[0]
    A.append(list(G))
    
    # A[1]
    A.append(list(G.edges))
    
    # A[n], 1 < n <= cutoff
    for n in range(2,cutoff + 1):
        A_n = []
        for e in A[1]:
            for f in A[n-1]:
               if e[1] == f[0]:
                   A_n.append((e[0],) + f)
        A.append(A_n)
    
    return A

def d(path):
    """
    Compute the differential of a path, assuming all paths are allowed.
    """
    coeffs = {}
    for i in range(len(path)):
        head = path[:i]
        tail = path[i+1:]
        
        # Don't include degeneracies
        if (len(head) == 0) or (len(tail) == 0) or head[-1] != tail[0]: 
            coeffs[head + tail] = (-1)**i
    
    return pd.Series(coeffs, name = path)

def D_full(A, n):
    """
    Compute differentials of all paths of length n
    
    A : list of list of paths (see docstring of A_n)
    """
    if len(A[n]) > 0:
        df = pd.concat([d(path) for path in A[n]], axis = 1)    
        return df.fillna(0)
    else:
        return pd.DataFrame()
    
def O_n(A):
    """
    Generate Omega chain complex for path homology 
    from list of list of paths, A (see docstring of A_n)

    O[n] contains generators for Omega_n,
    expressed as lists of linear combinations of paths in A[n]
    
    Linear combinations are expressed as dictionaries:
        - keys are elements of A[n]
        - values are the coefficients in the linear combination
    """
    O = []
    
    # O_0
    O.append([{(v,):1} for v in A[0]])
    
    # O_1
    O.append([{e:1} for e in A[1]])
    
    # O_n
    for n in range(2,len(A)):
        disallowed = D_full(A,n).drop(A[n-1], errors = 'ignore')
            
        if len(disallowed) > 0:
            O_n = []
            basis = Matrix(disallowed).nullspace()
            for v in basis:
                v = np.array(v, dtype = float)[:,0]
                v = pd.Series(v, index = A[n])
                v = v[v != 0]
                O_n.append(dict(v))
        else:
            O_n = [{p:1} for p in A[n]]            
        O.append(O_n)    
        
    return O
    
def pprint(coeffs_dict, sep = ''):
    """
    Pretty printing of linear combinations of paths.

    coeffs_dict stores linear combinations represented as a dictionary (see O_n docstring).
    
    Paths are represented as a string of vertices, separated by sep.    
    """
    result = ''
    for key in coeffs_dict:
        key_string = sep.join(key)
        coeff = coeffs_dict[key]
        sign = '+' if coeff > 0 else '-'
        val  = ' ' if abs(coeff) == 1 else f'{abs(coeff)} '
        result += f'{sign}{val}{key_string} '
    return result

def D(A, O, n, sep = ''):
    """
    Generate matrix for n^th differential in Omega chain complex,
    going from O[n] --> O[n-1]
    
    Matrix is expressed as a DataFrame, where:
        - columns = generators of O[n]
        - rows = generators of O[n-1]
    
    sep is used for pretty printing of generators
    """
    D_f = D_full(A,n)
    
    if n == 0:
        D_f.columns = [f"+ {c} " for c in D_f.columns]
        return D_f        
    
    if len(A[n-1]) > 0:
        index = pd.MultiIndex.from_tuples(A[n-1])
    else:
        index = []
    
    # Differentials of generators of O[n] in terms of A[n-1]
    D = pd.DataFrame(index = index)
    for coeffs in O[n]:
        v = pd.Series(0, index = index)
        for key in coeffs:
            v = v.add(D_f[key]*coeffs[key], fill_value = 0)
        D[pprint(coeffs, sep)] = v
    
    # Generators of O[n-1] in terms of A[n-1]
    B = pd.DataFrame(index = index)
    for coeffs in O[n-1]:
        v = pd.Series(0, index = index)
        v = v.add(pd.Series(coeffs), fill_value = 0)
        B[pprint(coeffs, sep)] = v    
        
    # Return differentials of generators of O[n] in terms of generators of O[n-1]
    Binv = Matrix(B).pinv()
    Binv = np.array(Binv, dtype = float)
    Binv = pd.DataFrame(Binv, index = B.columns, columns = B.index)
    
    return Binv @ D

def H_path_A(A, sep = ''):
    """
    Compute the path homology of a path set A.
    
    Returns:
        H : dimensions of path homology
        C : dimensions of Omega chain complex
        Diffs : differentials
        A : allowed paths
        O : generators of Omega
    """    
    O = O_n(A)
    cutoff = len(A) - 1

    Diffs = [D(A,O,i,sep) for i in range(cutoff + 1)]
    C = [df.shape[1] for df in Diffs]
    
    H = []
    
    # H_0
    ker = len(A[0])
    img = len(Matrix(Diffs[1]).columnspace())
    H.append(ker - img)
    
    # H_n, n > 0
    for n in range(1, cutoff):
        ker = len(Matrix(Diffs[n]).nullspace())
        img = len(Matrix(Diffs[n+1]).columnspace())
        H.append(ker - img)
  
    return H, C, Diffs, A, O    
    

def H_path(G, cutoff = 5, sep = ''):
    """
    Compute the path homology of G, up to cutoff.
    
    Returns:
        H : dimensions of path homology
        C : dimensions of Omega chain complex
        Diffs : differentials
        A : allowed paths
        O : generators of Omega
    """
    A = A_n(G, cutoff = cutoff)
    
    return H_path_A(A, sep = sep)


def edgelist_to_graph(edgelist, sep = '', vertices = None):
    """
    Generates a networkx DiGraph from a list of directed edges.
    Self-loops and parallel edges will be ignored.
    
    Examples of edgelists:
        - ['ab', 'ac', 'bd', 'bc'] (sep = '', i.e. no separator)
        - ['10:11', '10:12', '11:13', '12:13']  (sep = ":")
    
    sep is required if some vertices have length > 1
    
    Vertices will be inferred from the edgelist, but you may specify them:
        - have more control over the order of vertices
        - have "floating" vertices that don't belong to any edge
    
    Examples of vertices:
        - ['a', 'b', 'c', 'd']
        - ['10', '11', '12', '13']
        
    """
    G = nx.DiGraph()
    
    # Add vertices, if specified
    if vertices is not None:
        G.add_nodes_from(vertices)
    
    # Add edges
    for e in edgelist:
        if sep == '':
            e0, e1 = e
        else:
            e0, e1 = e.split(sep)
        
        # Check for self-loops. Only add non-self-loops
        if e0 != e1:
            G.add_edge(e0, e1)
        else:
            G.add_node(e0)
            
    return G




#%%
"""
IGNORE. Not updated; may not be compatible with other code.
# Path homology for arbitrary path complexes (not arising from graphs)

def H_path_A(A):
    O = O_n(A)
    H = []
    
    # H_0
    ker = len(A[0])
    img = len(Matrix(D(A,O,1)).columnspace())
    H.append(ker - img)
    
    # H_n
    for n in range(1,len(A)-1):
        ker = len(Matrix(D(A,O,n)).nullspace())
        img = len(Matrix(D(A,O,n+1)).columnspace())
        H.append(ker - img)
        
    return H

# Example
A = [
   ['a', 'b', 'c', 'd'],
   ['ab', ], #, 'bd', 'bc', 'cb'],
   ['bcd', 'acd', 'abd', 'abc'],
   #['acb', 'cbd'],
   #['acbd'],
   [],
   [],
   [],
   []
   ]

H = H_path_A(A)

"""