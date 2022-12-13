# -*- coding: utf-8 -*-

import networkx as nx
from path_homology import edgelist_to_graph, H_path, H_path_R

commuting_square = [
    'ab', 'bd', 'ac', 'cd', 
    ]

octahedron = [ # Example 1.27
    '01','21','03','23',
    'a0','a1','a2','a3',
    'b0','b1','b2','b3'
    ]

cube = [
    'ab', 'ac', 'bd', 'cd',
    'ef', 'eg', 'gh', 'fh',
    'ae', 'bf', 'cg', 'dh'
    ]

hollow_2x2x2_cube = [ # Example 1.30
    'ab', 'bc', 'de', 'ef', 'gh', 'hi',
    'ad', 'dg', 'be', 'eh', 'cf', 'fi',
    'jk', 'kl', 'op', 'pq', 
    'jm', 'mo', 'ln', 'nq',
    'rs', 'st', 'uv', 'vw', 'xy', 'yz',
    'ru', 'ux', 'sv', 'vy', 'tw', 'wz',
    'rj', 'ja', 'um', 'md', 'xo', 'og',
    'tl', 'lc', 'wn', 'nf', 'zq', 'qi',
    'sk', 'kb',
    'yp', 'ph',
]

exotic_sphere = [ # Theorem 6.1 of https://arxiv.org/pdf/1803.07497.pdf
        '12', '13', '14', '15',
        '26', '27', '36', '38', '47', '49', '58', '59',
        '60', '70', '80', '90'
    ]

G = edgelist_to_graph(commuting_square)

H, C, Diffs, R, O = H_path(G, cutoff = 5)

print('C:', C)
print('H:', H)



