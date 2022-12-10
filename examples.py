# -*- coding: utf-8 -*-

import networkx as nx
from path_homology import edgelist_to_graph, H_path, H_path_A

#%%
edgelist = [# Commuting square
    'ab', 'bd', 'ac', 'cd', 
    ]
#%%
edgelist = [# Two squares
    'ab', 'bd', 'ac', 'cd',
    'be', 'ef', 'df'
    ]
#%%
edgelist = [# 2-cycle
    'ab', 'ba'
    ]
#%%
edgelist = [# 3-cycle
    'ab', 'bc', 'ca',
    ]
#%%
edgelist = [# Suspension on 3-cycle, using * <-- * --> *
    'ab', 'bc', 'ca',
    'an', 'bn', 'cn',
    'as', 'bs', 'cs',
    ]

#%%
edgelist = [# Suspension on 3-cycle, using * --> * --> *
    'ab', 'bc', 'ca',
    'an', 'bn', 'cn',
    'sa', 'sb', 'sc',
    ]
#%%
edgelist = [# Suspension on larger cycle, using * --> * --> *
    'ab', 'bc', 'cz', 'ax', 'xy', 'yz',
    'an', 'bn', 'cn', 'xn', 'yn', 'zn',
    'sa', 'sb', 'sc', 'sx', 'sy', 'sz'
    ]
#%%
edgelist = [ # prism
    'ab', 'bc', 'ca',
    'de', 'ef', 'fd',
    'ad', 'be', 'cf',
    ]
#%%
edgelist = [ # pyramid
    'ab', 'ac', 'bd', 'cd',
    'ae', 'be', 'ce', 'de'
    ]
#%%
edgelist = [ # bi-pyramid
    'sa','sb','sc', 'sd',
    'ab','bd','ac', 'cd',
    'na','nb','nc', 'nd'
    ]
#%%
edgelist = [ # octahedron
    'sa','sb','sc', 'sd', 
    'ab','bd','ac', 'cd',
    'an','bn','cn', 'dn',
    ]
#%%
edgelist = [ # octahedron without equator
    'sa', 'sb', 'sc', 'sd',
    'an', 'bn', 'cn', 'dn',
    ]
#%%
edgelist = [
    'sa', 'sb', 'sc', 'sd',
    'an', 'bn', 'cn', 'dn',
    'ab', 'bc', 'cd', 'da'
    ]
#%%
edgelist = [
    'sa', 'sb', 'sc',
    'an', 'bn', 'cn',
    'ab', 'bc', 'ac'
    ]
#%%
edgelist = [ # Boundary of subdivided 3 simplex
    'ae', 'ag', 'af',
    'be', 'bi', 'bh',
    'cf', 'cj', 'ch',
    'dg', 'di', 'dj',
    'ek', 'en',
    'fl', 'fn',
    'gk', 'gl',
    'hm', 'hn',
    'ik', 'im',
    'jl', 'jm',
    # 'ko', 'lo', 'mo', 'no' # including the 3 simplex
    ]
#%%
edgelist = [ # Cube
    'ab', 'ac', 'bd', 'cd',
    'ef', 'eg', 'gh', 'fh',
    'ae', 'bf', 'cg', 'dh'
    ]
#%%
edgelist = [ # boundary of subdivided 3 cube
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
    'v@', '@e', 'm@', '@n', 'k@', '@p', # adding the center
]
#%%
edgelist = [ # Theorem 6.1 of https://arxiv.org/pdf/1803.07497.pdf
        '12', '13', '14', '15',
        '26', '27', '36', '38', '47', '49', '58', '59',
        '60', '70', '80', '90'
    ]
#%%
edgelist = [ # Larger example
            's1', 's2', 's3', 's4', 's5',
            '16', '17', '27', '28', '38', '39', '49', '40', '50', '56',
            '6n', '7n', '8n', '9n', '0n',
    ]

#%%
G = edgelist_to_graph(edgelist)

H, C, Diffs, A, O = H_path(G, cutoff = 5)

print('C:', C)
print('H:', H)



