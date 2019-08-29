#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 11:58:55 2018

@author: hielke
"""

from itertools import starmap, compress, chain
from operator import __or__


from pprint import PrettyPrinter as pp
pprint = pp().pprint

def gene_proximity_matrix(vals, gene_distance):
    length = len(vals)
    matrix = [[None for _ in range(length)] for _ in range(length)]
    for j, j_v in enumerate(vals):
        for i, i_v in enumerate(vals):
            coords = sorted(chain(j_v, i_v))
            res = abs(coords[1] - coords[2]) <= gene_distance
            matrix[i][j] = res
    for i in range(length):
        matrix[i][i] = False
    return matrix    
    


def proximity_matrix(vals, n=None):
    """
    Returns a 2D proximity matrix with the distance between two points.
    If n is set it returns this matrix as a boolean matrix where true denotes
    that the distance between the two points is equal or less than n.
    --
    vals: list
    n: int
    """
    length = len(vals)
    matrix = [[None for _ in range(length)] for _ in range(length)]
    for j, j_v in enumerate(vals):
        for i, i_v in enumerate(vals):
            res = abs(j_v - i_v)
            if n:
                 res = res <= n
            matrix[i][j] = res
    for i in range(length):
        matrix[i][i] = False
    return matrix

def create_clusters_from_proximity_matrix(D, clusters=None, skip=None):
    """
    D: a boolean proximity matrix 
    """
    if not skip:
        skip = set()
    if not clusters:
        clusters = list(map(lambda x : [x], range(len(D[0]))))
    
    for cur_id, vals in enumerate(D):
        if cur_id in skip:
            continue
        for other_id, v in enumerate(vals):
            if other_id in skip:
                continue
            if v:
                update_D(cur_id, other_id, D, clusters, skip)
                create_clusters_from_proximity_matrix(D, clusters, skip)
                break
            
    return list(compress(clusters, clusters))

def create_clusters_from_values(vals, n):
    """
    Returns the created cluster as indices of the input
    """
    D = proximity_matrix(vals, n=n)
    return create_clusters_from_proximity_matrix(D)


    

def update_D(cur_id, other_id, D, clusters, skip):
    """
    Removes the other_id and puts it in skip.
    Removes other id from the clusters
    Updates the proximity matrix with the newly formed cluster
    """
    D[cur_id] = list(starmap(__or__, zip(D[cur_id], D[other_id])))
    D[cur_id][cur_id] = False
    skip.add(other_id)
    clusters[cur_id].append(other_id)
    print(cur_id, other_id)
    print(clusters)
    clusters[other_id] = None
    
if __name__ == '__main__':
    A = [28, 3, 6, 8, 16, 18, 22, 30]
#    pprint(proximity_matrix(A))
#    pprint(proximity_matrix(A, n=5))
    pprint(create_clusters_from_values(A, n=5))