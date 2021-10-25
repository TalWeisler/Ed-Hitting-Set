import numpy as np
import math
import copy


def sunflowerAlgorithm(h):

    for i in range(h.edges - 1, -1, -1):
        if not (h.e_degree[i] == 0):
            h.is_dup(h.matrix[:, i], i - 1)

    for i in range(h.edges - 1, -1, -1):
        if h.e_degree[i] == 0 :  # empty edge or edge.degree < h.d
            h.delete_edge(i)

    a_kernel = copy.deepcopy(h)
    for i in range(a_kernel.edges - 1, -1, -1) :
        if not (a_kernel.e_degree[i]== h.d):   #empty edge or edge.degree < h.d
            a_kernel.delete_edge_2(i)


    if a_kernel.edges < math.factorial(a_kernel.d) * math.pow(a_kernel.k , a_kernel.d):
        return []

    else :
        (core, g) = sunflowerAlgorithmRec(a_kernel, np.array([0] * h.vertices), [])

        if sum(core) == 0:
            return None

        else:
            for i in range(len(g) - 1, -1, -1):
                h.delete_edge(g[i])
            h.add_edge(core)
            return sunflowerAlgorithm(h)

        return []


def sunflowerAlgorithmRec(a_kernel, core, g):

    # build G
    j = 0
    while sum(a_kernel.matrix[:, j]) == 0:
        j = j + 1
    g = [j]  # first group S1
    g_vec = a_kernel.matrix[:, j]

    for i in range(1, a_kernel.matrix.shape[1]):
        if not (a_kernel.e_degree[i] == 0):  # col is not deleted
            temp = g_vec
            temp = np.bitwise_and(temp, a_kernel.matrix[:, i])
            if sum(temp) == 0:  # the group is safe to G
                g_vec = g_vec + a_kernel.matrix[:, i]
                g.append(i)
            else:
                continue


    if a_kernel.k <= len(g):  # found sunflower with k-petals
        return core, g

    else:
        (v_index, v_edges) = findMaxV(a_kernel)
        for i in range(a_kernel.edges - 1, -1, -1):
            if not v_edges.__contains__(i):
                a_kernel.delete_edge_2(i)
        core[v_index] = 1
        a_kernel.delete_vertex(v_index)
        a_kernel.d = a_kernel.d - 1  # update d
        return sunflowerAlgorithmRec(a_kernel, core, g)


def findMaxV(a_kernel):
    max = 0
    v = 0
    for i in range(len(a_kernel.v_degree)):
        if a_kernel.v_degree[i] > max:
            max = a_kernel.v_degree[i]
            v = i
    return v, a_kernel.edges_with_v(v)






