import copy
import numpy as np
from itertools import combinations


"""
The following two pre-processing rules, dubbed domination rules, 
are used by previously known dHS kernelization algorithms. 
– Vertex Domination Rule: If x,y ∈ S are such that E(x) ⊂ E(y),then delete x.
– Edge Domination Rule: If e1,e2 ∈ C are such that V({e1}) ⊂ V({e2}) then delete e2
"""


def subset_vertex(h, v1, v2):
    for e in range(h.edges):
        if (h.matrix[v1][e] == 1) and (h.matrix[v2][e] == 0):
            return False
    return True


def check_vertex(h, v1):
    for v2 in range(h.vertices - 1, -1, -1):
        if (h.v_degree[v2] != 0) and (v2 != v1):
            if subset_vertex(h, v1, v2):
                return True
    return False


def vertex_domination_rule(h):
    # If x,y ∈ S are such that E(x) ⊂ E(y),then delete x.
    # print("*vertex_domination_rule*")
    change = 0
    for v in range(h.vertices - 1, -1, -1):
        if h.v_degree[v] != 0:
            if check_vertex(h, v):
                h.delete_vertex(v)
                change = 1
    return change


def subset_edge(h, e1, e2):
    if h.e_degree[e1] == h.e_degree[e2]:
        return False
    for v in range(h.vertices):
        if (h.matrix[v][e1] == 1) and (h.matrix[v][e2] == 0):
            return False
    return True


def is_dup(h, e1, e2):
    edge_e1 = h.matrix[:, e1]
    edge_e2 = h.matrix[:, e2]
    if sum(np.bitwise_xor(edge_e1,edge_e2)) == 0:
        return True
    return False


def check_edge(h, e1):
    for e2 in range(h.edges - 1, -1, -1):
        if e2 is not e1:
            if is_dup(h, e1, e2):
                return True
            if subset_edge(h, e2, e1):
                return True
    return False


def edge_domination_rule(h):
    # If e1,e2 ∈ C are such that V({e1}) ⊂ V({e2}) then delete e2
    # print("*edge_domination_rule*")
    change = 0
    for e in range(h.edges-1, -1, -1):
        if check_edge(h, e):
            h.delete_edge(e)
            change = 1
    return change


def get_vertex(h, e):
    for v in range(h.vertices):
        if h.matrix[v][e] == 1:
            return v
    return -1


def singleton(h):
    # print("*singleton*")
    s = []
    for e in range(h.edges-1, -1, -1):
        if h.e_degree[e] == 1:
            v = get_vertex(h, e)
            if v != -1:
                if v not in s:
                    s.append(v)
                    h.k -= 1
                h.delete_edge(e)
    return s


"""
High-Degree rule for d-HS: 
    If a subedge e ⊂ S satisfies:
    (i) |e| = d − 2 , and
    (ii) e is the pair-wise intersection of more than k edges.
    Then add e to C and delete all elements of C that contain e (as a subset).

If the high-degree rule does not apply to a dHS instance (S,C,k), then we shall say that (S,C,k) 
is reduced. This could be the result of applying the rule iteratively until it cannot be applied.
"""


def convert_edge(h, e):
    converted_e = np.array([0] * h.vertices)
    for v in e:
        converted_e[v] = 1
    return converted_e


def check_group(h, edges, e_checked):
    for e in edges:
        e_group = h.matrix[:, e]
        if sum(np.bitwise_and(e_group, e_checked)) != (h.d - 2):
            return False
    return True


def intersection(h, option):
    edges = []
    e_option = convert_edge(h, option)
    for e in range(h.edges):
        e_checked = h.matrix[:, e]
        if sum(np.bitwise_and(e_option, e_checked)) == (h.d - 2):
            if check_group(h, edges, e_checked):
                edges.insert(0,e)
    if len(edges) > h.k:
        for e in edges:
            h.delete_edge(e)
        h.add_edge(e_option)
        return 1
    return 0


def high_degree_rule(h):
    # print("*high_degree_rule*")
    change = 0
    combination = list(combinations(range(h.vertices), (h.d - 2)))
    for e in combination:
        check = intersection(h, e)
        change = check | change
    return change


"""
Weakly related: 
    Two edges are weakly related if their intersection contains at most d −2 elements 
    and neither of them is a subset of the other
"""


def is_not_subset_edges(e1, e2):
    num1 = 1
    num2 = 1
    different = np.bitwise_xor(e1,e2)
    for v in range(len(e1)):
        if different[v] == 1:
            if e1[v] == 1:
                num1 = 0
            elif e2[v] == 1:
                num2 = 0
            if num1 == 0 and num2 == 0:
                return True
    return False


def weakly_related(h, e1, e2):
    edge_e1 = h.matrix[:, e1]
    edge_e2 = h.matrix[:, e2]
    similar = np.bitwise_and(edge_e1, edge_e2)
    if sum(similar) <= (h.d - 2):
        if is_not_subset_edges(edge_e1, edge_e2):
            return True
    return False


"""
maximal set W  of weakly related edges - 
    The construction starts by placing all edges of size d−2 or less in W. 
    At each subsequent stage, an element e of C is selected and checked against the current elements of W.
    Then e is added to W if the elements of W ∪ {e} are weakly related.
"""


def check_union(h, w, e):
    for e2 in range(h.edges):
        if w[e2] == 1:
            if not weakly_related(h, e, e2):
                return False
    return True


def maximal_set(h):
    # print("*maximal_set*")
    w = np.array([0] * h.edges)
    check = np.array([0] * h.edges)
    for e in range(h.edges):
        if h.e_degree[e] <= (h.d - 2):
            w[e] = 1
        else:
            check[e] = 1
    for e in range(h.edges):
        if check[e] == 1:
            if check_union(h, w, e):
                w[e] = 1
    return w


"""
Let A be a solution of a reduced instance (S,C,k) of dHS, and let W be maximal set. 
Denote by We the set of edges of W that contain a subedge e properly. 
If e is a (d−2)-subedge then, by the high-degree rule,|We|<= k.
To achieve our sought upper bound, we apply the following reduction algorithm.

The High Occurrence Rule:
    For i=d−2 downto 1 do
        For each i-subedge e of W do 
            if |We| >k^(d−1−i),then
                Add e to W  <-אמור כבר להיות שם לא???
                Delete (from C and W) all edges containing e

Lemma 5. Let (S, C,k) be a reduced yes instance of dHS and let W be a maximal set of weakly related edges that result 
from applying the high-occurrence rule. Then |W |<= k^(d−1)
"""


def w_e(h, w, e1):
    ret = copy.deepcopy(w)
    for v in e1:
        edges_with_v = h.matrix[v]
        ret = np.bitwise_and(ret, edges_with_v)
    return ret


def high_occurrence_rule(h):
    w = maximal_set(h)
    # print("*high_occurrence_rule*")
    delete_from_w = [0]*h.edges
    for i in range(h.d - 2, 0, -1):
        comb = list(combinations(range(h.vertices), i))
        for e in comb:
            we = w_e(h, w, e)
            if sum(we) > pow(h.k, (h.d - 1 - i)):
                for d_e in range(len(we)):
                    if we[d_e] == 1:
                        delete_from_w[d_e] = 1
                        w[d_e] = 0
                h.add_edge(convert_edge(h, e))
                w = np.append(w, [1])
                delete_from_w = np.append(delete_from_w, [0])
    for x in range(h.edges - 1, -1, -1):
        if delete_from_w[x] == 1:
            h.delete_edge(x)
            w = np.delete(w, [x])
    return w


"""
Crown reduction
Let I be the complement of V(W) in S, and let H be the set of all (d − 1)-subedges of W . 
Then we observe the following:
    • Elements of C whose cardinality is less than d are placed in W .
    • I is an independent set.
    • Every edge e of size d that contains an element x of I consists of x and a (d − 1)-subedge H. 
        Otherwise, e could be added to W , violating W ’s maximality.
    • |V(H)|<=>dk^(d−1): each edge of size d of W gives rise to d (d − 1)-subedges.
As in the 3HS case, we proceed by constructing the bipartite graph GI,H . 
The rest of the work is identical to the crown construction and reduction that was used in the 3HS case
"""


def v_in_w(h, w, v):
    edges_with_v = h.matrix[v]
    edges_with_v = np.bitwise_and(edges_with_v, w)
    if sum(edges_with_v) > 0:
        return True
    return False


def find_neighbors(h, i_arr, e1):
    neighbors = [0] * h.vertices
    edge_e1 = h.matrix[:, e1]
    for e2 in range(h.edges):
        if h.e_degree[e2] == h.d:
            edge_e2 = h.matrix[:, e2]
            if sum(np.bitwise_and(edge_e1, edge_e2)) == (h.d - 1):
                diff = np.bitwise_xor(edge_e1, edge_e2)
                np.bitwise_or(neighbors, diff)
    neighbors = np.bitwise_and(neighbors, i_arr)
    return neighbors


def get_neighbor(s, e, v_num):
    for v in range(v_num):
        if s[e][v] == 1:
            return v
    return -1


def get_neighbors(s, group, v_num):
    neighbors = [0] * v_num
    for e in group:
        neighbors = np.bitwise_or(neighbors, s[e])
    return neighbors


def array2num(arr):
    res = []
    for x in range(arr):
        if arr[x] == 1:
            res.append(x)
    return res


def stable(s, group, neighbors):
    new_z = []
    new_nz = []
    for e in group:
        for v in neighbors:
            if s[e][v] == 1:
                new_z.append(e)
                new_nz.append(v)
                neighbors.remove(v)
    return [new_z, new_nz]


def hall(h, i_arr, h_arr, s):
    # print("*hall*")
    if len(h_arr) == 1:
        v = get_neighbor(s, h_arr[0], h.vertices)
        if v != -1:
            return [v]
        return []
    for size in range(1, len(h_arr), 1):
        comb = list(combinations(h_arr, size))
        for group in comb:
            neighbors = get_neighbors(s, group, h.vertices)
            if sum(neighbors) > 0:
                neighbors = array2num(neighbors)
                if len(neighbors) == size:
                    for z in group:
                        h_arr.remove(z)
                    for n_z in neighbors:
                        i_arr.remove(n_z)
                    return neighbors.extend(hall(h, i_arr, h_arr, s))
                if len(neighbors) > size:
                    s = stable(s, np.asarray(group), neighbors)
                    for z in s[0]:
                        h_arr.remove(z)
                    for n_z in s[1]:
                        i_arr.remove(n_z)
                    res = s[1]
                    res.extend(hall(h, i_arr, h_arr, s))
                    return res
    return []


def different_arrays(h, e1, e2):
    if e1 == e2:
        return -1
    if not (h.e_degree[e2] == h.d):
        return -1
    ret = -1
    edge_e1 = h.matrix[:, e1]
    edge_e2 = h.matrix[:, e2]
    for v in range(h.vertices):
        n = edge_e1[v] - edge_e2[v]
        if n == -1:
            return -1
        if n == 1:
            ret = v
    return ret


def construct_crown(h, w):
    # print("*construct_crown*")
    # based on Expansion lemma (q=1) and Hall's theorem
    i_arr = [1]*h.vertices
    h_arr = []
    s = []

    # complement of V(W) in S
    for v in range(h.vertices):
        if v_in_w(h, w, v):
            i_arr[v] = 0

    # (d − 1)-subedges of W
    for e1 in range(h.edges):
        if (w[e1] == 1) and (h.e_degree[e1] == h.d - 1):
            h_arr.append(e1)
            s.append(find_neighbors(h, i_arr, e1))
        else:
            arr = [0]*h.vertices
            s.append(arr)

    if len(i_arr) > len(h_arr):
        return hall(h, i_arr, h_arr, s)
    else:
        return None


"""
Algorithm for d-Hitting-Set:
Input: Instance (S, C, k) of DHS.
Output: Either NO if (S, C) has no hitting set of size ≤ k, or an instance (S', C', k')
        and a partial solution A such that k' = k − |A| and |S'| ≤ 5k^2 + k.
Begin
    Vertex Domination Rule
    Edge Domination Rule
    While C contains a singleton edge {x} add x to A, delete {x} from C and decrement k.
    Apply the high-degree rule
    Construct W - maximal set
    high occurrences rule 
    If |W| > k^2
        return NO
    Construct the simple bipartite graph G_(I,H)
    If |I| > |H|
        Call Construct Crown on G_(I,H) to obtain a crown (H', I', M)
        Remove all elements of I'
    Return the (possibly) new instance
End
"""


def repeat_steps(h, n, reduction):
    change = 1
    if n == 1:
        vertex_domination_rule(h)
    if n == 2:
        while change == 1:
            change = edge_domination_rule(h)
            if change == 1:
                change = vertex_domination_rule(h)
    if n == 3:
        while change == 1:
            x = len(reduction)
            reduction.extend(singleton(h))
            if h.k == 0 or h.edges == 0:
                return
            if len(reduction) > x:
                change = edge_domination_rule(h)
                if change == 1:
                    vertex_domination_rule(h)
                else:
                    change = vertex_domination_rule(h)
            else:
                change = 0
    if n == 4:
        while change == 1:
            if h.edges == 0 and h.k >= 0:
                return reduction
            if h.k <= 0:
                return None
            change = high_degree_rule(h)
            if change == 1 and h.k > 0 and h.edges > 0:
                x = len(reduction)
                reduction.extend(singleton(h))
                if h.k == 0 or h.edges == 0:
                    return
                if len(reduction) > x:
                    change = edge_domination_rule(h)
                    if change == 1:
                        vertex_domination_rule(h)
                    else:
                        change = vertex_domination_rule(h)
                else:
                    change = 0


def crown_decomposition_kernel(h):
    reduction = []
    repeat_steps(h, 1, [])  # vertex_domination_rule
    repeat_steps(h, 2, [])  # edge_domination_rule + vertex_domination_rule
    repeat_steps(h, 3, reduction)  # singleton + edge_domination_rule + vertex_domination_rule
    if h.edges == 0 and h.k >= 0:
        return reduction
    if h.k <= 0:
        return None
    repeat_steps(h, 4, reduction)  # high_degree_rule + singleton + edge_domination_rule + vertex_domination_rule
    if h.edges == 0 and h.k >= 0:
        return reduction
    if h.k <= 0:
        return None
    w = high_occurrence_rule(h)  # include maximal set
    if h.edges == 0 and h.k >= 0:
        return reduction
    if h.k <= 0:
        return None
    if sum(w) > pow(h.k, 2):
        return None
    y = construct_crown(h, w)
    if y is not None:
        for v in y:
            h.delete_vertex(v)
    return reduction

