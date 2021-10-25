import copy
import numpy as np
import random


class HittingSet:

    def __init__(self, v, e, d, k):
        self.vertices = v
        self.edges = e
        self.matrix = np.zeros((v, e), np.int8)
        self.v_degree = np.array([0] * v)
        self.e_degree = np.array([0] * e)
        self.d = d
        self.k = k
        self.values()

    def edges_with_v(self, v):
        # the set of all edges containing x
        s = set()
        for e in range(self.edges):
            if self.matrix[v][e] == 1:
                s.add(e)
        return s

    def vertices_of_e(self, edges):
        # for E ⊂ C, V(E) denotes the set of all vertices that are elements of elements of E
        s = set()
        for e in edges:
            for v in range(self.vertices):
                if self.matrix[v][e] == 1:
                    s.add(v)
        return s

    def values(self):

        for e in range(self.edges):
            # size = random.randint(1, self.d)
            size = self.d
            self.e_degree[e] = size
            arr = random.sample(range(0, self.vertices), size)
            print("E", e, " :", arr)
            for v in arr:
                self.matrix[v][e] = 1
                self.v_degree[v] += 1

    def print_hs(self):
        print("------------Hitting set Input------------")
        print(self.vertices, "vertices and", self.edges, "edges")
        print("vertices degree: ", self.v_degree)
        print("edges degree: ", self.e_degree)
        print("incidence matrix [vertices x edges]: \n", self.matrix)
        print("looking for ", self.d, "-Hitting set of size ≤", self.k)
        print("             Good-Luck!              ")
        print("-----------------------------------------")

    def print_edges(self):
        for e in range(self.edges):
            print("E", e, " : [ ", end="")
            for v in range(self.vertices):
                if self.matrix[v][e] == 1:
                    print(v, end=" ")
            print("]")

    def delete_vertex(self, v):
        # print("delete_vertex ", v)
        if self.v_degree[v] != 0:
            edges = self.edges_with_v(v)
            self.v_degree[v] = 0
            for e in edges:
                self.matrix[v][e] = 0
                self.e_degree[e] -= 1

    def delete_edge(self, e):
        # print("delete_edge ", e)
        vertices = self.vertices_of_e([e])
        self.matrix = np.delete(self.matrix, e, 1)
        self.edges -= 1
        self.e_degree = np.delete(self.e_degree, [e])
        for v in vertices:
            self.v_degree[v] -= 1

    def add_edge(self, e):
        # print("add_edge ", e)
        self.matrix = np.c_[self.matrix, e]
        self.e_degree = np.append(self.e_degree, [0])
        self.edges += 1
        for i in range(self.vertices):
            if e[i] == 1:
                self.v_degree[i] += 1
                self.e_degree[self.edges - 1] += 1

    def check_v(self, v, e_arr, left_e, k, sol):
        count = 0
        new_arr = copy.deepcopy(e_arr)
        for e in range(self.edges):
            if (self.matrix[v][e] == 1) and (e_arr[e] == 1):
                new_arr[e] = 0
                count += 1
        left_e -= count
        sol.append(v)
        if left_e == 0:
            return sol
        if (k-1) == 0:
            sol.remove(v)
            return None
        for new_v in range(v-1, -1, -1):
            if self.v_degree[new_v] > 0:
                res = self.check_v(new_v, new_arr, left_e, k-1, sol)
                if res is not None:
                    return res
        left_e += count
        sol.remove(v)
        return None

    def solve(self, sol):
        if self.edges == 0 and self.k >= 0:
            return sol
        if self.k <= 0:
            return None
        for v in range(self.vertices-1, -1, -1):
            if self.v_degree[v] > 0:
                res = self.check_v(v, [1]*self.edges, self.edges, self.k, sol)
                if res is not None:
                    return res
        return None

    def reset_edge(self, e):
        vertices = self.vertices_of_e([e])
        for v in vertices:
            self.matrix[v, e] = 0
        self.e_degree[e] = 0

    def delete_edge_2(self, e):
        vertices = self.vertices_of_e([e])
        for v in vertices:
            self.matrix[v, e] = 0
            self.v_degree[v] -= 1
        self.edges -= 1
        self.e_degree[e] = 0

    def is_dup(self, e, index):
        for i in range(index, -1, -1):
            if sum(np.bitwise_xor(e, self.matrix[:, i])) == 0:
                self.reset_edge(i)