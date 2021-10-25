import copy
import math

import numpy as np
import timeit

import Change
import Crown_Decomposition
import Sunflower
from Hitting_Set import HittingSet


def experiment(v, e, d, k, times):
    for i in range(times):
        h_crown = HittingSet(v, e, d, k)
        h_crown2 = copy.deepcopy(h_crown)
        # Crown decomposition
        # reduction = Crown_Decomposition.crown_decomposition_kernel(h_crown)
        reduction2 = Change.crown_decomposition_kernel(h_crown2)
        if reduction2 is not None:
            print(h_crown.solve(reduction2))
        else:
            print("There is no solution")


def check(s, d):
    res = 0
    for x in range(1, d-1, 1):
        num1 = math.factorial(s)
        num2 = math.factorial(x) * math.factorial(s-x)
        num = num1 / num2
        print(x, num)
        res += num
    print("check =", res)


check(30, 15)
experiment(20, 100, 8, 2, 1)