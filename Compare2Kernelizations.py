import copy
import numpy as np
import timeit
import Crown_Decomposition
import Sunflower
from Hitting_Set import HittingSet


def experiment(v, e, d, k, times):
    # 0 - time ker  1 - time of solve  2 -edges after ker   3 - solution
    results_crown = np.zeros((2, times))
    edges_crown = [0]*times
    solutions_crown = [""]*times
    results_sunflower = np.zeros((2, times))
    edges_sunflower = [0] * times
    solutions_sunflower = [""] * times
    for i in range(times):
        h_crown = HittingSet(v, e, d, k)
        h_sun = copy.deepcopy(h_crown)
        # Crown decomposition
        start = timeit.default_timer()
        reduction = Crown_Decomposition.crown_decomposition_kernel(h_crown)
        results_crown[0][i] = (timeit.default_timer() - start)
        edges_crown[i] = h_crown.edges
        if reduction is not None:
            start = timeit.default_timer()
            solutions_crown[i] = h_crown.solve(reduction)
            results_crown[1][i] = (timeit.default_timer() - start)
        else:
            results_crown[1][i] = 0
            solutions_crown[i] = None
        # Sunflower lemma
        start = timeit.default_timer()
        reduction = Sunflower.sunflowerAlgorithm(h_sun)
        results_sunflower[0][i] = (timeit.default_timer() - start)
        edges_sunflower[i] = h_sun.edges
        if reduction is not None:
            start = timeit.default_timer()
            solutions_sunflower[i] = h_sun.solve(reduction)
            results_sunflower[1][i] = (timeit.default_timer() - start)
        else:
            results_sunflower[1][i] = 0
            solutions_sunflower[i] = None
    print("Results:\n [time ker,time of solve,edges after ker,solution]*times")
    print("Crown Decomposition:\n", results_crown)
    print(edges_crown)
    print(solutions_crown)
    print("Sunflower Lemma:\n", results_sunflower)
    print(edges_sunflower)
    print(solutions_sunflower)

# run the experiment:
# Input :  S   C   d  k  times
experiment(20, 10, 3, 4, 1)
