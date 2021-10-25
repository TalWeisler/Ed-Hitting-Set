# Ed-Hitting-Set

Let C be a collection of subsets of a finite set S. A hitting set of C is a subset of S that has a non-empty intersection
with each element of C.

Given: A collection C of subsets (each subset size d) of a set S and a parameter k.

Question: Does C have a hitting set of size k or less?

Hitting Set is NP-complete

To solve this problem we use 2 kernelizations: 
1. The sunflower lemma
2. The crown decomposition

The algorithm compare between the 2 kernalizations. 
