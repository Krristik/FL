#!/usr/bin/env python3
import sys

def objective(X: list):
    x1, x2, x3 = X
    return (x3 + 2) * x1**2 * x2

def g1_func(X: list):
    return 1 - X[1]**3 * X[2] / (71785 * X[0]**4)

def g2_func(X: list):
    return 1 / (5108 * X[0]**2) + (4 * X[1]**2 - X[0] * X[1]) / (12566 * (X[1] * X[0]**3 - X[0]**4)) - 1

def g3_func(X: list):
    return 1 - 140.55 * X[0] / (X[1]** 2 * X[2])

def g4_func(X: list):
    return (X[0] + X[1]) / 1.5 - 1

def penalty_func(X: list):
    alpha = 1e6
    obj_value = objective(X)
    g1 = g1_func(X)
    g2 = g2_func(X)
    g3 = g3_func(X)
    g4 = g4_func(X)

    penalty = (
        max(0, g1)**2 +
        max(0, g2)**2 +
        max(0, g3)**2 +
        max(0, g4)**2
    )

    return obj_value + alpha * penalty



if __name__ == "__main__":
    X = list()
    for i in sys.argv[1:]:
        X.append(float(i))
    print(penalty_func(X))
