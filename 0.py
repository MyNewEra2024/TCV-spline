import numpy as np
import matplotlib.pyplot as plt


def compute_TSV_spline(P, t, beta, c):
    m = len(P) - 1
    T_plus = np.zeros_like(P)
    T_minus = np.zeros_like(P)

    for i in range(1, m):
        # Compute T[i]+
        T_plus[i] = (((1 - t) * (1 - beta) * (1 - c)) / 2) * (P[i + 1] - P[i]) + (
                    ((1 - t) * (1 + beta) * (1 + c)) / 2) * (P[i] - P[i - 1])

        # Compute T[i]-
        T_minus[i] = (((1 - t) * (1 - beta) * (1 + c)) / 2) * (P[i + 1] - P[i]) + (
                    ((1 - t) * (1 + beta) * (1 - c)) / 2) * (P[i] - P[i - 1])

    return T_plus, T_minus


def plot_TSV_spline(P, T_plus, T_minus):
    m = len(P) - 1
    for i in range(1, m):
        # Plotting the Hermite curve segment
        t = np.linspace(0, 1, 100)

        # Hermite basis functions
        H0 = 2 * t ** 3 - 3 * t ** 2 + 1
        H1 = -2 * t ** 3 + 3 * t ** 2
        H2 = t ** 3 - 2 * t ** 2 + t
        H3 = t ** 3 - t ** 2

        # Compute X using Hermite basis functions
        X = np.outer(H0, P[i]) + np.outer(H1, P[i + 1]) + np.outer(H2, T_minus[i]) + np.outer(H3, T_plus[i])

        plt.plot(X[:, 0], X[:, 1], 'b-')

    plt.scatter(P[:, 0], P[:, 1], c='r', zorder=5)
    plt.quiver(P[1:-1, 0], P[1:-1, 1], T_plus[1:-1, 0], T_plus[1:-1, 1], angles='xy', scale_units='xy', scale=1,
               color='g', zorder=10)
    plt.quiver(P[1:-1, 0], P[1:-1, 1], T_minus[1:-1, 0], T_minus[1:-1, 1], angles='xy', scale_units='xy', scale=1,
               color='y', zorder=10)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('TSV Spline')
    plt.grid(True)
    plt.show()


# Example usage:
P = np.array([
    [1, 2],
    [3, 1],
    [4, 3],
    [6, 2],
    [10,6]
])

t = 0.5
beta = 0.5
c = 0.5

T_plus, T_minus = compute_TSV_spline(P, t, beta, c)
plot_TSV_spline(P, T_plus, T_minus)

