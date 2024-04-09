import matplotlib.pyplot as plt
import numpy as np

def calculate_T_vectors(P, t, beta, c):
    # Вектор напрямку T+ and T-
    m = len(P)
    T_plus = np.zeros_like(P)
    T_minus = np.zeros_like(P)

    for i in range(1, m - 1):
        T_plus[i] = (((1 - t) * (1 - beta) * (1 - c)) / 2) * (P[i + 1] - P[i]) + (
                    ((1 - t) * (1 + beta) * (1 + c)) / 2) * (P[i] - P[i - 1])
        T_minus[i] = (((1 - t) * (1 - beta) * (1 + c)) / 2) * (P[i + 1] - P[i]) + (
                    ((1 - t) * (1 + beta) * (1 - c)) / 2) * (P[i] - P[i - 1])

    return T_plus, T_minus


def build_TSV_spline(P, T_plus, T_minus, num_points=100):
    m = len(P)
    #Створює рівномірно розподілений вектор від 0 до 1 num_points раз
    t = np.linspace(0, 1, num_points)
    #Двовимірний масив нулів
    S = np.zeros((num_points, 2))

    for i in range(m - 1):
        for j in range(num_points):
            # S - розраховується за допомогою формули кубічного сплайну
            S[j] = ((1 - t[j]) ** 3) * P[i] + (3 * (1 - t[j]) ** 2 * t[j]) * (P[i] + T_plus[i]) + (
                        3 * (1 - t[j]) * t[j] ** 2) * (P[i + 1] - T_minus[i + 1]) + (t[j] ** 3) * P[i + 1]
        plt.plot(S[:, 0], S[:, 1], 'b-')  # Візуалізація ліній
        plt.plot(P[:, 0], P[:, 1], 'ro')  # Візуалізуємо контрольні точки

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('TSV-Spline')
    plt.grid(True)
    plt.show()


# Приклад використання
P = np.array([[0, 0.4], [1.5, 0.5], [3, -0.3], [6, 0.7], [9, 0.3]])
t = 0.5
beta = 0.5
c = 0.5

T_plus, T_minus = calculate_T_vectors(P, t, beta, c)
build_TSV_spline(P, T_plus, T_minus)


