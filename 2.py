import numpy as np
import matplotlib.pyplot as plt

# Визначення точок Pi
points = np.array([[0, 0], [1, 2], [3, 1], [4, 3], [5, 0]])
m = len(points) - 1  # Кількість точок

# Визначення параметрів t, c, β
t = 0  # Натяг
c = 0  # Неперервність
β = 1  # Загальний нахил

# Обчислення векторів швидкостей T[i]+ та T[i]-
T_plus = []
T_minus = []

for i in range(1, m):
    P_prev = points[i - 1]
    P_curr = points[i]
    P_next = points[i + 1]

    T_plus.append(((1 - t) * (1 - β) * (1 - c) / 2) * (P_next - P_curr) +
                  ((1 - t) * (1 + β) * (1 + c) / 2) * (P_curr - P_prev))

    T_minus.append(((1 - t) * (1 - β) * (1 + c) / 2) * (P_next - P_curr) +
                   ((1 - t) * (1 + β) * (1 - c) / 2) * (P_curr - P_prev))

# Побудова ермітових кривих між кожною парою точок
def hermite_curve(P0, P1, T0, T1, num_points=100):
    u = np.linspace(0, 1, num_points)
    H00 = 2*u**3 - 3*u**2 + 1
    H10 = u**3 - 2*u**2 + u
    H01 = -2*u**3 + 3*u**2
    H11 = u**3 - u**2

    curve = np.outer(H00, P0) + np.outer(H10, T0) + \
            np.outer(H01, P1) + np.outer(H11, T1)
    return curve

# З'єднання кривих для створення ТСВ-сплайна
spline = np.empty((0, 2))
for i in range(m):
    P0 = points[i]
    P1 = points[i + 1]
    T0 = T_minus[i - 1] if i > 0 else np.zeros(2)
    T1 = T_plus[i] if i < m - 1 else np.zeros(2)
    spline = np.vstack((spline, hermite_curve(P0, P1, T0, T1)))

# Візуалізація ТСВ-сплайна
plt.plot(spline[:, 0], spline[:, 1], 'b', label='ТСВ-сплайн')
plt.plot(points[:, 0], points[:, 1], 'ro', label='Ключові точки')
plt.legend()
plt.show()
