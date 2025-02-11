import matplotlib

import numpy as np
import matplotlib.pyplot as plt

x_size = 100
y_size = 1991
z_size = x_size * y_size

x = np.linspace(0, 10, x_size+1)
y = np.linspace(0, 194.1, y_size+1)

z_H = np.loadtxt("z_H.txt", dtype=float)

# Переворачиваем z в нужную размерность
Z_H = z_H.reshape(y_size - 1, x_size - 1)  # Размер Z должен быть (1940, 100)

# Создаем сетку X, Y для pcolormesh
X, Y = np.meshgrid(x[:-1], y[:-1])  # Обрезаем x и y на 1 элемент, чтобы они совпали с размерностью Z

# Отображаем график
plt.pcolormesh(X, Y, Z_H, shading='flat')
plt.colorbar()
plt.savefig('z_H_.jpg')

plt.show()


