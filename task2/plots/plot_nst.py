import matplotlib

import numpy as np
import matplotlib.pyplot as plt

size = 101

x = np.linspace(0, 10, size)

y = np.loadtxt("z_V_nst.txt", dtype=float)

# Создаем сетку X, Y для pcolormesh
#X, Y = np.meshgrid(x, y)  # Обрезаем x и y на 1 элемент, чтобы они совпали с размерностью Z

# Отображаем график
plt.plot(x, y)
plt.savefig('z_V_nst_e.jpg')

plt.show()
