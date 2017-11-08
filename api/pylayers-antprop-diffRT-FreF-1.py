import matplotlib.pyplot as plt
import numpy as np
x = np.logspace(-4,2,400);
F = FreF(x)
plt.semilogx(x,,np.abs(F))
plt.grid()
