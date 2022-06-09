import numpy as np
import h5py
import matplotlib.pyplot as plt

f = h5py.File('C:\\Users\\DELL\\Desktop\\vector8.h5', 'r')
data = f["explicit solution"]
data1 = data[:]

x1 = np.linspace(0, 1, 5001)
y1 = np.sin(np.pi*x1)/np.pi**2

print(max(abs(y1-data1)))
# print(max(abs(y-data2))/max(abs(y-data1)))

plt.title('HPC,Grid=5000', size=20)
plt.xlabel('x', size=20)
plt.ylabel('T', size=20)
plt.plot(x1, data1, '--', label='Explicit')
plt.plot(x1, y1, '--', label='Analytic')
plt.legend(loc='upper right')
plt.show()
