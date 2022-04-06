import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

f = open("data/file.txt", "r")
data = f.read().split("\n")[:-1]
res = int(data[0])
data = data[1:]
for i in range(len(data)):
    data[i] = float(data[i])

N = res*5
M = res*15

"""for i in range(3*res, 8*res):
    for j in range(2*res, 3*res):
        data[j+i*N] = -2.0;"""

H = 1
h = H/res

x = np.linspace(h/2, 15*H - h/2, M)
y = np.linspace(h/2,  5*H - h/2, N)
x, y = np.meshgrid(x, y)

field = np.empty((N, M))
for i in range(M):
    for j in range(N):
        field[j,i] = data[j+i*N]

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(x, y, field, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()
