import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

f = open("data/file.txt", "r")
data = f.read().split("\n")
res = data[0].split(" ")
M = int(res[0])
N = int(res[1])
u = data[1].split(" ")[:-1]
v = data[2].split(" ")[:-1]
P = data[3].split(" ")[:-1]

for i in range(len(u)):
    u[i] = float(u[i])
for i in range(len(v)):
    v[i] = float(v[i])
for i in range(len(P)):
    P[i] = float(P[i])

H = 1
h = 5*H/N

x = np.linspace(0, 15*H, M+1)
y = np.linspace(-h/2,  5*H + h/2, N+2)
x, y = np.meshgrid(x, y)

field = np.empty((N+2, M+1))
for i in range(M+1):
    for j in range(N+2):
        field[j, i] = u[j+i*(N+2)]

print(len(x), len(x[0]))
print(len(y), len(y[0]))
print(len(field), len(field[0]))

print(u)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(x, y, field, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()
