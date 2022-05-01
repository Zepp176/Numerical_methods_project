import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotly.figure_factory as ff

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

H = 0.01
h = 5*H/N

x = np.linspace(0, 15*H, M+1)
y = np.linspace(0,  5*H, N+1)
x, y = np.meshgrid(x, y)

u = np.array(u).reshape(M+1, N+2).T
v = np.array(v).reshape(M+2, N+1).T

velocity = np.sqrt( ( (u[1:-1,1:] + u[1:-1,:-1])/2 )**2
                  + ( (v[1:,1:-1] + v[:-1,1:-1])/2 )**2 )

plt.figure(figsize=(18.5,5))
plt.pcolormesh(x, y, velocity, cmap=cm.coolwarm, shading='auto')
plt.fill([0.03, 0.08, 0.08, 0.03], [0.02, 0.02, 0.03, 0.03], c='k')
plt.axis('equal')
plt.title('velocity')
plt.colorbar()
plt.show()
#plt.savefig("figures/velocity_result_1.png", dpi=300)

vorticity = (-u[1:,:] + u[:-1,:] - v[:,:-1] + v[:,1:])/h
vorticity[2*N//5+1:3*N//5,3*N//5+1:8*N//5] = np.zeros((N//5-1, 5*N//5-1))

U = (u[1:-1,1:] + u[1:-1,:-1])/2
V = (v[1:,1:-1] + v[:-1,1:-1])/2
X = (x[1:,1:] + x[:-1,:-1])/2
Y = (y[1:,1:] + y[:-1,:-1])/2

plt.figure(figsize=(18.5,5))
plt.pcolormesh(x, y, vorticity, cmap=cm.coolwarm, shading='auto', vmax=100, vmin=-100)
plt.colorbar()
plt.streamplot(X, Y, U, V, color='k', density=[0.5, 1], linewidth=1)
plt.fill([0.03, 0.08, 0.08, 0.03], [0.02, 0.02, 0.03, 0.03], c='k')
plt.axis('equal')
plt.title('vorticity')
plt.show()
#plt.savefig("figures/vorticity_result_1.png", dpi=300)
