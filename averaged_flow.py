import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

case = 2 # only 1 or 2

dt_star = 0.062

filename = "D:/Numerical_Methods_data/case{}/step_{}.txt".format(case, 1)
f = open(filename, "r")
data = f.read().split("\n")
f.close()

res = data[0].split(" ")
M = int(res[0])
N = int(res[1])

P = data[3].split(" ")[:-1]
for j in range(len(P)):
    P[j] = float(P[j])
P = np.array(P).reshape(M, N).T

x = np.linspace(0, 15, M+1)
y = np.linspace(0,  5, N+1)
x, y = np.meshgrid(x, y)
x = (x[1:, 1:] + x[1:, :-1] + x[:-1, 1:] + x[:-1, :-1]) / 4
y = (y[1:, 1:] + y[1:, :-1] + y[:-1, 1:] + y[:-1, :-1]) / 4

u = np.zeros((N+2, M+1))
v = np.zeros((N+1, M+2))
P = np.zeros((N, M))

counter = 0

print("from {} to {}".format(int(20/dt_star+1), int(50/dt_star+1)+1))
for step in np.arange(int(20/dt_star+1), int(50/dt_star+1)+1):
    
    filename = "D:/Numerical_Methods_data/case{}/step_{}.txt".format(case, step)
    f = open(filename, "r")
    data = f.read().split("\n")
    f.close()
    
    u_temp = data[1].split(" ")[:-1]
    v_temp = data[2].split(" ")[:-1]
    P_temp = data[3].split(" ")[:-1]
    
    for j in range(len(u_temp)):
        u_temp[j] = float(u_temp[j])
    for j in range(len(v_temp)):
        v_temp[j] = float(v_temp[j])
    for j in range(len(P_temp)):
        P_temp[j] = float(P_temp[j])
    
    u_temp = np.array(u_temp).reshape(M+1, N+2).T
    v_temp = np.array(v_temp).reshape(M+2, N+1).T
    P_temp = np.array(P_temp).reshape(M, N).T
    
    u += u_temp
    v += v_temp
    P += P_temp
    
    counter += 1
    
    if counter % 5 == 0:
        print(step)

u /= counter
v /= counter
P /= counter

w_star = (-u[1:,:] + u[:-1,:] - v[:,:-1] + v[:,1:])*N/5.0
w_star[w_star > 10] = 10
w_star[w_star < -10] = -10

u = (u[1:-1, 1:] + u[1:-1, :-1]) / 2
v = (v[1:, 1:-1] + v[:-1, 1:-1]) / 2

# compute stagnation point
j = int(2.5*(N/5))
i = int(8.5*(N/5))
min_idx = np.argmin(np.abs((u[j, i:i+5*int(N/5)] + u[j+1, i:i+5*int(N/5)])/2)) + i
length = (min_idx-8*50)/50
print('length of recirculation zone: {}'.format(length))
 
plt.figure(figsize=(12, 4))
plt.streamplot(x, y, u, v, density=2, color='k', linewidth=0.2, arrowsize=0)
plt.pcolormesh(x, y, np.sqrt(u**2 + v**2), shading="gouraud", cmap=cm.Blues, vmax=2, vmin=0)
plt.colorbar()
plt.plot([min_idx/int(N/5)], [2.5], "xr")
plt.fill(np.array([3, 8, 8, 3]), np.array([2, 2, 3, 3]), c='y', zorder=2)
plt.axis('equal')
plt.gca().set_adjustable("box")
plt.gca().set_xlim(left=0, right=15)
plt.gca().set_ylim(bottom=0, top=5)
plt.xlabel("$x/H_{box}$")
plt.ylabel("$y/H_{box}$")
plt.title("$\\frac{||\\bar u(x,y)||}{U_{\infty}}$ and the streamlines of the time-averaged flow " + "(length of recirculation zone: {})".format(length))
plt.tight_layout()
plt.savefig("figures/time_averaged_flow/streamlines_avg.png", dpi=300)

x = np.linspace(0, 15, M+1)
y = np.linspace(0,  5, N+1)
x, y = np.meshgrid(x, y)

plt.figure(figsize=(12, 4))
levels = np.linspace(-10, 10, 30)
plt.contourf(x, y, w_star, levels=levels, cmap=cm.seismic)
plt.colorbar()
plt.contour(x, y, w_star, levels=levels, colors='k', linewidths=0.2, linestyles='solid')
plt.fill(np.array([3, 8, 8, 3]), np.array([2, 2, 3, 3]), c='y', zorder=10)
plt.axis('equal')
plt.gca().set_adjustable("box")
plt.gca().set_xlim(left=0, right=15)
plt.gca().set_ylim(bottom=0, top=5)
plt.xlabel("$x/H_{box}$")
plt.ylabel("$y/H_{box}$")
plt.title("dimensionless vorticity field $\omega^*$ of time-averaged flow")
plt.tight_layout()
plt.savefig("figures/time_averaged_flow/vorticity_avg.png", dpi=300)
