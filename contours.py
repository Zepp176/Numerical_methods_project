import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

case = 4
t_stars = [1,2,5,10,20,30,40,50]

dt_star = 0.062

for t_star in t_stars:
    step = int(t_star/dt_star+1)
    
    filename = "D:/Numerical_Methods_data/case{}/step_{}.txt".format(case, step)
    f = open(filename, "r")
    data = f.read().split("\n")
    f.close()
    
    res = data[0].split(" ")
    M = int(res[0])
    N = int(res[1])
    t_star = float(res[2])
    x_mesh = float(res[5])
    y_mesh = float(res[6])
    u = data[1].split(" ")[:-1]
    v = data[2].split(" ")[:-1]
    
    print(t_star)
    
    for j in range(len(u)):
        u[j] = float(u[j])
    for j in range(len(v)):
        v[j] = float(v[j])
    
    x = np.linspace(0, 15, M+1) + x_mesh
    y = np.linspace(0,  5, N+1) + y_mesh
    x, y = np.meshgrid(x, y)
    
    u = np.array(u).reshape(M+1, N+2).T
    v = np.array(v).reshape(M+2, N+1).T
    
    w_star = (-u[1:,:] + u[:-1,:] - v[:,:-1] + v[:,1:])*N/5.0
    w_star[w_star > 20] = 20
    w_star[w_star < -20] = -20
    
    plt.figure(figsize=(12, 4))
    levels = np.linspace(-20, 20, 30)
    plt.contourf(x, y, w_star, levels=levels, cmap=cm.seismic)
    plt.colorbar()
    plt.contour(x, y, w_star, levels=levels, colors='k', linewidths=0.2, linestyles='solid')
    plt.fill(np.array([3, 8, 8, 3]) + x_mesh, np.array([2, 2, 3, 3]) + y_mesh, c='y', zorder=10)
    plt.axis('equal')
    plt.gca().set_adjustable("box")
    plt.gca().set_xlim(left=0, right=15)
    plt.gca().set_ylim(bottom=0, top=5)
    plt.xlabel("$x/H_{box}$")
    plt.ylabel("$y/H_{box}$")
    plt.title("dimensionless vorticity field $\omega^*$ at $t^* = {:.1f}$".format(t_star))
    plt.tight_layout()
    plt.savefig("figures/vorticity/case{}/vorticity_contour_case{}_time{:.0f}.png".format(case, case, t_star), dpi=300)