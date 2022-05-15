import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

case = 1
t_stars = [12.5]

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
    u_mesh = float(res[3])
    v_mesh = float(res[4])
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
        
    # compute stagnation point
    j = int(2.5*(N/5))
    i = int(8.5*(N/5))
    min_idx = np.argmin(np.abs((u[j, i:i+5*int(N/5)] + u[j+1, i:i+5*int(N/5)])/2)) + i
    length = (min_idx-8*50)/50
    print('length of recirculation zone: {}'.format(length))
    
    u = (u[1:-1, 1:] + u[1:-1, :-1]) / 2
    v = (v[1:, 1:-1] + v[:-1, 1:-1]) / 2
    x = (x[1:, 1:] + x[1:, :-1] + x[:-1, 1:] + x[:-1, :-1]) / 4
    y = (y[1:, 1:] + y[1:, :-1] + y[:-1, 1:] + y[:-1, :-1]) / 4
    
    plt.figure(figsize=(12, 4))
    plt.streamplot(x, y, u - u_mesh, v - v_mesh, density=2, color='k', linewidth=0.2, arrowsize=0)
    plt.pcolormesh(x, y, np.sqrt(u**2 + v**2), shading="gouraud", cmap=cm.Blues, vmax=2, vmin=0)
    plt.colorbar()
    plt.plot([min_idx/int(N/5)], [2.5], "xr")
    plt.fill(np.array([3, 8, 8, 3]) + x_mesh, np.array([2, 2, 3, 3]) + y_mesh, c='y', zorder=10)
    plt.axis('equal')
    plt.gca().set_adjustable("box")
    plt.gca().set_xlim(left=0, right=15)
    plt.gca().set_ylim(bottom=0, top=5)
    plt.xlabel("$x/H_{box}$")
    plt.ylabel("$y/H_{box}$")
    plt.title("$\\frac{||u(x,y)||}{U_{\infty}}$ and the streamlines at " + "$t^* = {:.1f}$ (length of recirculation zone: {})".format(t_star, length))
    plt.tight_layout()
    plt.savefig("figures/streamlines/case{}/streamlines_case{}_time{:.0f}.png".format(case, case, t_star), dpi=300)