import numpy as np
import matplotlib.pyplot as plt

case = 4
steps = np.arange(1, 808, 2)

dt_star = 0.062
mesh_Re = np.empty(len(steps))
mesh_Re_v = np.empty(len(steps))

for i, step in enumerate(steps):
    
    filename = "D:/Numerical_Methods_data/case{}/step_{}.txt".format(case, step)
    f = open(filename, "r")
    data = f.read().split("\n")
    f.close()
    
    res = data[0].split(" ")
    M = int(res[0])
    N = int(res[1])
    
    u = data[1].split(" ")[:-1]
    v = data[2].split(" ")[:-1]
    
    for j in range(len(u)):
        u[j] = float(u[j])
    for j in range(len(v)):
        v[j] = float(v[j])
    
    u = np.array(u).reshape(M+1, N+2).T
    v = np.array(v).reshape(M+2, N+1).T
    
    mesh_Re[i] = np.max((np.abs((u[1:-1,1:] + u[1:-1,:-1])/2) + np.abs((v[1:,1:-1] + v[:-1,1:-1])/2)) * 500 / 50) 
    
    mesh_Re_v[i] = np.max(np.abs(-u[1:,:] + u[:-1,:] - v[:,:-1] + v[:,1:])*500/50)
    
    print(step)

#%%

plt.figure(figsize=(5, 4))
plt.plot(np.linspace(0, 50, len(mesh_Re)), mesh_Re)
plt.plot(np.linspace(0, 50, len(mesh_Re_v)), mesh_Re_v)
plt.title("Mesh Reynold's numbers for case {}".format(case))
plt.xlabel("$t^*$")
plt.grid()
plt.ylim(top=60, bottom=5)
plt.legend(["$Re_h$", "$Re_{h,\\omega}$"])
plt.tight_layout()
plt.savefig("figures/evolution_mesh_Re_case{}.png".format(case), dpi=300)


