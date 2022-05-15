import numpy as np
import matplotlib.pyplot as plt

def get_LD(u, v, P, M, N):
    C_L = 0.0
    C_D = 0.0
    res = int(N//5.0)
    Re = 500
    U_inf = 0.05
    
    # top face
    j = 3*res
    for i in range(3*res, 8*res):
        f = u[j+2, i] + u[j+2, i+1] - u[j, i] - u[j, i+1]
        C_D += 0.5 * f / Re
        C_L -= 2 * P[j, i] / (U_inf**2 * 50)
    
    # bottom face
    j = 2*res-1
    for i in range(3*res, 8*res):
        f = u[j, i] + u[j, i+1] - u[j+2, i] - u[j+2, i+1]
        C_D += 0.5 * f / Re
        C_L += 2 * P[j, i] / (U_inf**2 * 50)
    
    # left face
    i = 3*res-1
    for j in range(2*res, 3*res):
        f = v[j, i] + v[j+1, i] - v[j, i+2] - v[j+1, i+2]
        C_L += 0.5 * f / Re
        C_D += 2 * P[j, i] / (U_inf**2 * 50)
    
    # right face
    i = 8*res
    for j in range(2*res, 3*res):
        f = v[j, i+2] + v[j+1, i+2] - v[j, i] - v[j+1, i]
        C_L += 0.5 * f / Re
        C_D -= 2 * P[j, i] / (U_inf**2 * 50)
    
    return C_L, C_D

case = 4
steps = np.arange(1, 808, 2)

C_Ls = np.empty(len(steps))
C_Ds = np.empty(len(steps))

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
    P = data[3].split(" ")[:-1]
    
    for j in range(len(u)):
        u[j] = float(u[j])
    for j in range(len(v)):
        v[j] = float(v[j])
    for j in range(len(P)):
        P[j] = float(P[j])
    
    u = np.array(u).reshape(M+1, N+2).T
    v = np.array(v).reshape(M+2, N+1).T
    P = np.array(P).reshape(M, N).T
    
    C_Ls[i], C_Ds[i] = get_LD(u, v, P, M, N)
    
    if i % 5 == 0:
        print(step)

#%%

plt.figure(figsize=(5, 4))
plt.plot(np.linspace(0, 50, len(steps)), C_Ls)
plt.plot(np.linspace(0, 50, len(steps)), C_Ds)
plt.legend(["$C_L$", "$C_D$"])
plt.ylim(top=10, bottom=-7)
plt.xlabel("$t^*$")
plt.title("Lift and drag coefficients for case {}".format(case))
plt.grid()
plt.tight_layout()
plt.savefig("figures/LD_case{}.png".format(case), dpi=300)