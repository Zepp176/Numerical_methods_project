from matplotlib import pyplot as plt
from celluloid import Camera
import numpy as np
import matplotlib.cm as cm
import matplotlib.animation as animation

nStep = 120

fig = plt.figure(figsize=(14,5))
camera = Camera(fig)

for i in range(nStep):
    print("step {}/{}".format(i, nStep))
    f = open("data/video/step_{}.txt".format(i+1), "r")
    data = f.read().split("\n")
    f.close()
    res = data[0].split(" ")
    M = int(res[0])
    N = int(res[1])
    u = data[1].split(" ")[:-1]
    v = data[2].split(" ")[:-1]

    for i in range(len(u)):
        u[i] = float(u[i])
    for i in range(len(v)):
        v[i] = float(v[i])

    H = 0.01
    h = 5*H/N

    x = np.linspace(0, 15*H, M+1)
    y = np.linspace(0,  5*H, N+1)
    x, y = np.meshgrid(x, y)

    u = np.array(u).reshape(M+1, N+2).T
    v = np.array(v).reshape(M+2, N+1).T

    velocity = np.sqrt( ( (u[1:-1,1:] + u[1:-1,:-1])/2 )**2
                      + ( (v[1:,1:-1] + v[:-1,1:-1])/2 )**2 )

    plt.pcolormesh(x, y, velocity, cmap=cm.coolwarm, shading='auto', vmax=0.12, vmin=0)
    plt.fill([0.03, 0.08, 0.08, 0.03], [0.02, 0.02, 0.03, 0.03], c='k')
    plt.axis('equal')
    plt.tight_layout()

    camera.snap()

anim = camera.animate()
anim.save('figures/animation.mp4', fps=6)
