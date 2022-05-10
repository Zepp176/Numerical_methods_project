from matplotlib import pyplot as plt
from celluloid import Camera
import numpy as np
import matplotlib.cm as cm
import matplotlib.animation as animation

nb_frames = 400
filename = 'animation_test'
folder = 'video_test'

fig = plt.figure(figsize=(14,5))
camera = Camera(fig)

for i in range(nb_frames):
    print("step {}/{}".format(i, nb_frames))

    f = open("data/{}/step_{}.txt".format(folder, i+1), "r")
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

    for j in range(len(u)):
        u[j] = float(u[j])
    for j in range(len(v)):
        v[j] = float(v[j])

    x = np.linspace(0, 15, M+1) + x_mesh
    y = np.linspace(0,  5, N+1)
    x, y = np.meshgrid(x, y)

    u = np.array(u).reshape(M+1, N+2).T + u_mesh
    v = np.array(v).reshape(M+2, N+1).T

    velocity = np.sqrt( ( (u[1:-1,1:] + u[1:-1,:-1])/2 )**2
                      + ( (v[1:,1:-1] + v[:-1,1:-1])/2 )**2 )

    plt.pcolormesh(x, y, velocity, cmap=cm.coolwarm, shading='auto', vmax=0.12, vmin=0)
    plt.fill(np.array([3, 8, 8, 3]) + x_mesh, np.array([2, 2, 3, 3]), c='k')
    plt.axis('equal')
    plt.tight_layout()

    camera.snap()

anim = camera.animate()
anim.save('figures/{}.mp4'.format(filename), fps=24)
