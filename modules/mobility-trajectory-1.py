t = np.arange(0,10,0.01)
x = 2*t*np.cos(t)
y = 3*t*np.sin(t) 
pt =np.vstack((x,y)).T
traj = Trajectory(t,pt)
traj.plot()
plt.show()
