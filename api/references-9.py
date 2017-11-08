import matplotlib.pyplot as plt
df = DF(b=np.array([1,1],a=np.array([1,-1]))
N = 100 
s = np.zeros(N)
s[0] = 1
y = df.filter(s)
plt.stem(np.arange(N),x)
plt.stem(np.arange(N),y)
plt.show()
