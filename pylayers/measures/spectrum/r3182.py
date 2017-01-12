import Gpib
import numpy as np
import matplotlib.pyplot as plt 

inst=Gpib.Gpib(0,8) 
inst.write("*IDN?")
model = inst.read(100)
print "Model : ",model
inst.write("CF 868.MZ")
inst.write("SP 50.MZ")
inst.write("TPL")

for k in range(100):
    inst.write("TBA?")
    raw = inst.read(1001*2)
    data = np.frombuffer(raw,dtype='>u2')
    try:
        t = np.vstack((t,data))
    except:
        t = data

plt.imshow(data)
plt.show()


