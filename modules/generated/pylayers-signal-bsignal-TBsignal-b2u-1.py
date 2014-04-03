from pylayers.signal.bsignal import *
import matplotlib.pyplot as plt
x = np.array( [ 1, 3 , 6 , 11 , 18])
y = np.array( [ 0,1 ,-5, 8 , 10])
sb = TBsignal(x,y)
su20 = sb.b2u(20)
su100 = sb.b2u(100)
fi = plt.figure()
sb.stem()
fig,ax = su20.plot(color='k')
fig,ax = su100.plot(color='r')
ti = plt.title('b2u : sb(blue) su20(black) su200(red)')
