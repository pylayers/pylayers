from pylayers.simul.link import *
import pdb


# In[2]:

DL=DLink(L=Layout('Luebbers.ini'))


# In[3]:

# get_ipython().magic(u'matplotlib inline')
# DL.L.showG('i')


# In[7]:

DL.a = np.array(([25,21.67,2.]))
# DL.a = np.array(([37.5,52.2,2.]))
DL.b = np.array(([12.5,30.,2.]))
DL.fGHz=np.linspace(0.9,10,500)

DL.Aa=Antenna(typ='Omni')
DL.Ab=Antenna(typ='Omni')
# In[8]:

plt.ion()

# In[9]:
DL.eval(diffraction=True,force=True,ra_ceil_H=0)


# In[10]:

DL.R


# In[ ]:



