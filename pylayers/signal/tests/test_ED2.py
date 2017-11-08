from pylayers.signal.ED import *
#
fsGHz = 50
fcGHz = 4
BGHz = 3.5
Tns = 5
ed = ED(fsGHz=fsGHz,fcGHz=fcGHz,BGHz=BGHz,Tns=Tns)
# from 3 to 30 dB
SNR = np.arange(5,30,1)
snr = 10**(SNR/10.)
tE = snr*ed.scale
plt.figure(figsize=(10,10))
for E in tE:
    x = ed.roc(E)
    plt.plot(x[1],x[2],label='SNR ='+str(10*np.log10(E/ed.scale))+' dB')
plt.title('ROC curve',fontsize=22)
plt.legend()
plt.xlabel(r'$P_{fa}$',fontsize=18)
plt.ylabel(r'$P_d$',fontsize=18)
