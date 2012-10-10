#
#
# Values after calibration with S2 (3.5GHz)
#
#     WALLS  (PARTITION)   8    - 0.2 j      0.1 m
#     WOODEN_DOOR          2.5  - 0.03 j     0.1 m
#     INSULATION           2    - 0.5 j      0.2 m
#     PLASTERBOARD_14CM    2.5  - 0.3 j      0.14
#     PLASTERBOARD_7CM     2.5  - 0.3 j      0.07
#     CONCRETE_15CM3D      5.5  - 0.25 j     0.15
#     CONCRETE_20CM3D      5.5  - 0.25 j     0.2
#     CONCRETE_7CM3D       5.5  - 0.25 j     0.07
#     CONCRETE_25CM3D      5.5  - 0.25 j     0.25
#     3D_WINDOW_GLASSES    4.5  - 0.03 j     0.05
#     WOODEN_TABLES_DESK2  2.5  - 0.03 j     0.01
#     FLOOR                4.5  - 0.03 j     0.5
#     TOP_CEILING          4.5  - 0.6 j      0.5
#     CEILING              2.5  - 0.3 j      0.005
#     Metallics            1    - 1e6 j      0.005
#     PLASTIC              4    - 0.4 j      0.1
#
    #R1_A=furniture('Table',array([1.33,5.35]),120,0.75,1.4,0.7,0.03,'WOOD')
    #R1_B1=furniture('Desk',array([4.05,7.5]),90,0.75,1.22,0.8,0.03,'WOOD')
    #R1_B2=furniture('Desk',array([4.05,7.1]),0,0.75,0.85,0.4,0.03,'WOOD')
  
    #R2_A=furniture('Block of Tables',array([0.5,13.5]),0,0.75,3.30,1.2,0.03,'WOOD')
  
    #R6_A=furniture('1.3m long table',array([-5.47,8]),90,0.75,1.30,0.65,0.03,'WOOD')
    #R6_B=furniture('1.6m long table',array([-5.32,6]),90,0.75,1.60,0.8,0.03,'WOOD')
    #R6_C=furniture('0.9m high cupboard',array([-0.7,5]),0,0.9,0.60,0.6,0.01,'WOOD')
#    R6_D=furniture('refrigerator', np.array([-0.7,5.6]),0,0.9,0.60,0.6,0.005,'METAL')
#    R6_E=furniture('Drinking fountain', np.array([-0.59,6.2]),0,1.22,0.4,0.35,0.005,'METAL')
#    R6_F=furniture('1.7m high cupboard', np.array([-0.05,7.7]),90,1.7,0.9,0.42,0.005,'METAL')
#    R6_G=furniture('0.55m high printer', np.array([-5.88,6]),0,0.55,0.53,0.5,0.005,'METAL')
#    R6_HA=furniture('copy machine', np.array([-5.76,5]),0,1.17,1.2,0.6,0.005,'METAL')
#    R6_HB=furniture('copy machine', np.array([-5.46,5]),0,1.32,0.9,0.6,0.005,'METAL')
#    R6_IA=furniture('high printer', np.array([-4.3,4.85]),0,0.5,0.75,0.57,0.005,'METAL')
#    R6_IB=furniture('high printer', np.array([-4.3,4.85]),0,1.06,1.0,0.57,0.005,'METAL')
#  
#    R13_A = furniture('1.33-high cupboard', np.array([-6.17,13.05]),0,1.33,0.64,0.4,0.005,'METAL')
#    R13_B = furniture('2m-high cupboard', np.array([-6,13.5]),0,2,1.6,0.45,0.005,'METAL')
#    R13_C = furniture('1.7-high cupboard', np.array([-5.7,15.3]),90,1.7,0.9,0.42,0.005,'METAL')
#  
#  
#    R12_A1= furniture('1.98m high cupboards', np.array([-7.54,9.88]),0,1.98,1.2,0.42,0.005,'METAL')
#    R12_A2= furniture('1.98m high cupboards', np.array([-6.16,6.2]),90,1.98,1.2,0.42,0.005,'METAL')
#    #R12_B= furniture('Table', np.array([-8,7.4]),0,0.75,1.1,1.1,0.03,'WOOD')
#    #R12_C1= furniture('Desk', np.array([-6.17,4.92]),90,0.75,1.22,0.8,0.03,'WOOD')
#    #R12_C2= furniture('Desk', np.array([-7.82,6.16]),0,1.22,0.85,0.40,0.03,'WOOD')
#    R12_D= furniture('1.02m-high cupbpard', np.array([-8.73,7.28]),90,1.02,1.2,0.42,0.005,'METAL')
#  
#    #R7_A1= furniture('Groups of 4 desk', np.array([-10.70,12.4]),90,0.75,3.32,2.82,0.03,'WOOD')
#    #R7_A2= furniture('Groups of 4 desk', np.array([-6.92,12.4]),90,0.75,3.32,2.82,0.03,'WOOD')
#    #R7_B = furniture('Block of 6 desks', np.array([-9.49,5.42]),90,0.75,4.82,3.22,0.03,'WOOD')
#    R7_C1= furniture('1.98m-high cupboards', np.array([-13.71,11.69]),90,1.98,1.20,0.42,0.005,'METAL')
#    R7_C2= furniture('1.98m-high cupboards', np.array([-13.71,8.7]),90,1.98,1.20,0.42,0.005,'METAL')
#    R7_D = furniture('Group of 1.02m-high cupboards', np.array([-13.71,4.98]),90,1.02,3.62,0.42,0.005,'METAL')
#      
#    #R10_A = furniture('1.80m-longtable', np.array([-16.7,5]),0,0.75,1.80,0.80,0.03,'WOOD')
#    #R10_B = furniture('1.60m-longtable', np.array([-14.1,5]),90,0.75,1.60,0.80,0.03,'WOOD')
#    R10_C = furniture('0.70m-high cupboards', np.array([-14.1,8.43]),90,0.80,0.90,0.42,0.005,'METAL')
#    R10_D = furniture('1.98m-high cupboards', np.array([-15.52,9.80]),0,1.98,1.20,0.42,0.005,'METAL')
#      
#    #R11_A = furniture('Desk', np.array([-20.1,6.45]),0,0.75,2.10,1.04,0.03,'WOOD')
#    #R11_B= furniture('Table', np.array([-19.3,8.6]),0,0.75,1.1,1.1,0.03,'WOOD')
#    R11_C= furniture('1.98m-high cupbpard', np.array([-17.15,6.6]),90,1.98,1.2,0.42,0.005,'METAL')
#    R11_D1 = furniture('1.02m-high cupbpard', np.array([-17.15,5.4]),90,1.02,1.2,0.42,0.005,'METAL')
#    R11_D2 = furniture('1.02m-high cupbpard', np.array([-20.68,6.15]),90,1.02,1.2,0.42,0.005,'METAL')
#    R11_E1= furniture('0.70m-high cupbpard', np.array([-17.15,4.6]),90,0.70,0.80,0.42,0.005,'METAL')
#    R11_E2= furniture('0.70m-high cupbpard', np.array([-20.68,5.35]),90,0.70,0.80,0.42,0.005,'METAL')
#    R11_F = furniture('1.01m-high cupbpard', np.array([-20.68,4.64]),90,1.01,0.61,0.43,0.01,'METAL')
#    #R8_A = furniture('Block of tabel', np.array([-21.4,12.6]),0,0.75,5.57,2.30,0.03,'WOOD')
#    #R8_B = furniture('Table', np.array([-21.96,14.38]),90,0.75,1.40,0.70,0.03,'WOOD')
#  
#    #R9_A1 = furniture('Block of desks', np.array([-26.75,12.55]),0,0.75,3.2,3.2,0.03,'WOOD')
#    #R9_A2 = furniture('Block of desks', np.array([-25,5.52]),0,0.75,3.2,3.2,0.03,'WOOD')
#    R9_B1 = furniture('1.98m-high cupboards', np.array([-22.74,12.3]),90,1.98,1.20,0.42,0.005,'METAL')
#    R9_B2 = furniture('1.98m-high cupboards', np.array([-21.15,8.81]),90,1.98,1.20,0.42,0.005,'METAL')
#    R9_B3 = furniture('1.98m-high cupboards', np.array([-21.15,4.7]),90,1.98,1.20,0.42,0.005,'METAL')
#    #R9_C = furniture('Table', np.array([-27.71,11.11]),0,0.75,1.8,0.8,0.03,'WOOD')
#    R9_D = furniture('0.70m-high cupboards',np.array([-23,4.6]),0,0.70,0.80,0.43,0.005,'METAL')
#    R9_E1 = furniture('1.98m-high cupboards',np.array([-25.78,4.7]),0,1.33,0.64,0.40,0.005,'METAL')
#    R9_E2 = furniture('1.98m-high cupboards',np.array([-25.78,9.02]),0,1.33,0.64,0.40,0.005,'METAL')
#    R9_F = furniture('1.04m-high cupboards',np.array([-25.6,9.42]),0,1.04,0.42,0.30,0.005,'METAL')
#    R9_G = furniture('Printer',np.array([-23.29,11.7]),0,0.90,0.55,0.40,0.005,'METAL')
#  
#    Fur1 = R6_E
#    Fur2 = R6_F
#    Fur3 = R6_HB
#    Fur4 = R13_A
#    Fur5 = R13_B
#    Fur6 = R13_C
#    Fur7 = R12_A1
#    Fur8 = R12_A2
#    Fur9 = R7_C1
#    Fur10 = R7_C2
#    Fur11 = R10_D
#    Fur12 = R11_C
#    Fur13 = R9_B1
#    Fur14 = R9_B2
#    Fur15 = R9_B3
#    Fur16  = R9_E1
#    Fur17 = R9_E2
#  
#    FL = [Fur1,Fur2,Fur3,Fur4,Fur5,Fur6,Fur7,Fur8,Fur9,Fur10,Fur11,Fur12,Fur13,Fur14,Fur15,Fur16,Fur17]

#    ###
#    ### Dictionary of  the room-point
#    ###
#    # .. todo:: automatiser cela 
#
#    Group = {}
#    Group['R1']=np.array(range(297,333,1))
#    Group['R2']=np.array(range(333,377,1))
#    Group['R6']=np.array(range(255,297,1))
#    Group['R7']=np.array(range(133,232,1))
#    Group['R8']=np.array(range(82,122,1))
#    Group['R9']=np.hstack((np.array(range(1,82,1)),np.array(range(122,133,1))))
#    Group['R10']=np.array([])
#    Group['R11']=np.array([])
#    Group['R12']=np.array([])
#    Group['R13']=np.array([])
#    Group['Couloir right']=np.array(range(232,255,1))
#  
#    S   = Simul.Simul()
#    S.layout('sircut.str','matDB.mat','slabDB.ini')
#    S.load
#    Tx,Rx = ptSiradel()
#    #S.load('simul-siradel')
#    S.L.display['Node']=False
#    S.L.display['Scaled']=False
#    S.L.display['Thin']=True
#  
#
#    b    = S.L.boundary()
#    #zone = {}
#    #zone['xmin'] =b[0,0]
#    #zone['xmax'] =b[0,1]
#    #zone['ymin'] =b[1,0]
#    #zone['ymax'] =b[1,1]
#    #zone['Nx']   = 75
#    #zone['Ny']   = 25
#
#    #D=io.loadmat('M1-h1.mat')
#    #D = io.loadmat('newM1h1.mat')
#    D = io.loadmat('M1_essai_new.mat')
#    #
#    # Deal with Matlab file problem
#    #
#    #    Txs  = D['M1h1'][0][0].pTx
#    #    Emax = D['M1h1'][0][0].Emax
#    #    Etot = D['M1h1'][0][0].Etot
#    #    tau_rms =  D['M1h1'][0][0].tau_rms
#    #    tau_m =  D['M1h1'][0][0].tau_moy
#    #    toa_max =  D['M1h1'][0][0].toa_max
#    #    toa_th =  D['M1h1'][0][0].toa_th
#    #    toa_cum =  D['M1h1'][0][0].toa_cum
#    #    toa_new2 = D['M1h1'][0][0].toa_new2
#    #    tau0    =  D['M1h1'][0][0].distance/0.3
#    Txs  = D['M1']['pTx'][0][0]
#    Emax = D['M1']['Emax'][0][0]
#    Etot = D['M1']['Etot'][0][0]
#    Etau0 = D['M1']['Etau0'][0][0]
#    Efirst = D['M1']['Efirst'][0][0]
#    tau_rms =  D['M1']['tau_rms'][0][0]
#    tau_m =  D['M1']['tau_moy'][0][0]
#    toa_max =  D['M1']['toa_max'][0][0]
#    toa_th =  D['M1']['toa_th'][0][0]
#    toa_cum =  D['M1']['toa_cum'][0][0]
#    #toa_new2 = D['M1'][0][0].toa_new2
#    toa_win = D['M1']['toa_win'][0][0]
#    tau0    =  D['M1']['distance'][0][0]/0.3
#    LQI    = D['M1']['LQI1'][0][0]
#    index_conv = D['M1']['index_conv'][0][0]
#    los_cond = D['M1']['los_cond'][0][0][0][0]  
#  
#    u   = Txs[0:2,:].T
#    nd1 = find(LQI[0]<=10)
#    nd2 = find(LQI[1]<=10)
#    nd3 = find(LQI[2]<=10)
#    nd4 = find(LQI[3]<=10)
#
#    nu1 = find(LQI[0]>10)
#    nu2 = find(LQI[1]>10)
#    nu3 = find(LQI[2]>10)
#    nu4 = find(LQI[3]>10)
#  
#    ###
#    ### find the number of the point whose LQI > 10 for all the 4 receptions
#    ###
#    nu12   = np.array([val for val in nu1 if val in nu2])
#    nu34   = np.array([val for val in nu3 if val in nu4])
#    nu1234 = np.array([val for val in nu12 if val in nu34])
#  
#    etoa_max = abs(toa_max - tau0)
#    etoa_cum = abs(toa_cum - tau0)
#    etoa_th  = abs(toa_th - tau0)
#    #etoa_new2 = abs(toa_new2 - tau0)
#    etoa_win = abs(toa_win - tau0)  
#
#    rx1los = los_cond['Rx1'][0][0]['los']-1
#    rx2los = los_cond['Rx2'][0][0]['los']-1
#    rx3los = los_cond['Rx3'][0][0]['los']-1
#    rx4los = los_cond['Rx4'][0][0]['los']-1
#
#    rx1nlos = los_cond['Rx1'][0][0]['nlos']-1
#    rx2nlos = los_cond['Rx2'][0][0]['nlos']-1
#    rx3nlos = los_cond['Rx3'][0][0]['nlos']-1
#    rx4nlos = los_cond['Rx4'][0][0]['nlos']-1
#  
#    rx1nlos2 = los_cond['Rx1'][0][0]['nlos2']-1
#    rx2nlos2 = los_cond['Rx2'][0][0]['nlos2']-1
#    rx3nlos2 = los_cond['Rx3'][0][0]['nlos2']-1
#    rx4nlos2 = los_cond['Rx4'][0][0]['nlos2']-1
#
#    Urx1los   = np.array([val for val in nu1234 if val in rx1los])
#    Urx1nlos  = np.array([val for val in nu1234 if val in rx1nlos])
#    Urx1nlos2 = np.array([val for val in nu1234 if val in rx1nlos2])
#  
#    Urx2los   = np.array([val for val in nu1234 if val in rx2los])
#    Urx2nlos  = np.array([val for val in nu1234 if val in rx2nlos])
#    Urx2nlos2 = np.array([val for val in nu1234 if val in rx2nlos2])
#
#    Urx3los   = np.array([val for val in nu1234 if val in rx3los])
#    Urx3nlos  = np.array([val for val in nu1234 if val in rx3nlos])
#    Urx3nlos2 = np.array([val for val in nu1234 if val in rx3nlos2])
#
#    Urx4los   = np.array([val for val in nu1234 if val in rx4los])
#    Urx4nlos  = np.array([val for val in nu1234 if val in rx4nlos])
#    Urx4nlos2 = np.array([val for val in nu1234 if val in rx4nlos2])
####
#### Modèle MultiWALL Dans Le Batiment Siradel
####
#
#    Emax1 = Emax[0,:]
#    Emax2 = Emax[1,:]
#    Emax3 = Emax[2,:]
#    Emax4 = Emax[3,:]
#
#    Etot1 = Etot[0,:]
#    Etot2 = Etot[1,:]
#    Etot3 = Etot[2,:]
#    Etot4 = Etot[3,:]
#
#    Etau01 = Etau0[0,:]
#    Etau02 = Etau0[1,:]
#    Etau03 = Etau0[2,:]
#    Etau04 = Etau0[3,:]
#
#    Efirst1 = Efirst[0,:]
#    Efirst2 = Efirst[1,:]
#    Efirst3 = Efirst[2,:]
#    Efirst4 = Efirst[3,:]
#
#    Rx1 = Rx[1,0:2]
#    Rx2 = Rx[2,0:2]
#    Rx3 = Rx[3,0:2]
#    Rx4 = Rx[4,0:2]
#
#    D21  = Dgrid_points(u,Rx1)
#    D22  = Dgrid_points(u,Rx2)
#    D23  = Dgrid_points(u,Rx3)
#    D24  = Dgrid_points(u,Rx4)
#
##
## set frequency 
##
#
#    f      = 5.07
#
#    PL1    = OneSlopeMdl(D21,2,f)
#    PL2    = OneSlopeMdl(D22,2,f)
#    PL3    = OneSlopeMdl(D23,2,f)
#    PL4    = OneSlopeMdl(D24,2,f)
#
#    Lw1o,Lw1p  = Loss0_v2(S,u,f,Rx1)
#    Lw2o,Lw2p  = Loss0_v2(S,u,f,Rx2)
#    Lw3o,Lw3p  = Loss0_v2(S,u,f,Rx3)
#    Lw4o,Lw4p  = Loss0_v2(S,u,f,Rx4)
#
#
#
#
####
####  Calculate the diffraction loss due to the obstacle metal
####
#
#    LM = np.array([])
#
#    for i in range(302):
#        for j in range(4):
#            xt = u[i][0]
#            yt = u[i][1]
#            xr = Rx[:,0:2][j+1][0]
#            yr = Rx[:,0:2][j+1][1]
#            for k in range(17):
#                Slist = Interline(xt,yt,xr,yr,FL[k])
#                Lm = Loss_obstacle(Slist,xt,yt,xr,yr,FL[k])
#                LM = hstack((LM,Lm))
#
#    LI = np.array([])
#    LJ = np.array([])
#    NU = np.array([])
#    MM = LM.reshape(302,68)
#    for i in range (302):
#        G = MM[i].reshape(4,17)
#        for j in range(4):
#            Nu = find(G[j]!=0)
#            if len(Nu) >1:
#                print i,j,Nu
#                LI = hstack((LI,i))
#                LJ = hstack((LJ,j))
#                NU = hstack((NU,Nu))
#  
#    LT = np.array([])
#    for i in range(len(LI)):
#      
#      
#        xt = u[LI[i]][0]
#        yt = u[LI[i]][1]
#        xr = Rx[:,0:2][LJ[i]+1][0]
#        yr = Rx[:,0:2][LJ[i]+1][1]
#      
#
#        Lt = Loss_2obstacle(xt,yt,xr,yr,FL[int(NU[(2*i)])],FL[int(NU[(2*i+1)])])  
#        LT = np.hstack((LT,Lt))
#
#    Loss_lien = np.array([])
#    for i in range(302):
#        LSR = MM[i].reshape(4,17)
#        Ls = sum(LSR,axis=1)
#        Loss_lien = hstack((Loss_lien,Ls))
#      
#    Loss_lien = Loss_lien.reshape(302,4)
#  
#    for i in range(56):
#      
#        Loss_lien[LI[i],LJ[i]]= LT[i]
#
#    traj = np.array([1,2,3,5,6,7,8,9,10,11,14,19,20,37,38,43,44,49,50,55,56,61,62,67,68,83,
#        82,81,80,98,99,102,103,107,109,108,106,104,101,100,97,98,80,79,114,115,
#        119,121,123,131,135,140,143,148,147,146,153,154,159,158,157,168,174,178,
#        181,186,188,191,192,195,196,199,200,203,205,206,212,213,214,217,218,221,
#        222,225,226,227,224,223,220,219,216,215,213,212,209,210,228,229,232,234,
#        236,238,240,242,244,243,241,239,237,235,233,231,230,211,208,207,204,202,
#        201,198,197,194,193,190,189,185,182,177,174,169,166,161,151,147,144,139,
#        136,130,124,120,119,116,112,77,76,71,70,65,64,59,58,53,52,47,46,41,40,35,
#        34,33,32,31,30,29,26,25,24,23,22,21,36,20,19,18,17,16])
#
#
#
####
#### Pass loss calculate
####      
#
##        Pt    =  6.368
#    Pt = 0
#
#    Pregle = -1.7227573189526098
#    """
#    Power1o = Pt - PL1 - Lw1o
#    Power2o = Pt - PL2 - Lw2o
#    Power3o = Pt - PL3 - Lw3o
#    Power4o = Pt - PL4 - Lw4o
#    """
#    Power1o = Pt - PL1 - Lw1o - Loss_lien[:,0]
#    Power2o = Pt - PL2 - Lw2o - Loss_lien[:,1]
#    Power3o = Pt - PL3 - Lw3o - Loss_lien[:,2]
#    Power4o = Pt - PL4 - Lw4o - Loss_lien[:,3]
#    """  
#    Power1p = Pt - PL1 - Lw1p
#    Power2p = Pt - PL2 - Lw2p
#    Power3p = Pt - PL3 - Lw3p
#    Power4p = Pt - PL4 - Lw4p
#    """
#    Power1p = Pt - PL1 - Lw1p - Loss_lien[:,0]+Pregle
#    Power2p = Pt - PL2 - Lw2p - Loss_lien[:,1]+Pregle
#    Power3p = Pt - PL3 - Lw3p - Loss_lien[:,2]+Pregle
#    Power4p = Pt - PL4 - Lw4p - Loss_lien[:,3]+Pregle
#  
#    plot(Power1p[traj-1])
#    plot(Power2p[traj-1])
#    plot(Power3p[traj-1])
#    plot(Power4p[traj-1])
#
#    eMW_Ef1 = Power1p - Efirst1
#    eMW_Ef2 = Power2p - Efirst2
#    eMW_Ef3 = Power3p - Efirst3
#    eMW_Ef4 = Power4p - Efirst4
#
#    ### reglementation -1.7228dB
#
#    REG = -1.7228
#
#    OSM1 =  - PL1 + REG
#    OSM2 =  - PL2 + REG
#    OSM3 =  - PL3 + REG
#    OSM4 =  - PL4 + REG
#
#    MWM1 = -PL1-Lw1p + REG
#    MWM2 = -PL2-Lw2p + REG
#    MWM3 = -PL3-Lw3p + REG
#    MWM4 = -PL4-Lw4p + REG
#
#    MWM21 = -PL1-Lw1p-Loss_lien[:,0] + REG
#    MWM22 = -PL2-Lw2p-Loss_lien[:,1] + REG
#    MWM23 = -PL3-Lw3p-Loss_lien[:,2] + REG
#    MWM24 = -PL4-Lw4p-Loss_lien[:,3] + REG
#
#  
#    OSM_moy = mean(hstack((abs((OSM1 - Efirst1)[nu1234]),abs((OSM2 - Efirst2)[nu1234]),abs((OSM3 - Efirst3)[nu1234]),abs((OSM4 - Efirst4)[nu1234]))))
#    OSM_std = std(hstack((abs((OSM1 - Efirst1)[nu1234]),abs((OSM2 - Efirst2)[nu1234]),abs((OSM3 - Efirst3)[nu1234]),abs((OSM4 - Efirst4)[nu1234]))))  
#    MWM_moy = mean(hstack((abs((MWM1 - Efirst1)[nu1234]),abs((MWM2 - Efirst2)[nu1234]),abs((MWM3 - Efirst3)[nu1234]),abs((MWM4 - Efirst4)[nu1234]))))
#    MWM_std = std(hstack((abs((MWM1 - Efirst1)[nu1234]),abs((MWM2 - Efirst2)[nu1234]),abs((MWM3 - Efirst3)[nu1234]),abs((MWM4 - Efirst4)[nu1234]))))  
#    MWM2_moy = mean(hstack((abs((MWM21 - Efirst1)[nu1234]),abs((MWM22 - Efirst2)[nu1234]),abs((MWM23 - Efirst3)[nu1234]),abs((MWM24 - Efirst4)[nu1234]))))
#    MWM2_std = std(hstack((abs((MWM21 - Efirst1)[nu1234]),abs((MWM22 - Efirst2)[nu1234]),abs((MWM23 - Efirst3)[nu1234]),abs((MWM24 - Efirst4)[nu1234]))))  
#    """
#  
#    OSM_moy = mean(hstack(((OSM1 - Efirst1)[nu1234],(OSM2 - Efirst2)[nu1234],(OSM3 - Efirst3)[nu1234],(OSM4 - Efirst4)[nu1234])))
#    :x
#    MWM_moy = mean(hstack(((MWM1 - Efirst1)[nu1234],(MWM2 - Efirst2)[nu1234],(MWM3 - Efirst3)[nu1234],(MWM4 - Efirst4)[nu1234])))
#    MWM_std = std(hstack(((MWM1 - Efirst1)[nu1234],(MWM2 - Efirst2)[nu1234],(MWM3 - Efirst3)[nu1234],(MWM4 - Efirst4)[nu1234])))  
#    MWM2_moy = mean(hstack(((MWM21 - Efirst1)[nu1234],(MWM22 - Efirst2)[nu1234],(MWM23 - Efirst3)[nu1234],(MWM24 - Efirst4)[nu1234])))
#    MWM2_std = std(hstack(((MWM21 - Efirst1)[nu1234],(MWM22 - Efirst2)[nu1234],(MWM23 - Efirst3)[nu1234],(MWM24 - Efirst4)[nu1234])))  
#    """
#    """
##
## export the results of OSM, MWM and MWM2  
##
#
#    F = {}
#    F['OSM']  = np.array([])
#    F['MWM']  = np.array([])
#    F['MWM2'] = np.array([])
#    F['x_Tx'] = np.array([])
#    F['y_Tx'] = np.array([])
#    F['x_Rx'] = np.array([])
#    F['y_Rx'] = np.array([])
#    F['lqi'] = np.array([])  
#    F['id'] = np.array([])  
#    F['Efirst'] = np.array([])
#    F['distance'] = np.array([])  
#    for i in range(302):
#        F['OSM'] = np.hstack((F['OSM'],np.array([OSM1[i],OSM2[i],OSM3[i],OSM4[i]])))
#        F['MWM'] = np.hstack((F['MWM'],np.array([MWM1[i],MWM2[i],MWM3[i],MWM4[i]])))
#        F['MWM2'] = np.hstack((F['MWM2'],np.array([MWM21[i],MWM22[i],MWM23[i],MWM24[i]])))
#        F['x_Tx'] = np.hstack((F['x_Tx'],np.array([u[i][0],u[i][0],u[i][0],u[i][0]])))
#        F['y_Tx'] = np.hstack((F['y_Tx'],np.array([u[i][1],u[i][1],u[i][1],u[i][1]])))
#        F['x_Rx'] = np.hstack((F['x_Rx'],np.array([Rx[1][0],Rx[2][0],Rx[3][0],Rx[4][0]])))
#        F['y_Rx'] = np.hstack((F['y_Rx'],np.array([Rx[1][1],Rx[2][1],Rx[3][1],Rx[4][1]])))
#        F['lqi'] = np.hstack((F['lqi'],np.array([LQI[0][i],LQI[1][i],LQI[2][i],LQI[3][i]])))
#        F['id'] = np.hstack((F['id'],np.array([1,2,3,4])))
#        F['Efirst'] = hstack((F['Efirst'],np.array([Efirst1[i],Efirst2[i],Efirst3[i],Efirst4[i]])))
#        F['distance'] = hstack((F['distance'],np.array([tau0[0][i]*0.3,tau0[1][i]*0.3,tau0[2][i]*0.3,tau0[3][i]*0.3])))
#      
#    io.savemat('M1_h1_Multiwall',F)
#
#"""
#  
#  
