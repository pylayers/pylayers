#!/usr/bin/python
#-*- coding: utf-8 -*-
"""
.. currentmodule:: pylayers.measures.mesuwb


RAW_DATA Class
==============

.. autosummary::
    :toctree: generated/

     RAW_DATA.__init__

CAL_DATA Class
==============

.. autosummary::
    :toctree: generated/

    CAL_DATA.__init__
    CAL_DATA.plot
    CAL_DATA.getwave

Fdd Class
=========

.. autosummary::
    :toctree: generated/

     Fdd.__init__
     Fdd.plot

Tdd Class
=========

.. autosummary::
    :toctree: generated/

     Tdd.__init__
     Tdd.PL
     Tdd.show
     Tdd.show_span
     Tdd.box
     Tdd.plot

TFP Class
=========

.. autosummary::
    :toctree: generated/

     TFP.__init__
     TFP.append

FP Class
========

.. autosummary::
    :toctree: generated/

     FP.__init__

UWBMeasure
==========

.. autosummary::
    :toctree: generated/

    UWBMeasure.__init__
    UWBMeasure.info
    UWBMeasure.show
    UWBMeasure.Epercent
    UWBMeasure.toa_max2
    UWBMeasure.tau_Emax
    UWBMeasure.tau_moy
    UWBMeasure.tau_rms
    UWBMeasure.toa_new
    UWBMeasure.toa_win
    UWBMeasure.toa_max
    UWBMeasure.toa_th
    UWBMeasure.toa_cum
    UWBMeasure.taumax
    UWBMeasure.Emax
    UWBMeasure.Etot
    UWBMeasure.Efirst
    UWBMeasure.Etau0
    UWBMeasure.ecdf
    UWBMeasure.tdelay
    UWBMeasure.fp
    UWBMeasure.outlatex


Utility Functions
=================

.. autosummary::
    :toctree:

    mesname
    ptw1
    visibility
    trait

"""
import pdb
import doctest
import os.path
import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import pylayers.util.pyutil as pyu
import pylayers.util.geomutil as geu
import pylayers.signal.bsignal as bs
from   pylayers.util.project import *
from   scipy import io
from   scipy import signal, linspace, polyval, polyfit, stats
#from   pylayers.Simul import SimulEM
#
# Utility functions
#
def mesname(n, dirname, h=1):
    """ get the measurement data filename

    Parameters
    ----------
    n
         Tx point number
    dirname
        PATH to measurements directory
    h
        height index 1|2

    Notes
    -----
        This function handle filename from WHERE1 M1 measurement campaign
    """
    if (h == 1):
        mesdir = dirname + '/h1/'
        if n > 279:
            prefix = 'SIRADEL_08-08-01_P'
        else:
            prefix = 'SIRADEL_08-07-31_P'
    else:
        mesdir = dirname + '/h2/'
        prefix = 'SIRADEL_08-08-02_P'

    stn = str(n)
    if (len(stn) == 1):
        stn = '00' + stn
    if (len(stn) == 2):
        stn = '0' + stn

    filename = mesdir + prefix + stn + '.mat'
    return(filename)

def ptw1():
    """ return W1 Tx and Rx points
    """

    p0 = np.array([298626.5107, 2354073.213])

    Rx = np.zeros((5, 2))

    Rx[1, :] = np.array([298614.2383, 2354080.9762])
    Rx[2, :] = np.array([298607.736, 2354088.391])
    Rx[3, :] = np.array([298622.3689, 2354082.0733])
    Rx[4, :] = np.array([298617.4193, 2354088.4029])
    

    Rx[1, :] = Rx[1, :] - p0
    Rx[2, :] = Rx[2, :] - p0
    Rx[3, :] = Rx[3, :] - p0
    Rx[4, :] = Rx[4, :] - p0

    height = 1.2 * np.ones((5, 1))
    Rx = np.hstack((Rx, height))

    Tx = np.zeros((377, 2))

    Tx[1, :] = np.array([298599.1538, 2354087.9963]) - p0
    Tx[2, :] = np.array([298599.1662, 2354086.9940]) - p0
    Tx[3, :] = np.array([298599.1544, 2354085.9948]) - p0
    Tx[4, :] = np.array([298599.1397, 2354085.4917]) - p0
    Tx[5, :] = np.array([298599.6551, 2354085.4952]) - p0
    Tx[6, :] = np.array([298599.6567, 2354085.9977]) - p0
    Tx[7, :] = np.array([298600.1580, 2354086.0006]) - p0
    Tx[8, :] = np.array([298600.1656, 2354085.4991]) - p0
    Tx[9, :] = np.array([298600.6429, 2354085.5026]) - p0
    Tx[10, :] = np.array([298601.1453, 2354085.5170]) - p0
    Tx[11, :] = np.array([298601.6437, 2354085.5227]) - p0
    Tx[12, :] = np.array([298602.1446, 2354085.5180]) - p0
    Tx[13, :] = np.array([298602.6474, 2354085.5145]) - p0
    Tx[14, :] = np.array([298602.6432, 2354085.0153]) - p0
    Tx[15, :] = np.array([298602.1375, 2354085.0208]) - p0
    Tx[16, :] = np.array([298601.6324, 2354085.0098]) - p0
    Tx[17, :] = np.array([298601.1347, 2354084.9793]) - p0
    Tx[18, :] = np.array([298601.1560, 2354084.4782]) - p0
    Tx[19, :] = np.array([298601.6638, 2354084.5002]) - p0
    Tx[20, :] = np.array([298602.1413, 2354084.5114]) - p0
    Tx[21, :] = np.array([298602.6510, 2354084.5142]) - p0
    Tx[22, :] = np.array([298602.6596, 2354084.0093]) - p0
    Tx[23, :] = np.array([298602.1531, 2354083.9971]) - p0
    Tx[24, :] = np.array([298601.6619, 2354083.9895]) - p0
    Tx[25, :] = np.array([298601.1702, 2354083.9752]) - p0
    Tx[26, :] = np.array([298601.1892, 2354083.4697]) - p0
    Tx[27, :] = np.array([298601.6879, 2354083.4708]) - p0
    Tx[28, :] = np.array([298602.1734, 2354083.4953]) - p0
    Tx[29, :] = np.array([298602.1796, 2354083.0027]) - p0
    Tx[30, :] = np.array([298601.6950, 2354082.9746]) - p0
    Tx[31, :] = np.array([298601.1908, 2354082.9642]) - p0
    Tx[32, :] = np.array([298601.6730, 2354082.4642]) - p0
    Tx[33, :] = np.array([298602.1897, 2354082.5089]) - p0
    Tx[34, :] = np.array([298602.6956, 2354082.5094]) - p0
    Tx[35, :] = np.array([298602.6820, 2354083.0076]) - p0
    Tx[36, :] = np.array([298602.6668, 2354083.5180]) - p0
    Tx[37, :] = np.array([298603.1766, 2354083.5187]) - p0
    Tx[38, :] = np.array([298603.1619, 2354084.0124]) - p0
    Tx[39, :] = np.array([298603.1435, 2354084.5246]) - p0
    Tx[40, :] = np.array([298603.6412, 2354084.5477]) - p0
    Tx[41, :] = np.array([298603.6539, 2354084.0204]) - p0
    Tx[42, :] = np.array([298603.6735, 2354083.5422]) - p0
    Tx[43, :] = np.array([298604.1769, 2354083.5510]) - p0
    Tx[44, :] = np.array([298604.1560, 2354084.0272]) - p0
    Tx[45, :] = np.array([298604.1398, 2354084.5635]) - p0
    Tx[46, :] = np.array([298604.6394, 2354084.5710]) - p0
    Tx[47, :] = np.array([298604.6531, 2354084.0357]) - p0
    Tx[48, :] = np.array([298604.6723, 2354083.5647]) - p0
    Tx[49, :] = np.array([298605.1724, 2354083.5727]) - p0
    Tx[50, :] = np.array([298605.1552, 2354084.0452]) - p0
    Tx[51, :] = np.array([298605.1411, 2354084.5744]) - p0
    Tx[52, :] = np.array([298605.6407, 2354084.5896]) - p0
    Tx[53, :] = np.array([298605.6564, 2354084.0518]) - p0
    Tx[54, :] = np.array([298605.6725, 2354083.5854]) - p0
    Tx[55, :] = np.array([298606.1697, 2354083.5946]) - p0
    Tx[56, :] = np.array([298606.1504, 2354084.0704]) - p0
    Tx[57, :] = np.array([298606.1394, 2354084.5736]) - p0
    Tx[58, :] = np.array([298606.6116, 2354084.5855]) - p0
    Tx[59, :] = np.array([298606.6557, 2354084.0949]) - p0
    Tx[60, :] = np.array([298606.6505, 2354083.6042]) - p0
    Tx[61, :] = np.array([298607.1728, 2354083.6169]) - p0
    Tx[62, :] = np.array([298607.1490, 2354084.1052]) - p0
    Tx[63, :] = np.array([298607.1420, 2354084.5882]) - p0
    Tx[64, :] = np.array([298607.6405, 2354084.6156]) - p0
    Tx[65, :] = np.array([298607.6537, 2354084.1195]) - p0
    Tx[66, :] = np.array([298607.6718, 2354083.6162]) - p0
    Tx[67, :] = np.array([298608.1705, 2354083.6262]) - p0
    Tx[68, :] = np.array([298608.1498, 2354084.1300]) - p0
    Tx[69, :] = np.array([298608.1434, 2354084.6169]) - p0
    Tx[70, :] = np.array([298608.6389, 2354084.6207]) - p0
    Tx[71, :] = np.array([298608.6529, 2354084.1442]) - p0
    Tx[72, :] = np.array([298608.6908, 2354083.6284]) - p0
    Tx[73, :] = np.array([298609.1579, 2354083.6149]) - p0
    Tx[74, :] = np.array([298609.1469, 2354084.1487]) - p0
    Tx[75, :] = np.array([298609.1397, 2354084.6209]) - p0
    Tx[76, :] = np.array([298609.6398, 2354084.6236]) - p0
    Tx[77, :] = np.array([298609.6552, 2354084.1487]) - p0
    Tx[78, :] = np.array([298609.6819, 2354083.6253]) - p0
    Tx[79, :] = np.array([298610.1816, 2354083.6293]) - p0
    Tx[80, :] = np.array([298610.1512, 2354084.1489]) - p0
    Tx[81, :] = np.array([298610.1396, 2354084.6229]) - p0
    Tx[82, :] = np.array([298610.1414, 2354085.1562]) - p0
    Tx[83, :] = np.array([298609.6353, 2354085.1503]) - p0
    Tx[84, :] = np.array([298609.1376, 2354085.1306]) - p0
    Tx[85, :] = np.array([298608.6384, 2354085.1249]) - p0
    Tx[86, :] = np.array([298608.1365, 2354085.1167]) - p0
    Tx[87, :] = np.array([298607.6320, 2354085.1233]) - p0
    Tx[88, :] = np.array([298607.1375, 2354085.1285]) - p0
    Tx[89, :] = np.array([298606.6346, 2354085.1308]) - p0
    Tx[90, :] = np.array([298606.1340, 2354085.1187]) - p0
    Tx[91, :] = np.array([298605.6357, 2354085.1098]) - p0
    Tx[92, :] = np.array([298605.1339, 2354085.1129]) - p0
    Tx[93, :] = np.array([298605.1224, 2354085.6164]) - p0
    Tx[94, :] = np.array([298604.6349, 2354085.0960]) - p0
    Tx[95, :] = np.array([298604.6356, 2354085.6111]) - p0
    Tx[96, :] = np.array([298604.1292, 2354085.5975]) - p0
    Tx[97, :] = np.array([298604.1243, 2354086.1044]) - p0
    Tx[98, :] = np.array([298604.6182, 2354086.1130]) - p0
    Tx[99, :] = np.array([298604.6076, 2354086.6167]) - p0
    Tx[100, :] = np.array([298604.1310, 2354086.6027]) - p0
    Tx[101, :] = np.array([298604.1401, 2354087.0989]) - p0
    Tx[102, :] = np.array([298604.5986, 2354087.1199]) - p0
    Tx[103, :] = np.array([298611.6833, 2354085.6072]) - p0
    Tx[104, :] = np.array([298611.1745, 2354085.6015]) - p0
    Tx[105, :] = np.array([298611.1425, 2354086.0848]) - p0
    Tx[106, :] = np.array([298611.6452, 2354086.1118]) - p0
    Tx[107, :] = np.array([298611.6090, 2354086.6144]) - p0
    Tx[108, :] = np.array([298611.1054, 2354086.6012]) - p0
    Tx[109, :] = np.array([298611.1042, 2354087.1082]) - p0
    Tx[110, :] = np.array([298611.6032, 2354087.1149]) - p0
    Tx[111, :] = np.array([298612.1175, 2354087.1200]) - p0
    Tx[112, :] = np.array([298612.1255, 2354087.6165]) - p0
    Tx[113, :] = np.array([298611.6215, 2354087.6255]) - p0
    Tx[114, :] = np.array([298611.1133, 2354087.6264]) - p0
    Tx[115, :] = np.array([298612.1141, 2354088.1246]) - p0
    Tx[116, :] = np.array([298611.6216, 2354088.1302]) - p0
    Tx[117, :] = np.array([298611.1070, 2354088.1211]) - p0
    Tx[118, :] = np.array([298610.6204, 2354088.1582]) - p0
    Tx[119, :] = np.array([298610.6092, 2354088.6707]) - p0
    Tx[120, :] = np.array([298611.1083, 2354088.6447]) - p0
    Tx[121, :] = np.array([298611.6160, 2354088.6395]) - p0
    Tx[122, :] = np.array([298610.6813, 2354083.6292]) - p0
    Tx[123, :] = np.array([298610.6475, 2354084.1563]) - p0
    Tx[124, :] = np.array([298610.6448, 2354084.6357]) - p0
    Tx[125, :] = np.array([298611.1846, 2354084.6620]) - p0
    Tx[126, :] = np.array([298611.1994, 2354084.1390]) - p0
    Tx[127, :] = np.array([298611.1823, 2354083.6428]) - p0
    Tx[128, :] = np.array([298611.6880, 2354083.6406]) - p0
    Tx[129, :] = np.array([298611.6515, 2354084.1711]) - p0
    Tx[130, :] = np.array([298611.6407, 2354084.6489]) - p0
    Tx[131, :] = np.array([298612.1390, 2354084.6449]) - p0
    Tx[132, :] = np.array([298612.1540, 2354084.1705]) - p0
    Tx[133, :] = np.array([298612.6359, 2354084.6577]) - p0
    Tx[134, :] = np.array([298612.6395, 2354084.1697]) - p0
    Tx[135, :] = np.array([298612.6370, 2354083.6538]) - p0
    Tx[136, :] = np.array([298613.1405, 2354083.6868]) - p0
    Tx[137, :] = np.array([298613.1406, 2354084.1560]) - p0
    Tx[138, :] = np.array([298613.1327, 2354084.6577]) - p0
    Tx[139, :] = np.array([298613.1366, 2354085.1673]) - p0
    Tx[140, :] = np.array([298613.1184, 2354085.6717]) - p0
    Tx[141, :] = np.array([298613.1142, 2354086.1763]) - p0
    Tx[142, :] = np.array([298613.6457, 2354085.6697]) - p0
    Tx[143, :] = np.array([298613.6441, 2354085.1596]) - p0
    Tx[144, :] = np.array([298613.6326, 2354084.6580]) - p0
    Tx[145, :] = np.array([298613.6404, 2354084.1482]) - p0
    Tx[146, :] = np.array([298613.6499, 2354083.6886]) - p0
    Tx[147, :] = np.array([298613.1379, 2354083.1766]) - p0
    Tx[148, :] = np.array([298613.1497, 2354082.6790]) - p0
    Tx[149, :] = np.array([298613.1492, 2354082.1705]) - p0
    Tx[150, :] = np.array([298613.6630, 2354082.1454]) - p0
    Tx[151, :] = np.array([298614.1583, 2354081.6012]) - p0
    Tx[152, :] = np.array([298613.6540, 2354081.6356]) - p0
    Tx[153, :] = np.array([298613.1439, 2354081.6623]) - p0
    Tx[154, :] = np.array([298613.1530, 2354081.1617]) - p0
    Tx[155, :] = np.array([298613.6526, 2354081.1463]) - p0
    Tx[156, :] = np.array([298613.1480, 2354080.6561]) - p0
    Tx[157, :] = np.array([298613.1395, 2354080.1452]) - p0
    Tx[158, :] = np.array([298613.1540, 2354079.6386]) - p0
    Tx[159, :] = np.array([298613.1505, 2354079.1286]) - p0
    Tx[160, :] = np.array([298613.1415, 2354078.6210]) - p0
    Tx[161, :] = np.array([298613.6557, 2354078.6218]) - p0
    Tx[162, :] = np.array([298614.1436, 2354083.6585]) - p0
    Tx[163, :] = np.array([298614.1426, 2354084.1502]) - p0
    Tx[164, :] = np.array([298614.1348, 2354084.6603]) - p0
    Tx[165, :] = np.array([298614.1271, 2354085.1728]) - p0
    Tx[166, :] = np.array([298614.6380, 2354085.1642]) - p0
    Tx[167, :] = np.array([298614.6388, 2354084.6585]) - p0
    Tx[168, :] = np.array([298614.6445, 2354084.1448]) - p0
    Tx[169, :] = np.array([298614.6489, 2354083.6400]) - p0
    Tx[170, :] = np.array([298615.1552, 2354083.6422]) - p0
    Tx[171, :] = np.array([298615.1377, 2354084.1475]) - p0
    Tx[172, :] = np.array([298615.1324, 2354084.6569]) - p0
    Tx[173, :] = np.array([298615.1282, 2354085.1647]) - p0
    Tx[174, :] = np.array([298615.6277, 2354085.1596]) - p0
    Tx[175, :] = np.array([298615.6353, 2354084.6545]) - p0
    Tx[176, :] = np.array([298615.6411, 2354084.1476]) - p0
    Tx[177, :] = np.array([298615.6535, 2354083.6039]) - p0
    Tx[178, :] = np.array([298616.1607, 2354083.6049]) - p0
    Tx[179, :] = np.array([298616.1477, 2354084.1506]) - p0
    Tx[180, :] = np.array([298616.1348, 2354084.6555]) - p0
    Tx[181, :] = np.array([298616.1199, 2354085.1620]) - p0
    Tx[182, :] = np.array([298615.6512, 2354085.6585]) - p0
    Tx[183, :] = np.array([298615.6378, 2354086.1646]) - p0
    Tx[184, :] = np.array([298616.1548, 2354085.6584]) - p0
    Tx[185, :] = np.array([298616.6587, 2354085.6560]) - p0
    Tx[186, :] = np.array([298617.1646, 2354085.6780]) - p0
    Tx[187, :] = np.array([298617.1634, 2354086.1882]) - p0
    Tx[188, :] = np.array([298616.6570, 2354086.1598]) - p0
    Tx[189, :] = np.array([298616.1520, 2354086.1711]) - p0
    Tx[190, :] = np.array([298616.1759, 2354086.6578]) - p0
    Tx[191, :] = np.array([298616.6791, 2354086.6612]) - p0
    Tx[192, :] = np.array([298616.6891, 2354087.1669]) - p0
    Tx[193, :] = np.array([298616.1839, 2354087.1899]) - p0
    Tx[194, :] = np.array([298616.1933, 2354087.7160]) - p0
    Tx[195, :] = np.array([298616.6962, 2354087.6792]) - p0
    Tx[196, :] = np.array([298616.7053, 2354088.1869]) - p0
    Tx[197, :] = np.array([298616.1990, 2354088.2208]) - p0
    Tx[198, :] = np.array([298616.6932, 2354088.6978]) - p0
    Tx[199, :] = np.array([298616.5761, 2354085.1480]) - p0
    Tx[200, :] = np.array([298616.5881, 2354084.6199]) - p0
    Tx[201, :] = np.array([298616.6404, 2354084.1429]) - p0
    Tx[202, :] = np.array([298616.6655, 2354083.6028]) - p0
    Tx[203, :] = np.array([298617.1700, 2354083.5984]) - p0
    Tx[204, :] = np.array([298617.1451, 2354084.1559]) - p0
    Tx[205, :] = np.array([298617.1360, 2354084.6628]) - p0
    Tx[206, :] = np.array([298617.1415, 2354085.1648]) - p0
    Tx[207, :] = np.array([298617.6331, 2354085.1743]) - p0
    Tx[208, :] = np.array([298617.6413, 2354084.6687]) - p0
    Tx[209, :] = np.array([298617.6418, 2354084.1389]) - p0
    Tx[210, :] = np.array([298617.6736, 2354083.6153]) - p0
    Tx[211, :] = np.array([298618.2032, 2354083.5949]) - p0
    Tx[212, :] = np.array([298618.1399, 2354084.1416]) - p0
    Tx[213, :] = np.array([298618.1397, 2354084.6872]) - p0
    Tx[214, :] = np.array([298618.1316, 2354085.1840]) - p0
    Tx[215, :] = np.array([298618.6486, 2354085.2149]) - p0
    Tx[216, :] = np.array([298618.6310, 2354084.6815]) - p0
    Tx[217, :] = np.array([298618.6407, 2354084.1445]) - p0
    Tx[218, :] = np.array([298618.6744, 2354083.6374]) - p0
    Tx[219, :] = np.array([298619.1722, 2354083.6488]) - p0
    Tx[220, :] = np.array([298619.1376, 2354084.1545]) - p0
    Tx[221, :] = np.array([298619.1287, 2354084.6807]) - p0
    Tx[222, :] = np.array([298619.1330, 2354085.1822]) - p0
    Tx[223, :] = np.array([298619.1273, 2354085.6893]) - p0
    Tx[224, :] = np.array([298619.6187, 2354086.2068]) - p0
    Tx[225, :] = np.array([298619.6339, 2354085.6987]) - p0
    Tx[226, :] = np.array([298619.6466, 2354085.1823]) - p0
    Tx[227, :] = np.array([298619.6353, 2354084.6749]) - p0
    Tx[228, :] = np.array([298619.6417, 2354084.1745]) - p0
    Tx[229, :] = np.array([298619.6846, 2354083.6551]) - p0
    Tx[230, :] = np.array([298620.1394, 2354084.1739]) - p0
    Tx[231, :] = np.array([298620.1309, 2354084.6747]) - p0
    Tx[232, :] = np.array([298620.6375, 2354084.7072]) - p0
    Tx[233, :] = np.array([298620.6535, 2354084.1937]) - p0
    Tx[234, :] = np.array([298621.1562, 2354084.2068]) - p0
    Tx[235, :] = np.array([298621.1364, 2354084.7122]) - p0
    Tx[236, :] = np.array([298621.6385, 2354084.7303]) - p0
    Tx[237, :] = np.array([298621.6651, 2354084.2222]) - p0
    Tx[238, :] = np.array([298622.1601, 2354084.2389]) - p0
    Tx[239, :] = np.array([298622.1370, 2354084.7480]) - p0
    Tx[240, :] = np.array([298622.6385, 2354084.7663]) - p0
    Tx[241, :] = np.array([298622.6600, 2354084.2542]) - p0
    Tx[242, :] = np.array([298623.1648, 2354084.2689]) - p0
    Tx[243, :] = np.array([298623.1408, 2354084.7816]) - p0
    Tx[244, :] = np.array([298623.6409, 2354084.7926]) - p0
    Tx[245, :] = np.array([298623.6678, 2354084.2861]) - p0
    Tx[246, :] = np.array([298624.1386, 2354084.8036]) - p0
    Tx[247, :] = np.array([298624.1635, 2354084.2994]) - p0
    Tx[248, :] = np.array([298624.6636, 2354084.3205]) - p0
    Tx[249, :] = np.array([298624.6389, 2354084.8140]) - p0
    Tx[250, :] = np.array([298625.1366, 2354084.8288]) - p0
    Tx[251, :] = np.array([298625.1661, 2354084.3260]) - p0
    Tx[252, :] = np.array([298625.6665, 2354084.3414]) - p0
    Tx[253, :] = np.array([298625.6354, 2354084.8437]) - p0
    Tx[254, :] = np.array([298626.1636, 2354084.3514]) - p0
    Tx[255, :] = np.array([298623.1350, 2354083.2604]) - p0
    Tx[256, :] = np.array([298623.6459, 2354083.2540]) - p0
    Tx[257, :] = np.array([298624.1560, 2354083.2374]) - p0
    Tx[258, :] = np.array([298624.6663, 2354083.2462]) - p0
    Tx[259, :] = np.array([298624.6826, 2354082.7437]) - p0
    Tx[260, :] = np.array([298624.1749, 2354082.7208]) - p0
    Tx[261, :] = np.array([298623.6692, 2354082.7379]) - p0
    Tx[262, :] = np.array([298623.1639, 2354082.7469]) - p0
    Tx[263, :] = np.array([298623.1689, 2354082.2433]) - p0
    Tx[264, :] = np.array([298623.6787, 2354082.2317]) - p0
    Tx[265, :] = np.array([298624.1985, 2354082.2183]) - p0
    Tx[266, :] = np.array([298624.6866, 2354082.2383]) - p0
    Tx[267, :] = np.array([298625.1949, 2354082.2558]) - p0
    Tx[268, :] = np.array([298625.2143, 2354081.7460]) - p0
    Tx[269, :] = np.array([298624.7009, 2354081.7275]) - p0
    Tx[270, :] = np.array([298624.6954, 2354081.2262]) - p0
    Tx[271, :] = np.array([298625.2236, 2354081.2382]) - p0
    Tx[272, :] = np.array([298625.2363, 2354080.7305]) - p0
    Tx[273, :] = np.array([298624.7062, 2354080.7236]) - p0
    Tx[274, :] = np.array([298624.7133, 2354080.2191]) - p0
    Tx[275, :] = np.array([298625.2658, 2354080.2344]) - p0
    Tx[276, :] = np.array([298625.2716, 2354079.7283]) - p0
    Tx[277, :] = np.array([298624.7364, 2354079.7128]) - p0
    Tx[278, :] = np.array([298624.7416, 2354079.2061]) - p0
    Tx[279, :] = np.array([298625.2491, 2354079.2185]) - p0
    Tx[280, :] = np.array([298623.1637, 2354081.7384]) - p0
    Tx[281, :] = np.array([298622.6586, 2354081.7004]) - p0
    Tx[282, :] = np.array([298622.6515, 2354081.1890]) - p0
    Tx[283, :] = np.array([298622.1315, 2354081.1643]) - p0
    Tx[284, :] = np.array([298621.6273, 2354081.1723]) - p0
    Tx[285, :] = np.array([298621.6316, 2354080.6624]) - p0
    Tx[286, :] = np.array([298622.1501, 2354080.6556]) - p0
    Tx[287, :] = np.array([298622.6541, 2354080.6798]) - p0
    Tx[288, :] = np.array([298622.6727, 2354080.1768]) - p0
    Tx[289, :] = np.array([298622.1474, 2354080.1466]) - p0
    Tx[290, :] = np.array([298621.6630, 2354080.1541]) - p0
    Tx[291, :] = np.array([298621.6689, 2354079.6337]) - p0
    Tx[292, :] = np.array([298622.1509, 2354079.6419]) - p0
    Tx[293, :] = np.array([298622.6775, 2354079.6703]) - p0
    Tx[294, :] = np.array([298622.6726, 2354079.1622]) - p0
    Tx[295, :] = np.array([298622.1617, 2354079.1326]) - p0
    Tx[296, :] = np.array([298621.6529, 2354079.1229]) - p0
    Tx[297, :] = np.array([298626.6606, 2354084.3659]) - p0
    Tx[298, :] = np.array([298627.1570, 2354084.3670]) - p0
    Tx[299, :] = np.array([298627.1370, 2354084.8710]) - p0
    Tx[300, :] = np.array([298627.6319, 2354084.8747]) - p0
    Tx[301, :] = np.array([298627.6622, 2354084.3697]) - p0
    Tx[302, :] = np.array([298628.1326, 2354084.8789]) - p0
    Tx[303, :] = np.array([298628.1463, 2354084.3800]) - p0
    Tx[304, :] = np.array([298628.6313, 2354084.8911]) - p0
    Tx[305, :] = np.array([298628.6615, 2354084.3717]) - p0
    Tx[306, :] = np.array([298629.1306, 2354084.8933]) - p0
    Tx[307, :] = np.array([298629.1615, 2354084.3833]) - p0
    Tx[308, :] = np.array([298629.6307, 2354084.9125]) - p0
    Tx[309, :] = np.array([298629.6646, 2354084.3948]) - p0
    Tx[310, :] = np.array([298630.1340, 2354084.9249]) - p0
    Tx[311, :] = np.array([298630.1594, 2354084.4232]) - p0
    Tx[312, :] = np.array([298630.6318, 2354084.9324]) - p0
    Tx[313, :] = np.array([298630.6612, 2354084.4580]) - p0
    Tx[314, :] = np.array([298631.2187, 2354083.3883]) - p0
    Tx[315, :] = np.array([298631.2218, 2354082.3712]) - p0
    Tx[316, :] = np.array([298630.2159, 2354083.3792]) - p0
    Tx[317, :] = np.array([298629.2029, 2354083.3825]) - p0
    Tx[318, :] = np.array([298628.1967, 2354083.3782]) - p0
    Tx[319, :] = np.array([298627.2335, 2354082.3535]) - p0
    Tx[320, :] = np.array([298628.2356, 2354082.3791]) - p0
    Tx[321, :] = np.array([298629.2348, 2354082.3807]) - p0
    Tx[322, :] = np.array([298629.2770, 2354081.3799]) - p0
    Tx[323, :] = np.array([298628.2694, 2354081.3802]) - p0
    Tx[324, :] = np.array([298627.2663, 2354081.3675]) - p0
    Tx[325, :] = np.array([298628.3160, 2354080.3775]) - p0
    Tx[326, :] = np.array([298629.3220, 2354080.3841]) - p0
    Tx[327, :] = np.array([298630.3279, 2354080.3862]) - p0
    Tx[328, :] = np.array([298630.3704, 2354079.3924]) - p0
    Tx[329, :] = np.array([298629.3623, 2354079.3802]) - p0
    Tx[330, :] = np.array([298628.3342, 2354079.3824]) - p0
    Tx[331, :] = np.array([298628.3749, 2354078.3795]) - p0
    Tx[332, :] = np.array([298629.3982, 2354078.3895]) - p0
    Tx[333, :] = np.array([298630.6250, 2354085.4334]) - p0
    Tx[334, :] = np.array([298630.1148, 2354085.4346]) - p0
    Tx[335, :] = np.array([298629.5977, 2354085.4260]) - p0
    Tx[336, :] = np.array([298629.0921, 2354085.4069]) - p0
    Tx[337, :] = np.array([298628.5832, 2354085.3937]) - p0
    Tx[338, :] = np.array([298628.0748, 2354085.3920]) - p0
    Tx[339, :] = np.array([298627.5699, 2354085.3777]) - p0
    Tx[340, :] = np.array([298627.0581, 2354085.3496]) - p0
    Tx[341, :] = np.array([298627.0440, 2354085.8594]) - p0
    Tx[342, :] = np.array([298627.5491, 2354085.8708]) - p0
    Tx[343, :] = np.array([298628.0570, 2354085.8833]) - p0
    Tx[344, :] = np.array([298628.5701, 2354085.8841]) - p0
    Tx[345, :] = np.array([298629.0761, 2354085.8861]) - p0
    Tx[346, :] = np.array([298629.5895, 2354085.8885]) - p0
    Tx[347, :] = np.array([298630.1016, 2354085.8637]) - p0
    Tx[348, :] = np.array([298630.6107, 2354085.9368]) - p0
    Tx[349, :] = np.array([298630.5969, 2354086.4495]) - p0
    Tx[350, :] = np.array([298630.0616, 2354086.3773]) - p0
    Tx[351, :] = np.array([298629.5689, 2354086.3965]) - p0
    Tx[352, :] = np.array([298629.0546, 2354086.4048]) - p0
    Tx[353, :] = np.array([298628.5498, 2354086.3892]) - p0
    Tx[354, :] = np.array([298628.0564, 2354086.3965]) - p0
    Tx[355, :] = np.array([298627.5584, 2354086.3823]) - p0
    Tx[356, :] = np.array([298627.0596, 2354086.3679]) - p0
    Tx[357, :] = np.array([298627.0814, 2354087.9662]) - p0
    Tx[358, :] = np.array([298627.5875, 2354087.9724]) - p0
    Tx[359, :] = np.array([298630.0659, 2354086.8908]) - p0
    Tx[360, :] = np.array([298630.5976, 2354086.9659]) - p0
    Tx[361, :] = np.array([298630.6216, 2354087.4763]) - p0
    Tx[362, :] = np.array([298630.6231, 2354087.9851]) - p0
    Tx[363, :] = np.array([298630.1070, 2354087.9773]) - p0
    Tx[364, :] = np.array([298630.6289, 2354088.4926]) - p0
    Tx[365, :] = np.array([298630.6374, 2354089.0000]) - p0
    Tx[366, :] = np.array([298630.1284, 2354089.0167]) - p0
    Tx[367, :] = np.array([298630.1131, 2354088.4647]) - p0
    Tx[368, :] = np.array([298629.5976, 2354088.4960]) - p0
    Tx[369, :] = np.array([298629.6156, 2354089.0127]) - p0
    Tx[370, :] = np.array([298629.1024, 2354089.0193]) - p0
    Tx[371, :] = np.array([298629.0885, 2354088.4984]) - p0
    Tx[372, :] = np.array([298628.5778, 2354088.5128]) - p0
    Tx[373, :] = np.array([298628.6028, 2354089.0272]) - p0
    Tx[374, :] = np.array([298628.0993, 2354089.0330]) - p0
    Tx[375, :] = np.array([298628.0736, 2354088.4813]) - p0
    Tx[376, :] = np.array([298627.5638, 2354088.4844]) - p0

    return(Tx, Rx)

def visibility():
    """ determine visibility type of WHERE1 measurements campaign points

    Returns
    -------

    visi
        visibility status of each link

    Warning
    -------
    This should be done automatically in the future

    """

    R1 = []

    for i in range(1, 232):
            R1.append('NLOS2')

    for i in range(1, 24):
            R1.append('NLOS')

    for i in range(1, 43):
            R1.append('LOS')

    for i in range(1, 81):
            R1.append('NLOS2')

    R1[330] = 'NLOS'

    R2 = []

    for i in range(1, 82):
            R2.append('NLOS2')

    for i in range(1, 52):
            R2.append('NLOS')

    for i in range(1, 100):
            R2.append('LOS')

    for i in range(1, 3):
            R2.append('NLOS')

    for i in range(1, 144):
            R2.append('NLOS2')

    R2[121] = 'NLOS2'
    R2[122] = 'NLOS2'
    R2[123] = 'NLOS2'
    R2[124] = 'NLOS2'
    R2[125] = 'NLOS2'
    R2[127] = 'NLOS2'

    R3 = []

    for i in range(1, 106):
            R3.append('NLOS2')

    for i in range(1, 17):
            R3.append('NLOS')

    for i in range(1, 12):
            R3.append('NLOS2')

    for i in range(1, 100):
            R3.append('LOS')

    for i in range(1, 146):
            R3.append('NLOS2')

    R3[107] = 'NLOS2'
    R3[130] = 'NLOS'
    R3[131] = 'NLOS'
    R3[209] = 'NLOS2'
    R3[210] = 'NLOS2'
    R3[216] = 'NLOS2'
    R3[217] = 'NLOS2'
    R3[218] = 'NLOS2'
    R3[219] = 'NLOS2'
    R3[220] = 'NLOS2'
    R3[225] = 'NLOS2'
    R3[226] = 'NLOS2'
    R3[227] = 'NLOS2'
    R3[228] = 'NLOS2'
    R3[229] = 'NLOS2'
    R3[230] = 'NLOS2'

    R4 = []

    for i in range(1, 82):
            R4.append('NLOS')

    for i in range(1, 41):
            R4.append('LOS')

    for i in range(1, 12):
            R4.append('NLOS')

    for i in range(1, 6):
            R4.append('NLOS2')
    for i in range(1, 9):
            R4.append('NLOS')

    for i in range(1, 18):
            R4.append('NLOS2')

    for i in range(1, 70):
            R4.append('NLOS')

    for i in range(1, 146):
            R4.append('NLOS2')

    visi = []
    visi.append(R1)
    visi.append(R2)
    visi.append(R3)
    visi.append(R4)
    return visi

def trait(filename, itx=np.array([]), dico={}, h=1):
    """ evaluate various parameters for all measures in itx array

    Parameters
    ----------
    filename
    itx
    dico
    h
        height (default h1)

    Notes
    -----
        F['PL_FS_Rx1'] = -20*log10(((0.3/f)/(4*pi*M.de[0]*0.3)))
        F['PL_FS_Rx2'] = -20*log10(((0.3/f)/(4*pi*M.de[1]*0.3)))
        F['PL_FS_Rx3'] = -20*log10(((0.3/f)/(4*pi*M.de[2]*0.3)))
        F['PL_FS_Rx4'] = -20*log10(((0.3/f)/(4*pi*M.de[3]*0.3)))

    """

    F = {}
    #F['PL']   = {}
    F['Et'] = {}
    F['Es'] = {}
    F['Ef'] = {}
    F['dist'] = {}
    F['tau1'] = {}
    F['taumax'] = {}
    F['taurms'] = {}
    F['taumoy'] = {}
    F['los'] = {}
    F['lqi'] = {}
    for i in itx:
        iTx = dico[i]
        try:
            M = UWBMeasure(iTx, h)
            #f,pl = M.tdd.PL(2,6,0.02)
            #F['PL'] = appendlink(F['PL'],pl)
            etot = M.Etot()
            try:
                F['Tx'] = np.vstack((F['Tx'], M.tx))
            except:
                F['Tx'] = M.tx
            F['Et'] = appendlink(F['Et'], etot)
            emax = M.Emax()
            F['Es'] = appendlink(F['Es'], emax)
            dist = M.de * 0.3
            F['dist'] = appendlink(F['dist'], dist)
            efirst = M.Efirst()
            F['Ef'] = appendlink(F['Ef'], efirst)
            tau1 = M.toa_win()
            F['tau1'] = appendlink(F['tau1'], tau1)
            taurms = M.tau_rms()
            F['taurms'] = appendlink(F['taurms'], taurms)
            taumoy = M.tau_moy()
            F['taumoy'] = appendlink(F['taumoy'], taumoy)
            taumax = M.taumax()
            F['taumax'] = appendlink(F['taumax'], taumax)
            lqi = M.lqi
            F['lqi'] = appendlink(F['lqi'], lqi)
            los = M.type
            vlos = zeros(4)
            if los[0] == 'LOS':
                vlos[0] = 1
            else:
                vlos[0] = 0
            if los[1] == 'LOS':
                vlos[1] = 1
            else:
                vlos[1] = 0
            if los[2] == 'LOS':
                vlos[2] = 1
            else:
                vlos[2] = 0
            if los[3] == 'LOS':
                vlos[3] = 1
            else:
                vlos[3] = 0

            F['los'] = appendlink(F['los'], vlos)
        except:
            pass

    F['Rx'] = M.rx
    io.savemat(filename, F)

    return(F)

class RAW_DATA(PyLayers):
    """

    Members
    -------
    ch1
    ch2
    ch3
    ch4
    time
    timeTX
    tx

    """

    def __init__(self, d):
        #self.time = d.Time[0]
        #self.ch1 = d.CH1[0]
        #self.ch2 = d.CH2[0]
        #self.ch3 = d.CH3[0]
        #self.ch4 = d.CH4[0]
        #self.timetx = d.TimeTX[0]
        #self.tx = d.TX[0]

        self.time = d[0]
        self.ch1 = d[1]
        self.ch2 = d[2]
        self.ch3 = d[3]
        self.ch4 = d[4]
        self.timetx = d[5]
        self.tx = d[6]

class CAL_DATA(PyLayers):
    """
    CAL_DATA

    ch1
    ch2
    ch3
    ch4
    vna_att1
    vna_att2
    vna_freq

    """
    def __init__(self, d):
        """
            depending of the version of scipy
            io.loadma do not give the same output
        """
        #self.ch1 = d.CH1
        #self.ch2 = d.CH2
        #self.ch3 = d.CH3
        #self.ch4 = d.CH4
        #self.vna_freq = d.VNA_Freq
        #self.vna_att1 = d.VNA_Attenuator1
        #self.vna_att2 = d.VNA_Attenuator2

        self.ch1 = d[0]
        self.ch2 = d[1]
        self.ch3 = d[2]
        self.ch4 = d[3]
        self.vna_freq = d[4]
        self.vna_att1 = d[5]
        self.vna_att2 = d[6]

    def plot(self):
        plt.plot(self.ch1)
        plt.show()

    def getwave(self):
        """
        getwave
        """
        #
        # place a window on each channel
        #
        s1 = self.ch1
        s2 = self.ch2

class Fdd(PyLayers):
    """ Frequency Domain Deconv Data

    Attributes
    ----------
    ch1
        channel 1
    ch2
        channel 2
    ch3
        channel 3
    ch4
        channel 4
    freq
        frequency
    tx

    Methods
    -------
        plot

    """
    def __init__(self, d):
        """
        """

        #self.freq = d.Freq[0]*1e-9
        #self.ch1  = d.CH1[0]
        #self.ch2  = d.CH2[0]
        #self.ch3  = d.CH3[0]
        #self.ch4  = d.CH4[0]
        #self.tx   = d.TX[0]

        self.freq = d[0][0] * 1e-9
        self.ch1 = d[1][0]
        self.ch2 = d[2][0]
        self.ch3 = d[3][0]
        self.ch4 = d[4][0]
        self.tx = d[5][0]

    def plot(self, typ='moddB'):
        """ plot Fdd

        Parameters
        ----------

        typ : string
            'moddB' : modulus in dB ,
            'mod', : modulus in linear scale,
            'ang'   : unwraped phase in radians

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.project import *
            >>> from pylayers.measures.mesuwb import *
            >>> import matplotlib.pylab as plt
            >>> M  = UWBMeasure(1)
            >>> F  = M.fdd
            >>> fig = plt.figure()
            >>> F.plot('moddB')
            >>> plt.tight_layout()
            >>> fig = plt.figure()
            >>> F.plot('mod')
            >>> plt.tight_layout()
            >>> fig = plt.figure()
            >>> F.plot('ang')
            >>> plt.tight_layout()
            >>> plt.show()

        """


        f = self.freq
        if (typ == 'moddB'):
            v1 = 20 * np.log10(abs(self.ch1))
            v2 = 20 * np.log10(abs(self.ch2))
            v3 = 20 * np.log10(abs(self.ch3))
            v4 = 20 * np.log10(abs(self.ch4))
        if (typ == 'mod'):
            v1 = abs(self.ch1)
            v2 = abs(self.ch2)
            v3 = abs(self.ch3)
            v4 = abs(self.ch4)
        if (typ == 'ang'):
            v1 = np.unwrap(np.angle(self.ch1))
            v2 = np.unwrap(np.angle(self.ch2))
            v3 = np.unwrap(np.angle(self.ch3))
            v4 = np.unwrap(np.angle(self.ch4))

        plt.subplot(221)
        plt.plot(f, v1)
        plt.xlabel('Freq (GHz)')
        plt.title('CH1')
        plt.subplot(222)
        plt.plot(f, v2)
        plt.xlabel('Freq (GHz)')
        plt.title('CH2')
        plt.subplot(223)
        plt.plot(f, v3)
        plt.xlabel('Freq (GHz)')
        plt.title('CH3')
        plt.subplot(224)
        plt.plot(f, v4)
        plt.xlabel('Freq (GHz)')
        plt.title('CH4')
        plt.show()

class Tdd(PyLayers):
    """
    Time Domain Deconv Data

    Attributes
    ----------

    ch1
        signal on channel 1 (bs.TUsignal)
    ch2
        signal on channel 2 (bs.TUsignal)
    ch3
        signal on channel 3 (bs.TUsignal)
    ch4
        signal on channel 4 (bs.TUsignal)
    tx
        exitation waveform ( Impulse feeding Tx antenna) (bs.TUsignal)

    time is expressed in nano seconds

    Methods
    -------

    PL
        Calculate NB Path Loss on a given UWB frequency band
    box
        return min and max value
    show
        display the 4 channels
    show_span(delay,wide)
    """
    def __init__(self, d=[]):
        if d != []:
            #t = d.Time[0]*1e9
            #self.ch1  = bs.TUsignal(t,d.CH1[0])
            #self.ch2  = bs.TUsignal(t,d.CH2[0])
            #self.ch3  = bs.TUsignal(t,d.CH3[0])
            #self.ch4  = bs.TUsignal(t,d.CH4[0])
            #self.tx   = bs.TUsignal(t,d.TX)
            #self.tx   = bs.TUsignal(t,d.TX[0])
            #####
            ##### Definitely not obvious , but 
            ##### this order is the corect one regarding
            ##### to the measurements !!
            #####
            t = d[0][0] * 1e9
            self.ch3 = bs.TUsignal(t, d[1][0])
            self.ch4 = bs.TUsignal(t, d[2][0])
            self.ch1 = bs.TUsignal(t, d[3][0])
            self.ch2 = bs.TUsignal(t, d[4][0])
            #self.tx   = bs.TUsignal(t,d.TX)
            self.tx = bs.TUsignal(t, d[5][0])
        else:
            pass

    def PL(self, fmin, fmax, B):
        """  Calculate NB Path Loss on a given UWB frequency band

        Parameters
        ----------
        fmin
            start frequency GHz
        fmax
            stop frequency GHz
        B
            NB bandwith (GHz)

        Examples
        --------

       .. plot::
           :include-source:

            >>> from pylayers.util.project import *
            >>> from pylayers.measures.mesuwb import *
            >>> import matplotlib.pylab as plt
            >>> M  = UWBMeasure(1)
            >>> T  = M.tdd
            >>> freq,pl = T.PL(3,7,10)

        """
        S1 = self.ch1.fft()
        S2 = self.ch2.fft()
        S3 = self.ch3.fft()
        S4 = self.ch4.fft()
        df = S4.x[1] - S4.x[0]
        Stx = self.tx.fft()

        freq = np.arange(fmin, fmax, B)
        N = len(freq)
        i_start = np.nonzero(
            (S1.x >= fmin - df / 2.) & (S1.x < fmin + df / 2.))
        i_stop = np.nonzero(
            (S1.x >= fmax - df / 2.) & (S1.x < fmax + df / 2.))
        f = S1.x[i_start[0][0]:i_stop[0][0]]
        for k in range(N):
            u = np.nonzero((f >= fmin + k * B) & (f < fmin + (k + 1) * B))
            Et = 10 * np.log10(np.sum(Stx.y[u[0]] * np.conj(
                Stx.y[u[0]]))).astype('float')
            Er1 = 10 * np.log10(np.sum(
                S1.y[u[0]] * np.conj(S1.y[u[0]]))).astype('float')
            Er2 = 10 * np.log10(np.sum(
                S2.y[u[0]] * np.conj(S2.y[u[0]]))).astype('float')
            Er3 = 10 * np.log10(np.sum(
                S3.y[u[0]] * np.conj(S3.y[u[0]]))).astype('float')
            Er4 = 10 * np.log10(np.sum(
                S4.y[u[0]] * np.conj(S4.y[u[0]]))).astype('float')
            PL1 = Et - Er1
            PL2 = Et - Er2
            PL3 = Et - Er3
            PL4 = Et - Er4
            try:
                pl = np.vstack((pl, np.array([PL1, PL2, PL3, PL4])))
            except:
                pl = np.array([PL1, PL2, PL3, PL4])

        return(freq, pl)

    def show(self, delay=np.array([[0], [0], [0], [0]]), display=True,
             titre=['Rx1', 'Rx2', 'Rx3', 'Rx4'], col=['k', 'b', 'g', 'c'],
             xmin=0, xmax=200, C=0, NC=1,typ='v'):
        """ show the 4 Impulse Radio Impulse responses

        Parameters
        ----------

        delay
            delay values to be displayed vertically [tau0 - toa_th toa_cum toa_max , taum-taurms taum+taurms]
        display
            if display == False the first subplot is reserved for displaying the Layout

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.measures.mesuwb import *
            >>> import matplotlib.pylab as plt
            >>> ntx = 2
            >>> M  = UWBMeasure(ntx)
            >>> T  = M.tdd
            >>> fig = plt.figure()
            >>> t = plt.title('test Tdd.show  Tx='+str(ntx))
            >>> T.show()

        """

        fig = plt.gcf()
        ax = fig.get_axes()

        if (len(ax) < 4) | (NC == 2):
            createfig = True
        else:
            createfig = False

        if display:
            N = 4
            M = 0
        else:
            N = 5
            M = 1
        #plt.axis([0,200,-2,2])
        self.ch1.zlr(xmin, xmax)
        self.ch2.zlr(xmin, xmax)
        self.ch3.zlr(xmin, xmax)
        self.ch4.zlr(xmin, xmax)

        if createfig:
            sp1 = fig.add_subplot(N, NC, NC + C)
        elif N == 4:
            sp1 = ax[0]
        else:
            sp1 = ax[1]

        self.ch1.plot(color=col[0], vline=np.array([delay[0]]),ax=sp1,typ=typ)
        #              showlabel=[False, True], unit1='V', unit2='mV', ax=sp1)
        #self.ch1.plot(col=col[0],unit1='V',unit2='mV',ax=sp1,logx=False,logy=False)
        #plt.show()
        #axis([0,200,-2,2])
        #"sp1.add_title(titre[0])
        if createfig:
            sp2 = fig.add_subplot(N, NC, 2 * NC + C,sharex=sp1)
        elif N == 4:
            sp2 = ax[1]
        else:
            sp2 = ax[2]
        self.ch2.plot(color=col[1], vline=np.array([delay[1]]),ax=sp2,typ=typ)
        #              showlabel=[False, True], unit1='V', unit2='mV', ax=sp2)
        #self.ch2.plot(col=col[1],unit1='V',unit2='mV',ax=sp2,logx=False,logy=False)
        #plt.show()
        #plt.axis([0,200,-2,2])
        #plt.title(titre[1])
        if createfig:
            sp3 = fig.add_subplot(N, NC, 3 * NC + C,sharex=sp1)
        elif N == 4:
            sp3 = ax[2]
        else:
            sp3 = ax[3]
        self.ch3.plot(color=col[2], vline=np.array([delay[2]]),ax=sp3,typ=typ)
        #              showlabel=[False, True], unit1='V', unit2='mV', ax=sp3)
        #self.ch3.plot(col=col[2],unit1='V',unit2='mV',ax=sp3,logx=False,logy=False)
        #plt.show()
        #plt.title(titre[2])
        if createfig:
            sp4 = fig.add_subplot(N, NC, 4 * NC + C,sharex=sp1)
        elif N == 4:
            sp4 = ax[3]
        else:
            sp4 = ax[4]
        self.ch4.plot(color=col[3], vline=np.array([delay[3]]), ax=sp4,typ=typ)
        #              showlabel=[True, True], unit1='V', unit2='mV', ax=sp4)
        #self.ch4.plot(col=col[3],unit1='V',unit2='mV',ax=sp4,logx=False,logy=False)
        #plt.show()
        #plt.axis([0,200,-2,2])
        #plt.title(titre[3])
        if display:
            plt.show()

    def show_span(self, delay=np.array([0, 0, 0, 0]), wide=np.array([0, 0, 0, 0])):
        """ show span

        Examples
        --------

        .. plot::
            :include-source:

            >>> from pylayers.util.project import *
            >>> from pylayers.measures.mesuwb import *
            >>> M  = UWBMeasure(2)
            >>> T  = M.tdd
            >>> s1 = T.show_span()
            >>> plt.show()

        """
        fig = plt.figure()
        sp1 = fig.add_subplot(221)
        plt.plot(self.ch1.x, self.ch1.y)
        plt.axvspan(delay[0] - wide[0] / 2, delay[0] + wide[0] / 2,
                    facecolor='g', alpha=0.5)
        plt.xlabel('Time (ns)')
        plt.title('Rx1')
        sp2 = fig.add_subplot(222)
        plt.plot(self.ch2.x, self.ch2.y)
        plt.axvspan(delay[1] - wide[1] / 2, delay[1] + wide[1] / 2,
                    facecolor='g', alpha=0.5)
        plt.xlabel('Time (ns)')
        plt.title('Rx2')
        plt.sp3 = fig.add_subplot(223)
        plt.plot(self.ch3.x, self.ch3.y)
        plt.axvspan(delay[2] - wide[2] / 2, delay[2] + wide[2] / 2,
                    facecolor='g', alpha=0.5)
        plt.xlabel('Time (ns)')
        plt.title('Rx3')
        sp4 = fig.add_subplot(224)
        plt.plot(self.ch4.x, self.ch4.y)
        plt.axvspan(delay[3] - wide[3] / 2, delay[3] + wide[3] / 2,
                    facecolor='g', alpha=0.5)
        plt.xlabel('Time (ns)')
        plt.title('Rx4')

        plt.show()
        return(fig)

    def box(self):
        """ evaluate min and max of the 4 channels

        Examples
        --------

        >>> from pylayers.measures.mesuwb import *
        >>> M  = UWBMeasure(1)
        >>> T  = M.tdd
        >>> bo = T.box()

        """
        max1 = max(self.ch1.y)
        min1 = min(self.ch1.y)

        max2 = max(self.ch2.y)
        min2 = min(self.ch2.y)

        max3 = max(self.ch3.y)
        min3 = min(self.ch3.y)

        max4 = max(self.ch4.y)
        min4 = min(self.ch4.y)

        ma = max(max1, max2, max3, max4)
        mi = min(min1, min2, min3, min4)

        b = np.array([mi, ma])
        return(b)

    def plot(self, type='raw'):
        """
        type : raw | filter
        """

        b = self.box()
        t = self.time * 1e9

        tmin = min(t)
        tmax = max(t)
        ymin = b[0]
        ymax = b[1]
        ax = [tmin, tmax, ymin, ymax]

        s1 = self.ch1
        s2 = self.ch2
        s3 = self.ch3
        s4 = self.ch4
        if type == 'filter':
            [b, a] = signal.butter(5, 0.2)
            [w, h] = signal.freqz(b, a)
            plt.plot(w, abs(h))
            plt.show()
            s1 = signal.lfilter(b, a, s1)
            s2 = signal.lfilter(b, a, s2)
            s3 = signal.lfilter(b, a, s3)
            s4 = signal.lfilter(b, a, s4)

        plt.subplot(221)
        plt.plot(t, s1)
        plt.xlabel('Time (ns)')
        plt.title('CH1')
        plt.axis(ax)
        plt.subplot(222)
        plt.plot(t, s2)
        plt.xlabel('Time (ns)')
        plt.title('CH2')
        plt.axis(ax)
        plt.subplot(223)
        plt.plot(t, s3)
        plt.xlabel('Time (ns)')
        plt.title('CH3')
        plt.axis(ax)
        plt.subplot(224)
        plt.plot(t, s4)
        plt.xlabel('Time (ns)')
        plt.title('CH4')
        plt.axis(ax)
        plt. show()

class TFP(object):
    """
    Tx
    Rx
    distance
    los_cond

    Etot
    Emax
    Etau0
    tau_moy
    tau_rms
    Epercent
    toa_max
    toa_th
    toa_cum
    angular
    """

    def __init__(self):
        self.metadata = {}
        self.metadata['Tx'] = np.array([]).reshape(3, 0)
        self.metadata['Rx'] = np.array([]).reshape(3, 0)
        self.metadata['distance'] = np.array([])
        self.metadata['los_cond'] = np.parray([])

        self.chan_param = {}
        self.chan_param['Etot'] = np.array([])
        self.chan_param['Emax'] = np.array([])
        self.chan_param['Etau0'] = np.array([])
        self.chan_param['tau_moy'] = np.array([])
        self.chan_param['tau_rms'] = np.array([])
        self.chan_param['Epercent'] = np.array([])
        self.chan_param['toa_max'] = np.array([])
        self.chan_param['toa_th'] = np.array([])
        self.chan_param['toa_cum'] = np.array([])
        self.chan_param['toa_new'] = np.array([])
        self.chan_param['toa_win'] = np.array([])
        self.chan_param['toa_seu'] = np.array([])
        self.chan_param['angular'] = np.array([])

    def append(self, FP):
        """
        Parameters
        ----------
        FP 
        """
        tx = FP.metadata['Tx'].reshape(3, 1)
        rx = FP.metadata['Rx'].reshape(3, 1)
        self.metadata['Tx'] = hstack((self.metadata['Tx'], tx))
        self.metadata['Rx'] = hstack((self.metadata['Rx'], rx))
        self.metadata['distance'] = hstack(
            (self.metadata['distance'], FP.metadata['distance']))
        self.metadata['delay'] = hstack(
            (self.metadata['delay'], FP.metadata['delay']))
        self.metadata['los_cond'] = hstack(
            (self.metadata['los_cond'], FP.metadata['los_cond']))
        self.metadata['angular'] = hstack(
            (self.metadata['angular'], FP.metadata['angular']))

        self.chan_param['Etot'] = hstack(
            (self.chan_param['Etot'], FP.chan_param['Etot']))
        self.chan_param['Emax'] = hstack(
            (self.chan_param['Emax'], FP.chan_param['Emax']))
        self.chan_param['Etau0'] = hstack(
            (self.chan_param['Etau0'], FP.chan_param['Etau0']))
        self.chan_param['tau_moy'] = hstack(
            (self.chan_param['tau_moy'], FP.chan_param['tau_moy']))
        self.chan_param['tau_rms'] = hstack(
            (self.chan_param['tau_rms'], FP.chan_param['tau_rms']))
        self.chan_param['Epercent'] = hstack(
            (self.chan_param['Epercent'], FP.chan_param['Epercent']))
        self.chan_param['toa_max'] = hstack(
            (self.chan_param['toa_max'], FP.chan_param['toa_max']))
        self.chan_param['toa_th'] = hstack(
            (self.chan_param['toa_th'], FP.chan_param['toa_th']))
        self.chan_param['toa_cum'] = hstack(
            (self.chan_param['toa_cum'], FP.chan_param['toa_cum']))
        #self.chan_param['toa_cum_tmtm']=hstack((self.chan_param['toa_cum_tmtm'],FP.chan_param['toa_cum_tmtm']))
        #self.chan_param['toa_cum_tm']=hstack((self.chan_param['toa_cum_tm'],FP.chan_param['toa_cum_tm']))
        #self.chan_param['toa_cum_tmt']=hstack((self.chan_param['toa_cum_tmt'],FP.chan_param['toa_cum_tmt']))
        self.chan_param['toa_win'] = hstack(
            (self.chan_param['toa_win'], FP.chan_param['toa_win']))
        self.chan_param['toa_new'] = hstack(
            (self.chan_param['toa_new'], FP.chan_param['toa_new']))
        #self.chan_param['toa_th_tmtm']=hstack((self.chan_param['toa_th_tmtm'],FP.chan_param['toa_th_tmtm']))
        #self.chan_param['toa_th_tm']=hstack((self.chan_param['toa_th_tm'],FP.chan_param['toa_th_tm']))
        #self .chan_param['toa_th_tmt']=hstack((self.chan_param['toa_th_tmt'],FP.chan_param['toa_th_tmt']))

class FP(object):
    """
    Fingerprint

    F = FP(M,Rxid, alpha = 0.1,Tint=0.225,sym=0.25,nint=7,thlos=0.05,
        thnlos=0.65,thcum=0.15)
    Rxid : 1 | 2 | 3 | 4
    alpha : pdf percentile to suppress below and above (used in tau_rms and tau_moy)
    Tint : Integration time for Emax (ns)
    Tsym : Integration symetry default (0.25 | 0,75)
    nint : number of interval for toa_max
    thlos : threshold for los case  toa_th
    thnlos : threshold for nlos case toa_th
    thcum  : threshold for toa_cum

    This method build the Fingerprint from UWBMeasure M - Rx number k

    Attributes
    ----------
        metadata :  dictionnary
                    Tx
                    Rx
                    distance (m)
                    type
        chan_param : dictionnary
                    Etot
                    Emax
                    Etau0
                    tau_moy
                    tau_rms
                    Epercent
                    toa_max
                    toa_th
                    toa_cum
    """
    def __init__(self, M, k, alpha=0.1, Tint=0.225, sym=0.25, nint=8, thlos=0.05, thnlos=0.15, thcum=0.15, w=6):
        """

        Parameters
        ----------
            M : UWB Mesure
            k : Rx index
            alpha : quantile parameter
            Tint  : Integation time (ns)
            sym   : Symmetry
            nint  : int
            thlos : float
            w     : w
        """
        #
        # Metadata : general parameters links to measurement configuration
        #
        self.metadata = {}
        self.metadata['Tx'] = M.tx
        self.metadata['Rx'] = M.rx[k, :]
        self.metadata['distance'] = M.de[k - 1] * 0.3
        self.metadata['delay'] = M.de[k - 1]
        self.metadata['los_cond'] = M.type[k - 1]
        self.metadata['angular'] = M.ag[k - 1]
        #
        # Chan_param : radio channel features
        #
        self.chan_param = {}
        if k == 1:
            Etot = M.tdd.ch1.Etot(tau0=M.de[0], dB=True)
            Emax = M.tdd.ch1.Emax(Tint, sym, dB=True)
            #Etot = M.tdd.ch1.Etot(tau0=M.de[0],dB=False)
            #Emax = M.tdd.ch1.Emax(Tint,sym,dB=False)
            Etau0 = M.tdd.ch1.Etau0(tau0=M.de[0])
            tau_moy = M.tdd.ch1.tau_moy(alpha)
            tau_rms = M.tdd.ch1.tau_rms(alpha)
            Epercent = M.tdd.ch1.Epercent()
            toa_max = M.tdd.ch1.toa_max(nint)
            toa_th = M.tdd.ch1.toa_th(thlos, thnlos, visibility=M.type[0])
            toa_cum = M.tdd.ch1.toa_cum(thcum)
            #toa_cum_tmtm  = M.tdd.ch1.toa_cum_tmtm()
            #toa_cum_tm  = M.tdd.ch1.toa_cum_tm()
            #toa_cum_tmt  = M.tdd.ch1.toa_cum_tmt()
            #toa_new  = M.tdd.ch1.toa_new()
            toa_win = M.tdd.ch1.toa_win(w)
            #toa_th_tmtm = M.tdd.ch1.toa_th_tmtm()
            #toa_th_tm = M.tdd.ch1.toa_th_tm()
            #toa_th_tmt = M.tdd.ch1.toa_th_tmt()
        if k == 2:
            Etot = M.tdd.ch2.Etot(tau0=M.de[1], dB=True)
            Emax = M.tdd.ch2.Emax(Tint, sym, dB=True)
            #Etot = M.tdd.ch2.Etot(tau0=M.de[1],dB=False)
            #Emax = M.tdd.ch2.Emax(Tint,sym,dB=False)
            Etau0 = M.tdd.ch2.Etau0(tau0=M.de[1])
            tau_moy = M.tdd.ch2.tau_moy(alpha)
            tau_rms = M.tdd.ch2.tau_rms(alpha)
            Epercent = M.tdd.ch2.Epercent()
            toa_max = M.tdd.ch2.toa_max(nint)
            toa_th = M.tdd.ch2.toa_th(thlos, thnlos, visibility=M.type[1])
            toa_cum = M.tdd.ch2.toa_cum(thcum)
            #toa_cum_tmtm  = M.tdd.ch2.toa_cum_tmtm()
            #toa_cum_tm  = M.tdd.ch2.toa_cum_tm()
            #toa_cum_tmt  = M.tdd.ch2.toa_cum_tmt()
            #toa_new  = M.tdd.ch2.toa_new()
            toa_win = M.tdd.ch2.toa_win(w)
            #toa_th_tmtm = M.tdd.ch2.toa_th_tmtm()
            #toa_th_tm = M.tdd.ch2.toa_th_tm()
            #toa_th_tmt = M.tdd.ch2.toa_th_tmt()
        if k == 3:
            Etot = M.tdd.ch3.Etot(tau0=M.de[2], dB=True)
            Emax = M.tdd.ch3.Emax(Tint, sym, dB=True)
            #Etot = M.tdd.ch3.Etot(tau0=M.de[2],dB=False)
            #Emax = M.tdd.ch3.Emax(Tint,sym,dB=False)
            Etau0 = M.tdd.ch3.Etau0(tau0=M.de[2])
            tau_moy = M.tdd.ch3.tau_moy(alpha)
            tau_rms = M.tdd.ch3.tau_rms(alpha)
            Epercent = M.tdd.ch3.Epercent()
            toa_max = M.tdd.ch3.toa_max(nint)
            toa_th = M.tdd.ch3.toa_th(thlos, thnlos, visibility=M.type[2])
            toa_cum = M.tdd.ch3.toa_cum(thcum)
            #toa_cum_tmtm  = M.tdd.ch3.toa_cum_tmtm()
            #toa_cum_tm  = M.tdd.ch3.toa_cum_tm()
            #toa_cum_tmt  = M.tdd.ch3.toa_cum_tmt()
            #toa_new  = M.tdd.ch3.toa_new()
            toa_win = M.tdd.ch3.toa_win(w)
            #toa_th_tmtm = M.tdd.ch3.toa_th_tmtm()
            #toa_th_tm = M.tdd.ch3.toa_th_tm()
            #toa_th_tmt = M.tdd.ch3.toa_th_tmt()
        if k == 4:
            Etot = M.tdd.ch4.Etot(tau0=M.de[3], dB=True)
            Emax = M.tdd.ch4.Emax(Tint, sym, dB=True)
            #Etot = M.tdd.ch4.Etot(tau0=M.de[3],dB=False)
            #Emax = M.tdd.ch4.Emax(Tint,sym,dB=False)
            Etau0 = M.tdd.ch4.Etau0(tau0=M.de[3])
            tau_moy = M.tdd.ch4.tau_moy(alpha)
            tau_rms = M.tdd.ch4.tau_rms(alpha)
            Epercent = M.tdd.ch4.Epercent()
            toa_max = M.tdd.ch4.toa_max(nint)
            toa_th = M.tdd.ch4.toa_th(thlos, thnlos, visibility=M.type[3])
            toa_cum = M.tdd.ch4.toa_cum(thcum)
            #toa_cum_tmtm  = M.tdd.ch4.toa_cum_tmtm()
            #toa_cum_tm  = M.tdd.ch4.toa_cum_tm()
            #toa_cum_tmt  = M.tdd.ch4.toa_cum_tmt()
            #toa_new  = M.tdd.ch4.toa_new()
            toa_win = M.tdd.ch4.toa_win(w)
            #toa_th_tmtm = M.tdd.ch4.toa_th_tmtm()
            #toa_th_tm = M.tdd.ch4.toa_th_tm()
            #toa_th_tmt = M.tdd.ch4.toa_th_tmt()

        self.chan_param['Etot'] = Etot
        self.chan_param['Emax'] = Emax
        self.chan_param['Etau0'] = Etau0
        self.chan_param['tau_moy'] = tau_moy
        self.chan_param['tau_rms'] = tau_rms
        self.chan_param['Epercent'] = Epercent
        self.chan_param['toa_max'] = toa_max
        self.chan_param['toa_th'] = toa_th
        self.chan_param['toa_cum'] = toa_cum
        #self.chan_param['toa_cum_tmtm']=toa_cum_tmtm
        #self.chan_param['toa_cum_tm']=toa_cum_tm
        #self.chan_param['toa_cum_tmt']=toa_cum_tmt
        #self.chan_param['toa_new']=toa_new
        self.chan_param['toa_win'] = toa_win
        #self.chan_param['toa_th_tmtm']=toa_th_tmtm
        #self.chan_param['toa_th_tm']=toa_th_tm
        #self. chan_param['toa_th_tmt']=toa_th_tmt
#\hline
#\hline
#$E_{tot}$ & ~& & & \\
#\hline
# & & & & \\
# \hline
# & & & & \\
# \hline
# & & & &
#\\
#\hline
# & & & &
#\\
#\hline
# & & & &
#\\
#\hline
# & & & &
#\\
#\hline
# & & & &
#\\
#\hline
#\end{supertabular}
#\end{center}
#\bigskip


class UWBMeasure(PyLayers):
    """ UWBMeasure class

    Attributes
    ----------

    tdd
        Time domain deconv data
    fdd
        Freq domain deconv data
    Date_Time
    LQI
    Operators
    RAW_DATA
    CAL_DATAip
    Tx_height
    Tx_position

    Methods
    -------

    info()  :
    show()
    emax()
    etot()
    TDoA()
    Fingerprint()
    fp() : calculate fingerprint
    TOA()

    Notes
    -----

    """

    def __init__(self, nTx=1, h=1, display=False):
        """

        Parameters
        ----------

        nTx : int
            Tx index
        h : int
            height indicator (default 1 - 1.2m)  0 - 1.5 m

         Examples
         --------

         .. plot::
            :include-source:

            >>> from pylayers.measures.mesuwb import *
            >>> M1 = UWBMeasure(1)
            >>> M1.show()

        """
        super(UWBMeasure,self).__init__()

        self.validindex = [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
        23,24,25,26,27,28,29,30,32,33,34,35,36,37,38,39,40,41,42,
        43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,
        62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
        81,82,83,84,85,89,90,91,92,93,94,95,96,97,98,99,100,101,103,
        104,105,106,107,108,109,110,111,113,114,116,117,119,120,122,
        123,124,125,126,127,128,129,133,134,136,137,138,139,140,141,
        142,143,144,145,146,147,162,163,164,165,166,167,168,169,170,
        171,172,173,174,175,176,177,179,180,181,182,183,184,185,186,
        188,189,199,200,201,202,203,204,205,206,207,208,209,210,211,
        212,213,214,215,216,217,218,219,220,221,222,223,227,228,229,
        230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,
        245,246,247,248,249,250,251,252,253,258,259,266,267,268,269,
        270,271,272,273,274,275,276,277,278,279,297,298,299,300,301,
        302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,
        317,318,319,320,321,322,323,324,325,326,327,328,329,330,332,
        333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,
        348,349,350,351,352,353,354,355,356,360,361,362,363,364,365,
        366,367,368,369,370,371,372,373,374,375]

        # Raw data matlab file reading
        Tx, Rx = ptw1()
        visi = visibility()
        if not os.path.exists(mesdir):
            raise AttributeError('Incorrect Measure directory set in $MESDIR')
        filename = mesname(nTx, mesdir, h)
        if os.path.exists(filename) :
            b = io.loadmat(filename)
            # Conversion
            self.CAL_DATA = CAL_DATA(b['CAL_DATA'][0][0])
            self.RAW_DATA = RAW_DATA(b['RAW_DATA'][0][0])
            #T              = b['DECONV_DATA'][0][0].TIME_DOMAIN[0][0]
            #F              = b['DECONV_DATA'][0][0].FREQ_DOMAIN[0][0]
            #self.tdd       = Tdd(T)
            #self.fdd       = Fdd(F)
            self.tdd = Tdd(b['DECONV_DATA'][0][0][0][0][0])
            self.fdd = Fdd(b['DECONV_DATA'][0][0][1][0][0])
            self.Date_Time = b['Date_Time']
            self.LQI = {}
            #self.LQI['Method1_CH1']=b['LQI'][0][0].Method1_CH1[0][0]
            #self.LQI['Method1_CH2']=b['LQI'][0][0].Method1_CH2[0][0]
            #self.LQI['Method1_CH3']=b['LQI'][0][0].Method1_CH3[0][0]
            #self.LQI['Method1_CH4']=b['LQI'][0][0].Method1_CH4[0][0]
            #self.LQI['Method2_CH1']=b['LQI'][0][0].Method2_CH1[0][0]
            #self.LQI['Method2_CH2']=b['LQI'][0][0].Method2_CH2[0][0]
            #self.LQI['Method2_CH3']=b['LQI'][0][0].Method2_CH3[0][0]
            #self.LQI['Method2_CH4']=b['LQI'][0][0].Method2_CH4[0][0]

            self.LQI['Method1_CH1'] = b['LQI'][0][0][0][0][0]
            self.LQI['Method1_CH2'] = b['LQI'][0][0][1][0][0]
            self.LQI['Method1_CH3'] = b['LQI'][0][0][2][0][0]
            self.LQI['Method1_CH4'] = b['LQI'][0][0][3][0][0]
            self.LQI['Method2_CH1'] = b['LQI'][0][0][4][0][0]
            self.LQI['Method2_CH2'] = b['LQI'][0][0][5][0][0]
            self.LQI['Method2_CH3'] = b['LQI'][0][0][6][0][0]
            self.LQI['Method2_CH4'] = b['LQI'][0][0][7][0][0]

            lqi1 = self.LQI['Method1_CH1']
            lqi2 = self.LQI['Method1_CH2']
            lqi3 = self.LQI['Method1_CH3']
            lqi4 = self.LQI['Method1_CH4']
            self.lqi = np.array([lqi1, lqi2, lqi3, lqi4])

            self.Operators = b['Operators']
            self.Tx_height = b['Tx_height']
            self.Tx_position = b['Tx_position']
            # appending geometrical information
            if self.Tx_height == '120cm':
                d = 1.2
            else:
                d = 1.5
            self.rx = Rx
            self.tx = np.hstack((Tx[nTx, :], d))
            self.ntx = nTx
            # calculate delay
            d1 = pyu.delay(self.rx[1, :], self.tx)
            d2 = pyu.delay(self.rx[2, :], self.tx)
            d3 = pyu.delay(self.rx[3, :], self.tx)
            d4 = pyu.delay(self.rx[4, :], self.tx)
            self.de = np.array([d1, d2, d3, d4])
            #display measurements
            if display:
                self.tdd.show(self.de)
            #calculate angle
            ang1 = geu.angular(self.tx, self.rx[1, :])
            ang2 = geu.angular(self.tx, self.rx[2, :])
            ang3 = geu.angular(self.tx, self.rx[3, :])
            ang4 = geu.angular(self.tx, self.rx[4, :])
            self.ag = np.array([ang1, ang2, ang3, ang4])

            type1 = visi[0][nTx - 1]
            type2 = visi[1][nTx - 1]
            type3 = visi[2][nTx - 1]
            type4 = visi[3][nTx - 1]
            self.type = [type1, type2, type3, type4]
            self.valid = True
        else:
            raise AttributeError("non valid Rx point ")
            self.valid = False

#    def __del__(self):
#        if self.valid:
#            del self.F1
#            del self.F2
#            del self.F3
#            del self.F4
#            print "Detruit: ",self.Tx_position




         


    def info(self):
        print "Date_Time :", self.Date_Time
        print "Operators : ", self.Operators
        print "Tx_height :", self.Tx_height
        print "Tx_position :", self.Tx_position
        print "Tx : ", self.tx
        for i in range(4):
            print "------Tx" + str(i + 1) + " ------"
            print "delays     (ns):", self.de[i]
            print "range  (meters):", self.de[i] * 0.3
            print "visibility     :", self.type[i]
            print "angular (degree)  :", self.ag[i]
            if i == 0:
                print "LQI Meth1", self.LQI['Method1_CH1'], " (dB)"
                print "LQI Meth2", self.LQI['Method2_CH1'], " (dB)"
            if i == 1:
                print "LQI Meth1", self.LQI['Method1_CH2'], " (dB)"
                print "LQI Meth2", self.LQI['Method2_CH2'], " (dB)"
            if i == 2:
                print "LQI Meth1", self.LQI['Method1_CH3'], " (dB)"
                print "LQI Meth2", self.LQI['Method2_CH3'], " (dB)"
            if i == 3:
                print "LQI Meth1", self.LQI['Method1_CH4'], " (dB)"
                print "LQI Meth2", self.LQI['Method2_CH4'], " (dB)"

    #def pass_loss(self):
        #"""
        #calculate the path Loss for each subband
        #"""

        #y1 = self.fdd.ch1
        #y2 = self.fdd.ch2
        #y3 = self.fdd.ch3
        #y4 = self.fdd.ch4
        #f  = self.fdd.freq
        #yt = self.fdd.tx

        #PL1 = np.array([])
        #PL2 = np.array([])
        #PL3 = np.array([])
        #PL4 = np.array([])
        #delta = 0.02
        #df = 0.005

        #s1 = 20*log10(abs(y1))
        #s2 = 20*log10(abs(y2))
        #s3 = 20*log10(abs(y3))
        #s4 = 20*log10(abs(y4))
        #st = 20*log10(abs(yt))

        #for i in range (399):

            #es1 = s1[((delta*(i+1)-delta/2)/df):((delta*(i+1)+delta/2)/df)]
            #es2 = s2[((delta*(i+1)-delta/2)/df):((delta*(i+1)+delta/2)/df)]
            #es3 = s3[((delta*(i+1)-delta/2)/df):((delta*(i+1)+delta/2)/df)]
            #es4 = s4[((delta*(i+1)-delta/2)/df):((delta*(i+1)+delta/2)/df)]
            #est = st[((delta*(i+1)-delta/2)/df):((delta*(i+1)+delta/2)/df)]

            #cumS1 = cumsum(es1)
            #cumS2 = cumsum(es2)
            #cumS3 = cumsum(es3)
            #cumS4 = cumsum(es4)
            #cumSt = cumsum(est)

            #S1 = cumS1[-1]
            #S2 = cumS2[-1]
            #S3 = cumS3[-1]
            #S4 = cumS4[-1]
            #St = cumSt[-1]

            #Pr1 = S1*df
            #Pr2 = S2*df
            #Pr3 = S3*df
            #Pr4 = S4*df
            #Pt  = St*df

            #pl1  = Pt-Pr1
            #pl2  = Pt-Pr2
            #pl3  = Pt-Pr3
            #pl4  = Pt-Pr4

            #PL1 = hstack((PL1,pl1))
            #PL2 = hstack((PL1,pl2))
            #PL3 = hstack((PL1,pl3))
            #PL4 = hstack((PL1,pl4))

        #return PL1,PL2,PL3,PL4
    def show(self, delay=np.array([[0], [0], [0], [0]]), display=True,
             col=['k', 'b', 'g', 'c'], xmin=0, xmax=100, C=0, NC=1,typ='v'):
        """ show measurement in time domain

        Parameters
        ----------
        delay : np.array(1,4)
        display
            optional
        col
            optional
        xmin
            optional
        xmax
            optional
        C
            optional
         NC
            optional

        """
        titre = []
        self.tdd.show(delay,
                      display=display,
                      titre=['Rx1 ' + 'Tx' + str(self.ntx),
                             'Rx2', 'Rx3', 'Rx4'],
                      col=col,
                      xmin=xmin,
                      xmax=xmax,
                      C=C,
                      NC=NC,
                     typ='v')

    def Epercent(self):
        """
        """

        epercent1 = self.tdd.ch1.Epercent()
        epercent2 = self.tdd.ch2.Epercent()
        epercent3 = self.tdd.ch3.Epercent()
        epercent4 = self.tdd.ch4.Epercent()
        epercent = np.array([epercent1, epercent2, epercent3, epercent4])

        return epercent

    #def TDoA(self):
        #"""
        #"""

        #y1 = self.tdd.ch1.TDoA()
        #y2 = self.tdd.ch2.TDoA()
        #y3 = self.tdd.ch3.TDoA()
        #y4 = self.tdd.ch4.TDoA()

                #c12 = correlate(abs(y1),abs(y2),mode='full')
                #c13 = correlate(abs(y1),abs(y3),mode='full')
                #c14 = correlate(abs(y1),abs(y4),mode='full')
                #n12 = find(c12==max(c12))
                ##n12 = len(y1)-1-n12
                #t12 = 0.005*(n12-int(len(c12)/2))
                ##tdoa12 = abs((self.de[1]-self.de[0])/0.3)
                ##err12 = t12-tdoa12
                #n13 = find(c13==max(c13))
                ##n13 = len(y1)-1-n13
                #t13 = 0.005*(n13-int(len(c13)/2))
                ##tdoa13 = abs((self.de[2]-self.de[0])/0.3)
                ##err13 = t13-tdoa13
                #n14 = find(c14==max(c14))
                ##n14 = len(y1)-1-n14
                #t14 = 0.005*(n14-int(len(c14)/2))
                ##tdoa14 = abs((self.de[3]-self.de[0])/0.3)
                ##err14 = t14-tdoa14
                #tdoa = np.array([t12,t13,t14])
                #return tdoa
    def toa_max2(self):
        """ calculate toa_max (meth2)
        """
        toa_max21 = self.tdd.ch1.toa_max2()
        toa_max22 = self.tdd.ch2.toa_max2()
        toa_max23 = self.tdd.ch3.toa_max2()
        toa_max24 = self.tdd.ch4.toa_max2()

    def tau_Emax(self):
        """
        Calculate the delay of energy peak
        """
        tau_Emax1 = self.tdd.ch1.tau_Emax()
        tau_Emax2 = self.tdd.ch2.tau_Emax()
        tau_Emax3 = self.tdd.ch3.tau_Emax()
        tau_Emax4 = self.tdd.ch4.tau_Emax()
        tau_Emax = np.array([tau_Emax1, tau_Emax2,
                             tau_Emax3, tau_Emax4])

        return tau_Emax

    def tau_moy(self, display=False):
        """
            calculate mean excess delay
        """

        taum1 = self.tdd.ch1.tau_moy()
        taum2 = self.tdd.ch2.tau_moy()
        taum3 = self.tdd.ch3.tau_moy()
        taum4 = self.tdd.ch4.tau_moy()

        taum = np.array([taum1, taum2, taum3, taum4])
        if display:
            self.tdd.show(taum)
        return taum

    def tau_rms(self, display=False):
        """ calculate the RMS delay spread
        """

        taurms1 = self.tdd.ch1.tau_rms()
        taurms2 = self.tdd.ch2.tau_rms()
        taurms3 = self.tdd.ch3.tau_rms()
        taurms4 = self.tdd.ch4.tau_rms()

        taurms = np.array([taurms1, taurms2, taurms3, taurms4])
        if display:
            self.tdd.show_span(self.tau_moy(display=False), taurms)
        return taurms

    def toa_new(self, display=False):
        """ descendant threshold based toa estimation
        """
        toa1 = self.tdd.ch1.toa_new()
        toa2 = self.tdd.ch2.toa_new()
        toa3 = self.tdd.ch3.toa_new()
        toa4 = self.tdd.ch4.toa_new()

        toa = np.array([toa1, toa2, toa3, toa4])
        if display:
            self.tdd.show(toa)
        return toa

    def toa_win(self, n=9, display=False):
        """ descendant threshold based toa estimation

        Parameters
        ----------
        n  : key parameter n = 9
        display : False

        """
        toa1 = self.tdd.ch1.toa_win(w=n)
        toa2 = self.tdd.ch2.toa_win(w=n)
        toa3 = self.tdd.ch3.toa_win(w=n)
        toa4 = self.tdd.ch4.toa_win(w=n)

        toa = np.array([toa1, toa2, toa3, toa4])
        if display:
            self.tdd.show(toa)
        return toa

    def toa_max(self, n, display=False):
        """
            descendant threshold based toa estimation
        """

        toa1 = self.tdd.ch1.toa_max(nint=n)
        toa2 = self.tdd.ch2.toa_max(nint=n)
        toa3 = self.tdd.ch3.toa_max(nint=n)
        toa4 = self.tdd.ch4.toa_max(nint=n)

        toa = np.array([toa1, toa2, toa3, toa4])
        if display:
            self.tdd.show(toa)
        return toa

    def toa_th(self, r, k, display=False):
        """
            threshold based toa estimation using energy pic
        """

        toa1 = self.tdd.ch1.toa_th(visibility=self.type[0], thlos=r, thnlos=k)
        toa2 = self.tdd.ch2.toa_th(visibility=self.type[1], thlos=r, thnlos=k)
        toa3 = self.tdd.ch3.toa_th(visibility=self.type[2], thlos=r, thnlos=k)
        toa4 = self.tdd.ch4.toa_th(visibility=self.type[3], thlos=r, thnlos=k)

        toa = np.array([toa1, toa2, toa3, toa4])

        if display:
            self.tdd.show(toa)
        return toa

    #def toa_th(self,display=False):

        #"""
        #threshold based toa estimation using energy pic
        #"""

        #toa1  = self.tdd.ch1.toa_th_tmtm()
        #toa2  = self.tdd.ch2.toa_th_tmtm()
        #toa3  = self.tdd.ch3.toa_th_tmtm()
         #toa4  = self.tdd.ch4.toa_th_tmtm()

        #toa   = np.array([toa1,toa2,toa3,toa4])
                #if display ==True:
            #self.tdd.show(toa)
        #return toa

    def toa_cum(self, n, display=False):
        """ threshold based toa estimation using cumulative energy
        """

        toa1 = self.tdd.ch1.toa_cum(th=n)
        toa2 = self.tdd.ch2.toa_cum(th=n)
        toa3 = self.tdd.ch3.toa_cum(th=n)
        toa4 = self.tdd.ch4.toa_cum(th=n)

        toa = np.array([toa1, toa2, toa3, toa4])
        if display:
            self.tdd.show(toa)
        return toa

    #def toa_cum(self,display=False):
        #"""
        #threshold based toa estimation using cumulative energy
        #"""

        #toa1 = self.tdd.ch1.toa_cum_tmt()
        #toa2 = self.tdd.ch2.toa_cum_tmt()
        #toa3 = self.tdd.ch3.toa_cum_tmt()
        #toa4 = self.tdd.ch4.toa_cum_tmt()

        #toa  =  np.array([toa1,toa2,toa3,toa4])
        #if display:
            #self.tdd.show(toa)
        #return toa

    def taumax(self):
        """

        """
        tmax1 = self.tdd.ch1.taumax()
        tmax2 = self.tdd.ch2.taumax()
        tmax3 = self.tdd.ch3.taumax()
        tmax4 = self.tdd.ch4.taumax()
        taumx = np.array([tmax1, tmax2, tmax3, tmax4])
        return  taumx

    def Emax(self, Tint=1, sym=0.25, dB=True):
        """
            Emax(Tint=0.225,sym=0.25)
            calculate maximum energy
            Time integration 0.225 ns , sym=0.25
        """
        taumax = self.taumax()
        Emax1 = self.tdd.ch1.Ewin(taumax[0], Tint=Tint, sym=sym, dB=dB)
        Emax2 = self.tdd.ch2.Ewin(taumax[1], Tint=Tint, sym=sym, dB=dB)
        Emax3 = self.tdd.ch3.Ewin(taumax[2], Tint=Tint, sym=sym, dB=dB)
        Emax4 = self.tdd.ch4.Ewin(taumax[3], Tint=Tint, sym=sym, dB=dB)

        emax = np.array([Emax1, Emax2, Emax3, Emax4])
        return emax

    def Etot(self):
        """
            Etot

            Calculate total energy for the 4 channels
        """
        de0 = self.de[0] + 0.7
        de1 = self.de[1] + 0.7
        de2 = self.de[2] + 0.7
        de3 = self.de[3] + 0.7
        Etot1 = self.tdd.ch1.Etot(de0, de0 + 75)
        Etot2 = self.tdd.ch2.Etot(de1, de1 + 75)
        Etot3 = self.tdd.ch3.Etot(de2, de2 + 75)
        Etot4 = self.tdd.ch4.Etot(de3, de3 + 75)

        etot = np.array([10 * log10(Etot1), 10 * log10(Etot2), 10 *
                         log10(Etot3), 10 * log10(Etot4)])

        return etot

    def Efirst(self, Tint=1, sym=0.25, dB=True):
        """
        """
        # ???
        #sig = 1/(2*np.sqrt(22))

        #phi1 = self.tdd.ch1.correlate(self.tdd.tx,True)
        #phi2 = self.tdd.ch2.correlate(self.tdd.tx,True)
        #phi3 = self.tdd.ch3.correlate(self.tdd.tx,True)
        #phi4 = self.tdd.ch4.correlate(self.tdd.tx,True)

        #Efirst1 = phi1.corrgauss(sig).Efirst_loc(12,sum(self.tdd.tx.y * self.tdd.tx.y)* self.tdd.tx.dx())
        #Efirst2 = phi2.corrgauss(sig).Efirst_loc(12,sum(self.tdd.tx.y * self.tdd.tx.y)* self.tdd.tx.dx())
        #Efirst3 = phi3.corrgauss(sig).Efirst_loc(12,sum(self.tdd.tx.y * self.tdd.tx.y)* self.tdd.tx.dx())
        #Efirst4 = phi4.corrgauss(sig).Efirst_loc(12,sum(self.tdd.tx.y * self.tdd.tx.y)* self.tdd.tx.dx())
        toa = self.toa_win()
        Ef1 = self.tdd.ch1.Ewin(tau=toa[0], Tint=Tint, sym=sym, dB=dB)
        Ef2 = self.tdd.ch2.Ewin(tau=toa[1], Tint=Tint, sym=sym, dB=dB)
        Ef3 = self.tdd.ch3.Ewin(tau=toa[2], Tint=Tint, sym=sym, dB=dB)
        Ef4 = self.tdd.ch4.Ewin(tau=toa[3], Tint=Tint, sym=sym, dB=dB)

        efirst = np.array([Ef1, Ef2, Ef3, Ef4])
        return  efirst

    def Etau0(self, Tint=1, sym=0.25, dB=True):
        """
        calculate the energy around delay tau0
        """
        Etau01 = self.tdd.ch1.Ewin(tau=self.de[0], Tint=Tint, sym=sym, dB=dB)
        Etau02 = self.tdd.ch2.Ewin(tau=self.de[1], Tint=Tint, sym=sym, dB=dB)
        Etau03 = self.tdd.ch3.Ewin(tau=self.de[2], Tint=Tint, sym=sym, dB=dB)
        Etau04 = self.tdd.ch4.Ewin(tau=self.de[3], Tint=Tint, sym=sym, dB=dB)

        etau0 = np.array([Etau01, Etau02, Etau03, Etau04])

        return etau0

    def ecdf(self, Tnoise=10, rem_noise=True, in_positivity=False, display=False, normalize=True, delay=0):
        """
        ecdf

        calculate energy cumulative density function
        """

        ecdf1, var1 = self.tdd.ch1.ecdf(delay=self.de[0])
        ecdf2, var2 = self.tdd.ch2.ecdf(delay=self.de[1])
        ecdf3, var3 = self.tdd.ch3.ecdf(delay=self.de[2])
        ecdf4, var4 = self.tdd.ch4.ecdf(delay=self.de[3])

        ecdf = np.array([ecdf1, ecdf2, ecdf3, ecdf4])
        var = np.array([var1, var2, var3, var4])

        return ecdf, var

    def tdelay(self):
        """
        Build an array with delay values

        tau0
        """

        tau0 = self.de
        tau_moy = np.array([self.F1.chan_param['tau_moy'],
                            self.F2.chan_param['tau_moy'],
                            self.F3.chan_param['tau_moy'],
                            self.F4.chan_param['tau_moy']])

        tau_rms = np.array([self.F1.chan_param['tau_rms'],
                            self.F2.chan_param['tau_rms'],
                            self.F3.chan_param['tau_rms'],
                            self.F4.chan_param['tau_rms']])

        toa_cum = np.array([self.F1.chan_param['toa_cum'],
                            self.F2.chan_param['toa_cum'],
                            self.F3.chan_param['toa_cum'],
                            self.F4.chan_param['toa_cum']])

        toa_th = np.array([self.F1.chan_param['toa_th'],
                           self.F2.chan_param['toa_th'],
                           self.F3.chan_param['toa_th'],
                           self.F4.chan_param['toa_th']])

        toa_max = np.array([self.F1.chan_param['toa_max'],
                            self.F2.chan_param['toa_max'],
                            self.F3.chan_param['toa_max'],
                            self.F4.chan_param['toa_max']])

        t1 = vstack((tau0, toa_th, toa_cum, toa_max, tau_moy -
                     tau_rms, tau_moy + tau_rms))
        t2 = t1.T
        return(t2)

    def fp(self, alpha=0.1, type='all'):
        """
            Build fingerprint

            alpha is a quantile parameter for taum and taurms calculation

        """
        self.F1 = FP(self, 1, alpha=alpha)
        self.F2 = FP(self, 2, alpha=alpha)
        self.F3 = FP(self, 3, alpha=alpha)
        self.F4 = FP(self, 4, alpha=alpha)

    def outlatex(self, S):
        """
        Measurement output latex

        M.outlatex(S)

        S : a Simulation object

        """
        tt = self.tdelay()
        tt = np.array([[tt[0, 0]], [tt[1, 0]], [tt[2, 0]], [tt[3, 0]]])
        self.show(tt, display=False)
        fig = plt.gcf()
        sp1 = fig.add_subplot(5, 1, 1)
        S.indoor.show(fig, sp1)
        sp1.plot(self.rx[1:, 0], self.rx[1:, 1], 'ob')
        sp1.annotate('Rx1', xy=(self.rx[1, 0], self.rx[1, 1]))
        sp1.annotate('Rx2', xy=(self.rx[2, 0], self.rx[2, 1]))
        sp1.annotate('Rx3', xy=(self.rx[3, 0], self.rx[3, 1]))
        sp1.annotate('Rx4', xy=(self.rx[4, 0], self.rx[4, 1]))
        sp1.plot(hstack((self.tx[0], self.tx[0])),
                 hstack((self.tx[1], self.tx[1])), 'or')
        title(self.Tx_position)
        sz = fig.get_size_inches()
        fig.set_size_inches(sz * 2.3)
        ext1 = '.pdf'
        ext2 = '.tex'
        ext3 = '.png'
#
        if self.tx[2] < 1.3:
            h = 1
        else:
            h = 2
#
               #filename='Tx'+self.Tx_position+'h'+str(h)
        filename = str(h) + self.Tx_position
        filename = filename.replace('P', '')
        fig.savefig('./figures/' + filename + ext3, orientation='portrait')
        fig.clf()
        plt.close(fig)
        plt.clf()
        #fig.savefig('./figures/'+filename+ext1,orientation='portrait')
        #clf()  # very important to close the figure otherwise memory error
#        entete='\\clearpage\\setcounter{page}{1}\\pagestyle{Standard} \n {\\centering  '+ self.Tx_position
#        position = ' $T_{x}=[%5.2f,%5.2f,%5.2f]$ \\par} {\\centering \\par} \n' % (self.tx[0] , self.tx[1] , self.tx[2])
#        tfigure = '\\begin{figure}[htp]\n\\centering\n \\includegraphics[width=14.078cm,height=14.131cm]{'+filename+ext1+'}\n \\end{figure}\n'
#
#        l1='\\begin{center}\n'
#        l11='\\end{center}\n'
#        l2='\\tablehead{}\n'
#        l3='\\begin{supertabular}{|l|l|l|l|l|}\n'
#        l33='\\end{supertabular}\n'
#        l4='\\hline\n'
#        l5='Parameter & $Rx_1$ & $Rx_2$ &$Rx_3$ &$Rx_4$ \\\\\n'
#
#        d1 = '%5.2f' % self.F1.metadata['distance']
#        d2 = '%5.2f' % self.F2.metadata['distance']
#        d3 = '%5.2f' % self.F3.metadata['distance']
#        d4 = '%5.2f' % self.F4.metadata['distance']
#
#        t1 = '%5.2f' % self.F1.metadata['delay']
#        t2 = '%5.2f' % self.F2.metadata['delay']
#        t3 = '%5.2f' % self.F3.metadata['delay']
#        t4 = '%5.2f' % self.F4.metadata['delay']
#
#        ldistance='distance (m) & ' + d1 +'&' + d2 +'&' + d3 +'&' + d4 +' \\\\\n'
#        ldelay   ='$\\tau_0$ (ns) & ' + t1 +'&' + t2 +'&' + t3 +'&' + t4 +' \\\\\n'
#
#        Etot1 = '%5.3f' % self.F1.chan_param['Etot']
#        Etot2 = '%5.3f' % self.F2.chan_param['Etot']
#        Etot3 = '%5.3f' % self.F3.chan_param['Etot']
#        Etot4 = '%5.3f' % self.F4.chan_param['Etot']
#        lEtot='$E_{tot}$ (dBJ) & ' + Etot1 +'&' + Etot2 +'&' + Etot3 +'&' + Etot4 +' \\\\\n'
#
#        Emax1 = '%5.3f' % self.F1.chan_param['Emax']
#        Emax2 = '%5.3f' % self.F2.chan_param['Emax']
#        Emax3 = '%5.3f' % self.F3.chan_param['Emax']
#        Emax4 = '%5.3f' % self.F4.chan_param['Emax']
#        lEmax='$E_{max}$ (dBJ)   & ' + Emax1 +' &' + Emax2 +'&' + Emax3 +' &' + Emax4 +' \\\\\n'
#
#        tau_moy1 = '%5.3f' % self.F1.chan_param['tau_moy']
#        tau_moy2 = '%5.3f' % self.F2.chan_param['tau_moy']
#        tau_moy3 = '%5.3f' % self.F3.chan_param['tau_moy']
#        tau_moy4 = '%5.3f' % self.F4.chan_param['tau_moy']
#        ltau_moy='$\\tau_{m}$ (ns) & ' + tau_moy1 +'&' + tau_moy2 +'&' + tau_moy3 +'&' + tau_moy4 +' \\\\\n'
#
#
#        tau_rms1 = '%5.3f' % self.F1.chan_param['tau_rms']
#        tau_rms2 = '%5.3f' % self.F2.chan_param['tau_rms']
#        tau_rms3 = '%5.3f' % self.F3.chan_param['tau_rms']
#        tau_rms4 = '%5.3f' % self.F4.chan_param['tau_rms']
#        ltau_rms='$\\tau_{rms}$ (ns) & ' + tau_rms1 +'&' + tau_rms2 +'&' + tau_rms3 +'&' + tau_rms4 +' \\\\\n'
#
#        toa_max1 = '%5.3f' % self.F1.chan_param['toa_max']
#        toa_max2 = '%5.3f' % self.F2.chan_param['toa_max']
#        toa_max3 = '%5.3f' % self.F3.chan_param['toa_max']
#        toa_max4 = '%5.3f' % self.F4.chan_param['toa_max']
#        ltoa_max='$\\hat{\\tau}_{0}^{max}$ (ns) & ' + toa_max1 +'&' + toa_max2 +'&' + toa_max3 +'&' + toa_max4 +' \\\\\n'
#
#        toa_th1 = '%5.3f' % self.F1.chan_param['toa_th']
#        toa_th2 = '%5.3f' % self.F2.chan_param['toa_th']
#        toa_th3 = '%5.3f' % self.F3.chan_param['toa_th']
#        toa_th4 = '%5.3f' % self.F4.chan_param['toa_th']
#        ltoa_th='$\\hat{\\tau}_{0}^{th}$ (ns) & ' + toa_th1 +'&' + toa_th2 +'&' + toa_th3 +'&' + toa_th4 +' \\\\\n'
#
#
#        toa_cum1 = '%5.3f' % self.F1.chan_param['toa_cum']
#        toa_cum2 = '%5.3f' % self.F2.chan_param['toa_cum']
#        toa_cum3 = '%5.3f' % self.F3.chan_param['toa_cum']
#        toa_cum4 = '%5.3f' % self.F4.chan_param['toa_cum']
#        ltoa_cum='$\\hat{\\tau}_{0}^{cum}$ (ns) & ' + toa_cum1 +'&' + toa_cum2 +'&' + toa_cum3 +'&' + toa_cum4 +' \\\\\n'
#
#        fd = open('./figures/'+filename+ext2,'w')
#        fd.write(entete+position)
#        fd.write(tfigure)
#        fd.write(l1)
#        fd.write(l2)
#        fd.write(l3)
#        fd.write(l4)
#        fd.write(l5)
#        fd.write(l4)
#        fd.write(ldistance)
#        fd.write(l4)
#        fd.write(ldelay)
#        fd.write(l4)
#        fd.write(l4)
#        fd.write(lEtot)
#        fd.write(l4)
#        fd.write(lEmax)
#        fd.write(l4)
#        fd.write(ltau_moy)
#        fd.write(l4)
#        fd.write(ltau_rms)
#        fd.write(l4)
#        fd.write(ltoa_th)
#        fd.write(l4)
#        fd.write(ltoa_cum)
#        fd.write(l4)
#        fd.write(ltoa_max)
#        fd.write(l4)
#        fd.write(l33)
#        fd.write(l11)
#        fd.close()


if __name__ == "__main__":
    doctest.testmod()
