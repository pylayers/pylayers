#-*- coding:Utf-8 -*-
from scipy.spatial import Delaunay

from pylayers.gis.layout import *
import shapely.geometry as sh
L = Layout('WHERE2_2.ini',force=True)
#L.dumpr()
L.build('t')


def polyplot(poly):
    fig,ax=L.showG('s')
    color=['r','b','g']*10
    for ip, p in enumerate(poly):
        fig,ax = p.plot(fig=fig,ax=ax,color=color[ip],alpha =0.5)



ch = L.ma.convex_hull
P = ch.difference(L.ma)
polys = []
if isinstance(P,sh.MultiPolygon):
    for p in P:
        if p.area > 5e-2:
            polys.append(geu.Polygon(p))
            polys[-1].setvnodes(L)


ncy = max(L.Gt.nodes())+1
for p in polys:
    uaw = np.where(p.vnodes == 0)
    for aw in uaw :
        print p.vnodes[aw-1][0], p.vnodes[aw+1][0]
        awid = L.add_segment(p.vnodes[aw-1][0], p.vnodes[aw+1][0], name='AIR')
        p.vnodes[aw] = awid
        # G = nx.subgraph(L.Gs,p.vnodes)
        # G.pos = {}
        # G.pos.update({l: L.Gs.pos[l] for l in p.vnodes})
        # cy  = cycl.Cycle(G)
        # L.Gt.add_node(ncy,cycle=cy)
        # L.Gt.pos[ncy] = tuple(cy.g)
        # L.Gt.node[ncy]['polyg'] = p
        # L.Gt.node[ncy]['isopen'] = True
        # L.Gt.node[ncy]['indoor'] = False
        # import ipdb
        # ipdb.set_trace()
        # for k in L.Gt.nodes():
        #     if (k != ncy) and (k != 0):
        #         print k
        #         vnodes0 = np.array(L.Gt.node[ncy]['cycle'].cycle)
        #         vnodes1 = np.array(L.Gt.node[k]['cycle'].cycle)
        #         #
        #         # Connect Cycles if they share at least one segments
        #         #
        #         intersection_vnodes = np.intersect1d(vnodes0, vnodes1)

        #         if len(intersection_vnodes) > 1:
        #             segment = intersection_vnodes[np.where(intersection_vnodes>0)]
        #             L.Gt.add_edge(ncy, k,segment= segment)

        # for v in filter(lambda x: x>0,p.vnodes):
        #     # add new ncycle to Gs for the new airwall
        #     # that new airwall always separate the new created cycle
        #     # and the outdoor cycle
        #     if v == awid :
        #         L.Gs.node[awid]['ncycles']=[ncy,0]
        #     # other wise update the cycles seen by semengts
        #     else :
        #         cy = L.Gs.node[v]['ncycles'].pop()
        #         # if the pop cycle is the outdoor cycle, replace it with the new cycle
        #         if cy == 0:
        #             L.Gs.node[v]['ncycles'].append(ncy)
        #         # else replace old value with [pos cycle , new cycle]
        #         else:
        #             L.Gs.node[v]['ncycles']=[cy,ncy]
        ncy=ncy+1
# tcc, nn = L.ma.ptconvex()
# utconvex = np.nonzero(tcc == -1)[0]
# # all possible diffracting point (in and out of cycle)
# utsconvex = np.nonzero(abs(tcc) == 1)[0]
# if len(utconvex) != 0:
#     # get points ID in the cycle
#     uus = filter(lambda x: x<0,L.ma.vnodes)
#     # get point convex ID
#     uc = np.array(uus)[utconvex]
#     ucs = np.array(uus)[utsconvex]
#     puc = array(map(lambda x: L.Gs.pos[x], uc))
#     pucs = array(map(lambda x: L.Gs.pos[x], ucs))
#     uc = np.hstack((uc,uc[0]))
#     ucs = np.hstack((ucs,ucs[0]))
#     pucs = np.vstack((pucs,pucs[0]))

#     if len(ucs) >2:
#         trid=Delaunay(pucs)
#         tri =trid.simplices
#         aucs = np.arange(len(ucs))
#         # filter tri in the cycle
#         kt = []
#         pkt = []
#         polys = []
#         naw = []
#         for t in tri:
#             ts = geu.Polygon(pucs[t])
#             #check if inside the original polygon
#             # U = L.Gt.node[n]['polyg'].contains(ts)
#             U = L.ma.intersection(ts)
#             ats = ts.area
#             if U.area < (1*ats/100):
#                 #pkt.append(pucs[t])
#                 #if 2 convex nodes are directly following in the list 
#                 # of all convex nodes (ucs) they are already connected
#                 # otherwise, an airwall has to be create.
#                 # 
#                 # ucs[t]

#                 daucs = np.diff(aucs[t])
#                 # search where an airwall is required
#                 # ncp : not connected points
#                 ncp = np.where(daucs != 1)[0]
#                 for i in ncp:
#                     # keep trace of created airwalls, because some 
#                     # of them will be destroyed in step 3
#                     print ucs[t][i],ucs[t][i+1]
#                     naw.append(L.add_segment(ucs[t][i],ucs[t][i+1],name='AIR'))
#                 kt.append(t) 
#                 polys.append(ts)
#         # L.showG('s',labels=True, nodelist=uc)
#         # polyplot(polys)
#         cpolys = []
#         nbpolys = len(polys)
#         while polys !=[]:
#             # import ipdb
#             # ipdb.set_trace()
#             p = polys.pop(0)
#             for ip2,p2 in enumerate(polys):
#                 conv=False
#                 inter = p.intersection(p2)
#                 # if 2 triangles have a common segment
#                 pold = p
#                 if isinstance(inter,sh.LineString):
#                     p = p + p2
#                     if p.isconvex():
#                         polys.pop(ip2)
#                         polys.insert(0,p)
#                         conv=True
#                         break
#                     else:
#                         # if pold not in cpolys:
#                         #     cpolys.append(pold)
#                         p = pold
#             # if (ip2 >= len(polys)):# and (conv):
#             # if conv :
#             #     if p not in cpolys:
#             #         cpolys.append(p)
#             if not conv:#else:
#                 if pold not in cpolys:
#                     cpolys.append(pold)
#             if len(polys) == 0:
#                 cpolys.append(p)
#         ####
#         #### 4. ensure the correct vnode numerotaion of the polygons
#         #### and remove unecessary airwalls
#         ncpol = []
#         vnodes=[]
#         # for p in cpolys:
#         #     ptmp = geu.Polygon(L.Gt.node[n]['polyg'].intersection(p))
#         #     ptmp.setvnodes(L)
#         #     ncpol.append(ptmp)
#         #     vnodes.extend(ptmp.vnodes)
#         # # air walls to be deleted (because origin Delaunay triangle
#         # # has been merged )
#         # daw = filter(lambda x: x not in vnodes,naw)
#         # [L.del_segment(d) for d in daw]
# L.showG('s')
# polyplot(cpolys)
