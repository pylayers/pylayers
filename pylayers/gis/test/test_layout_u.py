import unittest
from pylayers.gis.layout import *
L1 = Layout('defstr.lay')
class TestLayout(unittest.TestCase):

    def test_add_fnod(self):
        L = Layout('defstr.lay')
        L.add_fnod(p=(10,10))
        self.assertEqual(L.Np,13)

    def test_add_furniture(self):
        L = Layout('defstr.lay')
        L.add_furniture(name='R1_C',
                         matname='PARTITION', origin=(5., 5.),
                         zmin=0.,
                         height=1.,
                         width=1.,
                         length=1.,
                         angle=0.)

    def test_add_nfpe(self):
        L = Layout('defstr.lay')
        L.add_nfpe(-8,7,6)
        self.assertEqual(L.Np,13)

    def test_angleonlink(self):
        data1 = L1.angleonlink(np.array([2,2.5]),np.array([8,4]))
        data2 = L1.angleonlink3(np.array([2,2.5,1.5]),np.array([8,4,1.5]))
        print(data1)
        print(data2)
        # The 3d intersection releaves ambiguity on vertical segments
        # Problem at jonction between seg here 2 segments are found 
        #
        data1 = L1.angleonlink(np.array([2,2]),np.array([8,4]))
        data2 = L1.angleonlink3(np.array([2,2,1.5]),np.array([8,4,1.5]))
        print(data1)
        print(data2)

    def test_boundary(self):
        L = Layout('defstr.lay')
        L.boundary()

    def test_build(self):
        L = Layout('defstr.lay')
        L.build()

    def test_cleanup(self):
        L = Layout('defstr.lay')
        L.add_fnod(p=(10,10))
        L.cleanup()
        self.assertEqual(L.Np,12)

    def test_load(self):
        L = Layout('defstr.lay')
        self.assertEqual(L.Np,12)
        self.assertEqual(L.Ns,15)

    def test_check(self):
        bc,ds = L1.check()
        self.assertTrue(bc)

    def test_check2(self):
        L = Layout('defstr.lay')
        L.build()
        tseg = L.check2()
        L.build()
        L.check_Gi()

    def test_have_subseg(self):
        self.assertTrue(L1.have_subseg(1))
        self.assertTrue(L1.have_subseg(2))
        self.assertTrue(L1.have_subseg(3))
        self.assertFalse(L1.have_subseg(4))
        self.assertFalse(L1.have_subseg(5))
        self.assertFalse(L1.have_subseg(6))
        self.assertFalse(L1.have_subseg(7))
        self.assertFalse(L1.have_subseg(8))
        self.assertFalse(L1.have_subseg(9))
        self.assertFalse(L1.have_subseg(10))
        self.assertFalse(L1.have_subseg(11))

    def test_add_pons(self):
        L = Layout('defstr.lay')
        L.add_pons(1,alpha=0.6)
        self.assertEqual(L.Np,13)

    def test_isseg(self):
        self.assertTrue(L1.isseg(-8,-7))

    def test_ispoint(self):
        pto = np.array([[0],[0]])
        num = L1.ispoint(pto)
        self.assertEqual(num,-1)

    def test_seg_intersection(self):
        pt1 = L1.Gs.pos[-8]
        pt2 = L1.Gs.pos[-7]
        liseg,lipsh = L1.seg_intersection(**{'ta':pt1,'he':pt2})

    def test_clip(self):
        seglist =  L1.clip(2,8,2,4)
        self.assertEqual(sum(seglist),10)

    def test_cy2pt(self):
        L = Layout('defstr.lay')
        L.build()
        pt = L.cy2pt(2)

    def test_geomfile(self):
        L1.geomfile()

    def test_DLRosm(self):
        L = Layout('DLR.osm')

if __name__ == '__main__':
    unittest.main()

# add_pnod
# add_pons
# add_segment
# angleonlinkold
# build
# buildGi
# buildGr
# buildGt
# buildGv
# buildGw
# _convex_hull
# cy2pt
# cycleinline
# _delaunay
# del_points
# del_segment
# diag
# distwall
# dumpr
# dumpw
# ed2nd
# edit_seg
# exportosm
# extrseg
# facet3D
# facets3D
# filterGi
# _find_diffractions
# find_edgelist
# g2npy
# geomfile
# getangles
# get_diffslab
# get_paths
# get_points
# get_Sg_pos
# get_zone
# have_subseg
# importosm
# importres
# importshp
# info
# info_segment
# intercy
# _interlist
# isindoor
# ispoint
# isseg
# layerongrid
# layeronlink
# load
# loadfur
# load_modif
# ls
# mask
# _merge_polygons
# merge_segment
# nd2seg
# numseg
# off_overlay
# offset_index
# onseg
# outputGi
# outputGi_mp
# outputGi_new
# plot
# plot_segments
# pltlines
# pltpoly
# pltvnodes
# point_touches_seg
# polysh2geu
# pt2cy
# pt2ro
# ptGs2cy
# ptin
# randTxRx
# repair
# room2nodes
# room2segments
# rotate
# save
# scl_overlay
# seg2pts
# seg2ro
# seginframe
# seginframe2
# seginline
# segpt
# seguv
# show
# show3
# _show3
# showG
# _showGi
# showGs
# _showGt
# _showGv
# show_layer
# show_nodes
# show_seg1
# show_segment
# showSig
# signature
# subseg
# thwall
# translate
# _triangle
# _triangle_old
# updateshseg
# _updGsncy
# visilist
# visi_papb
# _visual_check
# waypointGw
# wedge
# wedge2
