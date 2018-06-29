import unittest
from pylayers.gis.layout import *

L1 = Layout('defstr.lay')
class TestLayout(unittest.TestCase):

    def test_load(self):
        self.assertEqual(L1.Np,12)
        self.assertEqual(L1.Ns,15)

    def test_check(self):
        bc,ds = L1.check()
        self.assertTrue(bc)

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

if __name__ == '__main__':
    unittest.main()
