import unittest
from Planck import *

class TestPlanckLFIHFI(unittest.TestCase):

    def setUp(self):
        self.Planck = Planck()
        self.lfi = self.Planck.inst['LFI']
        self.hfi = self.Planck.inst['HFI']

    def test_Planck(self):
        self.assertEqual(len(self.Planck.ch), 71)
        
    def test_LFI(self):
        self.assertEqual(self.lfi.name , 'LFI')
        self.assertEqual(len(self.lfi.ch), 22 )
        self.assertEqual(len(self.lfi.d), 44 )
        self.assertEqual(self.lfi['LFI25M'].tag, 'LFI25M')
        self.assertEqual(self.lfi['LFI25M'].arm, 'M')
        self.assertEqual(self.lfi['LFI25M'].RCA, 25)
        self.assertEqual(self.lfi['LFI27M'][0].tag, 'LFI27M-00')
        self.assertEqual(self.lfi['LFI27M'][0].ch.tag, 'LFI27M')
        self.assertEqual(self.lfi['LFI28S'].tag, 'LFI28S')
        self.assertEqual(self.lfi.ch[0].tag, 'LFI18M')
        self.assertEqual(self.lfi.ch[0].tag, 'LFI18M')

    def test_HFI(self):
        self.assertEqual(len(self.hfi.ch), 49)
        self.assertEqual(self.hfi['217-8a'].tag, '217-8a')
        self.assertEqual(self.hfi['545-4'].tag, '545-4')
