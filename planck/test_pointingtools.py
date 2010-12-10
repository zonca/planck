import unittest

from pointingtools import *

def AHF_path_to_od(path):
    return int(path.split('/')[-2])

class TestPointingTools(unittest.TestCase):


    def test_AHF_btw_od(self):
        obt = [1629623011.4397736, 1629883390.4398956]
        files = AHF_btw_OBT(obt)
        ods = map(AHF_path_to_od, files)
        self.assertEqual(ods, [100, 101, 102, 103])

if __name__ == '__main__':
    unittest.main()
