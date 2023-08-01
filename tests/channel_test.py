import unittest

from fatools.lib.fautil.mixin2 import FSAMixIn


class GetSaturatedRtimesTest(unittest.TestCase):

    #mocks up a channels and alleles
    class Channel:
        def __init__(self, alleles):
            self.alleles = alleles
        def get_alleles(self):
            return self.alleles
        
    class Allele:
        def __init__(self, rfu, rtime):
            self.rfu = rfu
            self.rtime = rtime

    high_alleles = [Allele(29000,100),Allele(29001,200), Allele(35000,300), Allele(40000,400)]
    low_alleles = [Allele(28999,100),Allele(3000,200), Allele(1000,300), Allele(500,400)]  

    def setUp(self):
        self.fsa = FSAMixIn()

    def test_get_saturated_rtimes_empty_channels(self):
        self.fsa.channels = []
        saturated_rtimes = self.fsa.get_saturated_rtimes()
        self.assertEqual(saturated_rtimes, [])

    def test_get_saturated_rtimes_no_saturated_alleles(self):
        self.fsa.channels = [self.Channel(self.low_alleles)]
        saturated_rtimes = self.fsa.get_saturated_rtimes()
        self.assertEqual(saturated_rtimes, [])

    def test_get_saturated_rtimes_some_saturated_allele(self):
        self.fsa.channels = [self.Channel(self.low_alleles),self.Channel(self.high_alleles)]
        saturated_rtimes = self.fsa.get_saturated_rtimes()
        self.assertEqual(saturated_rtimes, [100,200,300,400])

if __name__ == "__main__":
    unittest.main()