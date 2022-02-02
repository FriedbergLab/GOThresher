#!/usr/bin/python
import unittest
import os
from lib import debias
print(os.getcwd()+os.sep+"goa_yeast.gaf")
data = debias.chooseProteinsBasedOnPublications(os.getcwd()+os.sep+"goa_yeast.gaf", 500, 0)

for i in data:
    print(i)
# t1_dic, all_protein_t1 = create_benchmark.read_gaf(cwd+"//tests//1")
# t2_dic, all_protein_t2 = create_benchmark.read_gaf(cwd+"//tests//2")
# NK_dic, LK_dic = create_benchmark.analyze(t1_dic, t2_dic, all_protein_t1)

'''
def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestDebias)
    return suite


class TestDebias(unittest.TestCase):
    def test_not_include(self):
        for protein in ["A2P2R3", "A2P2R4", "A2P2R5", "A2P2R6"]:
            for ontology in NK_dic:
                self.assertNotIn(protein, NK_dic[ontology], "This protein should not be in NK ")
        for protein in ["A2P2R3", "A2P2R4", "A2P2R5", "A2P2R6"]:
            for ontology in LK_dic:
                self.assertNotIn(protein, LK_dic[ontology], "This protein should not be in LK ")
'''

if __name__ == '__main__':
    unittest.main()
