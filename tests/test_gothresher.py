#!/usr/bin/python
import unittest
import os

from Bio.UniProt import GOA
from gothresher import gothresher

input_file = "example_data" + os.sep + "goa_exampleYeast.gaf"
# print(input_file)

# Convert from GAF to required format
gaf_output = GOA._gaf20iterator(open(input_file, "r"))
data = gothresher.convert_from_gaf_to_required_format(gaf_output)
gaf_output.close()

prot_list = []
for key in data.keys():
    prot_list.append(data[key]["DB_Object_Symbol"])

# Choose proteins based on publications
new_data = gothresher.choose_proteins_based_on_publications(data, 50, None)
new_prot_list = []
for key in new_data.keys():
    new_prot_list.append(new_data[key]["DB_Object_Symbol"])

# Count proteins
prot_count = gothresher.count_proteins(data)

# GO terms frequency
go_freq = gothresher.freq_go_term(data)

# Choose GO based on aspect
data_aspect = gothresher.choose_go_based_on_aspect(data, "P")

# Choose GO based on assigned by
data_ab = gothresher.choose_go_based_on_assigned_by(data, "SGD", None)
# print(data.keys() - data_ab.keys())


def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestGOThresher)
    return suite


class TestGOThresher(unittest.TestCase):
    def test_not_include_data(self):
        for entry in ["anntn46881", "anntn31699", "anntn17651"]:
            self.assertNotIn(entry, new_data.keys(), "This annotation should not be in new_data")
            # print(data["anntn17651"])

    def test_include_data(self):
        self.assertIn("LIA1", data["anntn46881"]["DB_Object_Symbol"], "protein LIA1 should be present in data")
        self.assertIn("TIF4632", data["anntn17651"]["DB_Object_Symbol"], "protein TIF4632 should should be present in "
                                                                         "data")
        self.assertNotIn("TRT2", new_prot_list, "protein TRT2 should not be present in new_data")
        self.assertNotIn("SUP16", new_prot_list, "protein SUP16 should not be present in new_data")

    def test_count_protein(self):
        self.assertEqual(prot_count, 6189)

    def test_go_freq_include(self):
        for entry in ["GO:0005635", "GO:0031304", "GO:0001228"]:
            self.assertIn(entry, go_freq.keys())

    def test_go_freq_equal(self):
        self.assertEqual(go_freq["GO:0001228"], 96)
        self.assertEqual(go_freq["GO:0005635"], 74)
        self.assertEqual(go_freq["GO:0031304"], 4)

    def test_not_include_data_aspect(self):
        for entry in ["anntn26561", "anntn33749", "anntn36384"]:
            self.assertNotIn(entry, data_aspect.keys(), "This annotation should not be in new_data")

    def test_include_data_aspect(self):
        # data_aspect["anntn47066"], data_aspect["anntn11857"]
        self.assertIn("URB2", data_aspect["anntn47066"]["DB_Object_Symbol"])
        self.assertIn("GO:0006364", data_aspect["anntn11857"]["GO_ID"])

    def test_not_include_data_ab(self):
        for entry in ["anntn30185", "anntn28782", "anntn16907"]:
            self.assertNotIn(entry, data_ab.keys(), "This annotation should not be in new_data")

    def test_include_data_ab(self):
        self.assertIn("URB2", data_ab["anntn47066"]["DB_Object_Symbol"])
        self.assertIn("GO:0006364", data_ab["anntn11857"]["GO_ID"])


if __name__ == '__main__':
    unittest.main()
