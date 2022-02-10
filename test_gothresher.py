#!/usr/bin/python
import unittest
import os

from Bio.UniProt import GOA
from gothresher import gothresher

input_file = os.getcwd() + os.sep + ".." + os.sep + "example_data" + os.sep + "goa_exampleYeast.gaf"
print(input_file)

gaf_output = GOA._gaf20iterator(open(input_file, "r"))
data = gothresher.convert_from_gaf_to_required_format(gaf_output)
gaf_output.close()

prot_list = []
for key in data.keys():
    prot_list.append(data[key]["DB_Object_Symbol"])

new_data = gothresher.choose_proteins_based_on_publications(data, 50, None)
new_prot_list = []
for key in new_data.keys():
    new_prot_list.append(new_data[key]["DB_Object_Symbol"])


def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestChooseProteinsBasedOnPublications)
    return suite


class TestChooseProteinsBasedOnPublications(unittest.TestCase):
    def test_not_include(self):
        for entry in ["anntn46881", "anntn31699", "anntn17651"]:
            self.assertNotIn(entry, new_data.keys(), "This annotation should not be in new_data")
            # print(data["anntn17651"])

    def test_include_data(self):
        self.assertIn("LIA1", data["anntn46881"]["DB_Object_Symbol"], "protein LIA1 should be present in data")

        self.assertIn("TIF4632", data["anntn17651"]["DB_Object_Symbol"], "protein TIF4632 should should be present in "
                                                                         "data")
        self.assertNotIn("TRT2", new_prot_list, "protein TRT2 should not be present in new_data")

        self.assertNotIn("SUP16", new_prot_list, "protein SUP16 should not be present in new_data")


if __name__ == '__main__':
    unittest.main()
