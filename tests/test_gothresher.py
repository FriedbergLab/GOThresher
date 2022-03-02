#!/usr/bin/python
import datetime
import unittest
import os
import sys
import io
from Bio.UniProt import GOA
from gothresher import gothresher
import numpy
import hashlib

# script_dir = os.path.dirname(sys.argv[0])
# input_file = os.path.join(os.getcwd(), script_dir, "goa_exampleYeast.gaf")

input_file = os.path.abspath("goa_exampleYeast.gaf")

# print(input_file)
# print("scriptdir", script_dir)

# Convert from GAF to required format
gaf_output = GOA._gaf20iterator(open(input_file, "r"))
gothresher.init_globals()
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

# Choose proteins based on evidence codes
data_ec = gothresher.choose_proteins_based_on_evidence_codes(data, ["IPI", "IBA"], None)
# print(data_ec.keys())

# Choose proteins based on references
data_ref = gothresher.choose_proteins_based_on_references(data, ["PMID:9799240"], None)
# print(data_ref)

# Choose annotations based on PL Information Content
data_plic = gothresher.calculate_phillip_lord_information_content(data, 10, None)

# Create protein to go mapping
mapping_obj = gothresher.create_protein_to_go_mapping(data)


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

    def test_not_include_data_ec(self):
        for entry in ["anntn11857", "anntn28782", "anntn36384"]:
            self.assertNotIn(entry, data_ec.keys(), "This annotation should not be in new_data")

    def test_include_data_ec(self):
        self.assertIn("URB2", data_ec["anntn47066"]["DB_Object_Symbol"])
        self.assertIn("GO:0000981", data_ec["anntn46610"]["GO_ID"])

    def test_data_ref(self):
        got_list = []
        for k in data_ref.keys():
            got_list.append(data[k]['DB:Reference'])
        for entry in got_list:
            self.assertNotIn(entry, ["PMID:11927560", "PMID:8663399", "PMID:9472021"])
        self.assertEqual(len(numpy.unique(got_list)), 1)

    def test_print_details_about_data(self):
        captured_out = io.StringIO()
        sys.stdout = captured_out
        gothresher.print_details_about_data(data)
        sys.stdout = sys.__stdout__
        expected_out = "Total number of annotations in the provided Database  47072\n"+\
                       "Total number of unique proteins in the provided Database  6189\n"+\
                       "Total number of unique references in the provided Database  11599\n"
        self.assertEqual(expected_out, captured_out.getvalue())

    def test_data_before_date(self):
        # Choose annotations based on date
        data_date = gothresher.choose_annotations_based_on_date(data, "2011-7-22", None)
        date_list = []
        for k in data_date.keys():
            date_list.append(data[k]['Date'])
        for entry in date_list:
            self.assertLessEqual(entry, datetime.date(2011, 7, 22))

    def test_data_after_date(self):
        # Choose annotations based on date
        data_date = gothresher.choose_annotations_based_on_date(data, None, "2011-7-22")
        date_list = []
        for k in data_date.keys():
            date_list.append(data[k]['Date'])
        for entry in date_list:
            self.assertGreaterEqual(entry, datetime.date(2011, 7, 22))

    def test_data_include_pl_ic(self):
        for entry in ["anntn45985", "anntn39155", "anntn37332"]:
            self.assertIn(entry, data_plic.keys(), "This annotation should be in new_data")

    def test_prot_to_go_mapping(self):
        self.assertIsInstance(mapping_obj[0], dict)
        self.assertIsInstance(mapping_obj[1], list)

    def test_prot_to_go_mapping_include(self):
        for entry in ["SGD_S000003026", "SGD_S000005098", "SGD_S000000845"]:
            self.assertIn(entry, mapping_obj[0].keys())


if __name__ == '__main__':
    unittest.main()
