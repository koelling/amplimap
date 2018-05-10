#make sure we can import from package directory
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) 

import unittest
import string, math
from amplimap.naive_mapper import *

class NaiveMapper(unittest.TestCase):
    def test_align(self):
        ref = "ACGTAAAAAAAAAAAAAAAG"

        #match
        self.assertEqual(align_and_find_cigar('ACGT', ref), ( 0, [ (0, 4) ] ))
        #match at end
        self.assertEqual(align_and_find_cigar('AAAAG', ref), ( 15, [ (0, 5) ] ))
        #self.assertEqual(align_and_find_cigar('AAAAG', ref, reverse = True), ( 15, [ (0, 5) ] ))
        #match with deletion
        self.assertEqual(align_and_find_cigar('AGTAAA', ref), ( 0, [ (0, 1), (2, 1), (0, 5) ] ))
        #match with insertion (but actually ends up getting matched as A skip and C mismatch)
        self.assertEqual(align_and_find_cigar('ACCGT', ref), ( 0, [(4, 1), (0, 4)] ))
        #proper match with insertion 
        self.assertEqual(align_and_find_cigar('ACGGTAAA', ref), ( 0, [ (0, 2), (1, 1), (0, 5) ] ))
        #soft clip start
        self.assertEqual(align_and_find_cigar('GGGACGT', ref), ( 0, [ (4, 3), (0, 4) ] ))
        #soft clip start and insert
        self.assertEqual(align_and_find_cigar('GGGACGGTA', ref), ( 0, [ (4, 3), (0, 2), (1, 1), (0, 3) ] ))
        #start offset
        self.assertEqual(align_and_find_cigar('TAAAAA', ref), ( 3, [ (0, 6) ] ))
        #self.assertEqual(align_and_find_cigar('TAAAAA', ref, reverse = True), ( 3, [ (0, 6) ] ))
        #start offset and change
        self.assertEqual(align_and_find_cigar('TAAAGAA', ref), ( 3, [ (0, 7) ] ))
        #start offset and insert (need at least five to be better than mismatch!)
        self.assertEqual(align_and_find_cigar('TAAAGGGGGAA', ref), ( 3, [ (0, 4), (1, 5), (0, 2) ] ))


if __name__ == '__main__':
    unittest.main()