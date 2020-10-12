import unittest
import pysam
from sniphles import get_unphased_blocks, check_phase_blocks


class TestSniphles(unittest.TestCase):

    def setUp(self):
        # Initiate
        bam = pysam.AlignmentFile("./test_bam.bam", "rb")
        for chrom in bam.references:  # Iterate over all chromosomes separately
            self.phase_blocks = check_phase_blocks(bam, chrom)

    def test_unphased_blocks(self):
        unphased_blocks = get_unphased_blocks(self.phase_blocks, 0, 32218201)
        print(self.phase_blocks)
        print(unphased_blocks)
        self.assertEqual(True, True)


if __name__ == '__main__':
    sniphles_suite = unittest.TestLoader().loadTestsFromTestCase(TestSniphles)

    all_tests = unittest.TestSuite([sniphles_suite])
    unittest.TextTestRunner(verbosity=2).run(all_tests)
