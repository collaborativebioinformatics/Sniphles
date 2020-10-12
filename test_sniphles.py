import unittest
from sniphles import get_unphased_blocks, PhaseBlock


class TestSniphles(unittest.TestCase):

    def setUp(self):
        # Initiate
        phased_intervals = [[31400000, 31600000], [31800000, 31850000], [32000000, 32300000]]
        self.phase_blocks = [PhaseBlock(
            id='test-id',
            start=interval[0],
            end=interval[1],
            phase=[1, 2],
            status='biphasic'
        ) for interval in phased_intervals]

    def test_unphased_blocks(self):
        unphased_blocks = get_unphased_blocks(self.phase_blocks, 0, 33000000)
        self.assertEqual([[0, 3140000], [3160000, 31800000], [31850000, 32000000], [32300000, 33000000]],
                         unphased_blocks)


if __name__ == '__main__':
    sniphles_suite = unittest.TestLoader().loadTestsFromTestCase(TestSniphles)

    all_tests = unittest.TestSuite([sniphles_suite])
    unittest.TextTestRunner(verbosity=2).run(all_tests)
