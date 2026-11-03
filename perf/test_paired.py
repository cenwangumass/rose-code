"""Decision-rule checks; run: python3 -m unittest discover -s perf -p test_paired.py."""
import unittest
from paired import summarize


def samples(ab, aa):
    return [dict(kind=kind, order='AB' if i % 2 == 0 else 'BA', log_ratio=value)
            for kind, values in [('AB', ab), ('AA', aa)]
            for i, value in enumerate(values)]


class PairedScreenTests(unittest.TestCase):
    def test_identical_timers_do_not_show_gain(self):
        self.assertFalse(summarize(samples([0.] * 32, [0.] * 32))['positive_screen'])

    def test_stable_gain_with_quiet_controls(self):
        self.assertTrue(summarize(samples([-.02] * 32, [0.] * 32))['positive_screen'])

    def test_execution_order_confound_is_rejected(self):
        ab = [-.2 if i % 2 == 0 else .1 for i in range(32)]
        self.assertFalse(summarize(samples(ab, [0.] * 32))['positive_screen'])

    def test_biased_same_binary_control_is_rejected(self):
        self.assertFalse(summarize(samples([-.02] * 32, [-.01] * 32))['positive_screen'])
