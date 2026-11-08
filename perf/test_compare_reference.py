import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('compare_reference', Path(__file__).with_name('compare_reference.py'))
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)


class BenchmarkContracts(unittest.TestCase):
    def test_valid_nonoptimal_exits_are_timed(self):
        for status, method, code in [('numerical_failure', 'dual', 3), ('suboptimal', 'barrier', 3),
                                     ('iteration_limit', 'pdhg', 3), ('infeasible', 'dual', 0),
                                     ('infeasible', 'barrier', 3), ('optimal', 'primal', 0)]:
            bench.validate_rose_exit({'status': status}, code, method)
            with self.assertRaises(ValueError):
                bench.validate_rose_exit({'status': status}, 3 if code == 0 else 0, method)
        with self.assertRaises(ValueError):
            bench.validate_rose_exit({'status': 'optimal'}, 1, 'primal')

    def test_method_specific_iteration_limits(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'bench.txt'
            p.write_text('g model.mps --method pdhg --pdhg-iter-limit 20000\n'
                         'b model.mps --method barrier\n'
                         'd model.mps --method dual --profile release --max-iterations 17\n')
            g, b, d = bench.rows(p)
            self.assertEqual(g['reference_extra'], {'OutputFlag': 0, 'PDHGIterLimit': 20000})
            self.assertEqual(b['reference_extra'], {'OutputFlag': 0, 'BarIterLimit': 1000})
            self.assertEqual(d['reference_profile'], 'dual-release')
            self.assertEqual(d['reference_extra']['IterationLimit'], 17)


if __name__ == '__main__':
    unittest.main()
