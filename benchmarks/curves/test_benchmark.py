"""Small harness checks; these do not launch a benchmark suite or Sage."""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parents[1]/'src'))
from mpmath import mp
import generapy
from curve_cases import catalog, x_transform
from validation import basis_change, compare, symplectic
from benchmark_curves import report, worker, REPO, acquire_run_lock


class HarnessChecks(unittest.TestCase):
    def test_clustered_exact_input_with_guard_digits(self):
        case = catalog()['hyper-g2-close-00001']
        for digits in (30, 50):
            records = []
            for extra in (15, 30):
                with mp.workdps(digits+extra):
                    request = dict(case=case,engine='generapy',bits=mp.prec,repo=str(REPO))
                data = worker(request,'unused',10)
                self.assertEqual(data['status'],'ok',data)
                self.assertEqual(data['working_bits'], request['bits'])
                records.append(data)
            self.assertTrue(compare(mp,*records,digits)['passed'])

    def test_concurrent_runs_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'lock'
            with acquire_run_lock(path):
                with self.assertRaises(RuntimeError):
                    acquire_run_lock(path)
            with acquire_run_lock(path):
                pass

    def test_basis_alignment_is_exact(self):
        source = [[[0,0,'1']],[[1,0,'1']]]
        target = [[[0,0,'2']],[[0,0,'1'],[1,0,'1']]]
        self.assertEqual(basis_change(source,target),[[Q(2),Q(0)],[Q(1),Q(1)]])
        with self.assertRaises(ValueError):
            basis_change(source, [[[0,1,'1']],[[1,0,'1']]])

    def test_no_fitted_scale_or_nonsymplectic_map(self):
        reference = dict(genus=1,basis=[[[0,0,'1']]],periods=[[['1','0'],['0','1']]])
        sample = dict(genus=1,basis=[[[0,0,'2']]],periods=[[['2','0'],['0','2']]])
        self.assertTrue(compare(mp,sample,reference,20)['passed'])
        sample['periods'] = reference['periods']
        self.assertFalse(compare(mp,sample,reference,20)['passed'])
        self.assertFalse(symplectic([[2,0],[0,2]]))
        self.assertTrue(symplectic([[0,1],[-1,0]]))

    def test_catalog_and_transform(self):
        cases = catalog()
        self.assertGreaterEqual(len(cases),26)
        self.assertNotEqual(cases['general-trig-g3']['terms'],cases['general-trig-dense-g3']['terms'])
        self.assertEqual(x_transform({(2,0):1},1,2),{(0,0):Q(1,4),(1,0):Q(-1,2),(2,0):Q(1,4)})
        self.assertTrue(cases['singular-g1-double-point']['expected_rejection'])

    def test_worker_fraction_import_path_and_timeout(self):
        case = catalog()['hyper-g1-lemniscatic']
        request = dict(case=case,engine='generapy',bits=70,repo=str(REPO))
        data = worker(request,'unused',10)
        self.assertEqual(data['status'],'ok',data)
        self.assertEqual(data['actual_engine'],'hyperelliptic')
        self.assertTrue(Path(data['generapy_path']).is_relative_to(REPO))
        self.assertEqual(worker(request,'unused',0.00001)['status'],'timeout')

    def test_report_never_ranks_failed_accuracy_and_labels_hyperelliptic(self):
        rows = [dict(case='a',digits=20,implementation='generapy',actual_engine='hyperelliptic',
                     status='ok',seconds=1,accuracy={'status':'failed'}),
                dict(case='a',digits=20,implementation='sage',status='ok',seconds=2,
                     accuracy={'status':'passed'})]
        text = report(dict(samples=rows,references=[]))
        self.assertIn('hyperelliptic',text)
        self.assertNotIn('| a | 20 | hyperelliptic | 2.000 |',text)


if __name__=='__main__':
    unittest.main()
