import importlib.util
import pathlib
import unittest

import numpy as np
import pandas as pd

modfile = pathlib.Path(__file__).resolve().parents[1]/'02_spatial_neighborhood_qc.py'
spec = importlib.util.spec_from_file_location('spatial_qc', modfile)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class NeighborhoodTests(unittest.TestCase):
    def example(self):
        return pd.DataFrame({
            'cell_id':['a','b','c','d','e','f'],
            'x_um':[0,10,100,110,200,220],
            'y_um':[0,0,0,0,0,0],
            'CK':[5,1,6,1,4,0.1],
            'CD8':[0.1,1,0.1,0.1,1,0.02],
            'Autofluorescence':[1,2,3,4,5,6],
        })

    def test_counts_and_neighborhood(self):
        r, per_ck, g, l, masks = mod.summarize(self.example(), radius=30, permutations=19, seed=5)
        self.assertEqual(r['counts']['ck_only'],2)
        self.assertEqual(r['counts']['cd8_only'],1)
        self.assertEqual(r['counts']['dual'],1)
        self.assertEqual(r['counts']['low_low'],2)
        self.assertEqual(r['observed']['ck_only_with_any_cd8_only_within_radius_fraction'],.5)
        self.assertEqual(sorted(per_ck['n_CD8_only_within_radius'].tolist()),[0,1])
        self.assertEqual(len(g),19)
        self.assertEqual(len(l),19)
        self.assertEqual(r['af_audit']['dual']['n'],1)

    def test_duplicate_id_rejected(self):
        df = self.example()
        df.loc[1,'cell_id'] = 'a'
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            mod.summarize(df)

    def test_single_group_rejected(self):
        df = self.example()
        df['CD8'] = 0.
        with self.assertRaisesRegex(ValueError, 'one candidate group'):
            mod.summarize(df)


if __name__ == '__main__':
    unittest.main()
