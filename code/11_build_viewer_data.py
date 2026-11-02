#!/usr/bin/env python3
"""Build docs/viewer_data.js for the variant browser (https://taehojo.github.io/rarevariants/).

Input, both in data/:
  Table_S2_variant_data.csv  the 9,943 analyzed variants (allele count >= 3, deduplicated, AlphaGenome-scored),
                             one retained stratum record per variant, as in Supplementary Table S3
  Table_S1_gene_list.csv     the 85 genes with the cell-type category used in the paper

Per gene the browser shows every analyzed variant, the variant count by ancestry stratum of the retained record,
the percentage of case-enriched variants (case-only variants counted as case-enriched, as in Table 2 and
Supplementary Table S2), the median case-control ratio of variants with a finite ratio, and the gene-level
interaction ratio of Supplementary Table S6 (RNA-seq; variants with a finite case-control ratio; high-effect =
score at or above the gene's 90th percentile; genes with at least 10 such variants and at least 3 in each group).

Usage: python code/11_build_viewer_data.py [--data data] [--out docs/viewer_data.js]
"""
import argparse
import json
import math
import os

import numpy as np
import pandas as pd

GVC_CAUSAL = {'ABCA7', 'ABI3', 'ADAM10', 'APBB3', 'APP', 'BIN1', 'CASP7', 'CR1', 'EPHA1', 'MS4A6A', 'PICALM',
              'PILRA', 'PLCG2', 'PSEN1', 'PSEN2', 'RIN3', 'SORL1', 'SPI1', 'TREM2'}
MODALITIES = [('rna_seq_effect', 'RNA-seq'), ('cage_effect', 'CAGE'), ('dnase_effect', 'DNase'),
              ('chip_histone_effect', 'Histone'), ('chip_tf_effect', 'ChIP-TF'),
              ('splice_junctions_effect', 'Splice-Jn'), ('splice_sites_effect', 'Splice-Site'),
              ('splice_site_usage_effect', 'Splice-Usage')]
POPULATIONS = ['NHW', 'AA', 'Asian', 'Hispanic']


def r(x, nd):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    if isinstance(x, float) and math.isinf(x):
        return 'inf'
    return round(float(x), nd)


def gene_ir(v):
    """Supplementary Table S6 interaction ratio for one gene."""
    v = v[np.isfinite(v['cc_ratio']) & v['rna_seq_effect'].notna()]
    if len(v) < 10:
        return None
    high = v['rna_seq_effect'] >= v['rna_seq_effect'].quantile(0.90)
    if high.sum() < 3 or (~high).sum() < 3:
        return None
    ce = v['cc_ratio'] > 1
    hp, lp = ce[high].mean() * 100, ce[~high].mean() * 100
    if lp == 0:
        return None
    return hp / lp, hp, lp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default='data')
    ap.add_argument('--out', default=os.path.join('docs', 'viewer_data.js'))
    a = ap.parse_args()
    var = pd.read_csv(os.path.join(a.data, 'Table_S2_variant_data.csv'))
    genes_tab = pd.read_csv(os.path.join(a.data, 'Table_S1_gene_list.csv'))
    cell_type = dict(zip(genes_tab['gene_name'], genes_tab['Cell_Type']))

    genes = {}
    for g in sorted(genes_tab['gene_name']):
        v = var[var['gene_name'] == g].copy()
        finite = v['cc_ratio'][np.isfinite(v['cc_ratio'])]
        ir = gene_ir(v)
        rows = []
        for _, x in v.iterrows():
            d = {'id': x['variant_id'], 'p': int(x['pos']), 'ref': x['REF'], 'alt': x['ALT'],
                 'cc': r(x['cc_ratio'], 2), 'enr': x['enrichment'], 'pop': x['population']}
            for col, short in MODALITIES:
                val = r(x[col], 2)
                if val is not None:
                    d[short] = val
            rows.append(d)
        rows.sort(key=lambda d: max([d[s] for _, s in MODALITIES if s in d] or [0]), reverse=True)
        n_ce = int(v['enrichment'].isin(['case_enriched', 'case_only']).sum())
        genes[g] = {
            'name': g, 'chr': int(v['chr_num'].iloc[0]) if len(v) else None,
            'start': int(v['pos'].min()) if len(v) else 0, 'end': int(v['pos'].max()) if len(v) else 0,
            'cellType': cell_type.get(g, 'Ubiquitous'), 'gvc': 'Causal' if g in GVC_CAUSAL else 'GWAS Loci',
            'ir': ir[0] if ir else None, 'irHighPct': ir[1] if ir else None, 'irLowPct': ir[2] if ir else None,
            'pctCE': round(100 * n_ce / len(v), 1) if len(v) else None,
            'meanCC': r(finite.mean(), 3) if len(finite) else None,
            'medianCC': r(finite.median(), 3) if len(finite) else None,
            'pops': {p: {'total': int((v['population'] == p).sum())} for p in POPULATIONS},
            'variants': rows, 'nVar': int(len(v)),
        }

    by_ct, by_gvc = {}, {}
    for g in genes.values():
        for key, dct in ((g['cellType'], by_ct), (g['gvc'], by_gvc)):
            dct.setdefault(key, {'genes': 0, 'variants': 0})
            dct[key]['genes'] += 1
            dct[key]['variants'] += g['nVar']
    pct = {}
    for col, short in MODALITIES:
        vals = np.sort(var[col].dropna().round(2).values)
        n = len(vals)
        pct[short] = {k: r(vals[int(n * q)], 2) for k, q in (('p50', .5), ('p75', .75), ('p90', .9), ('p95', .95), ('p99', .99))}
        pct[short]['max'] = r(vals[-1], 2)
    summary = {'totalGenes': len(genes), 'populations': POPULATIONS, 'modalities': [s for _, s in MODALITIES],
               'modPctiles': pct, 'byCellType': by_ct, 'byGVC': by_gvc}

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, 'w') as f:
        f.write('// Generated by code/11_build_viewer_data.py from data/Table_S2_variant_data.csv\n')
        f.write('// (the 9,943 analyzed variants) and data/Table_S1_gene_list.csv; all 8 AlphaGenome modalities\n\n')
        f.write('var VIEWER_DATA={\ngenes:')
        f.write(json.dumps(genes, separators=(',', ':')))
        f.write(',\nsummary:')
        f.write(json.dumps(summary, separators=(',', ':')))
        f.write('\n};\n')
    print('genes', len(genes), '| variants', sum(g['nVar'] for g in genes.values()),
          '| genes with an interaction ratio', sum(g['ir'] is not None for g in genes.values()),
          '| %.0f KB' % (os.path.getsize(a.out) / 1024))


if __name__ == '__main__':
    main()
