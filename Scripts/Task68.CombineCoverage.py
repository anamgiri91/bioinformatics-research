"""Aggregate only fragment-coverage eligibility, never methylation outcomes."""
import argparse,csv,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--work',required=True)
p.add_argument('--output',required=True);a=p.parse_args()
rows=list(csv.DictReader(open(a.manifest),delimiter='\t'));counts=None;reference=None;missing=[]
for r in rows:
    prefix=Path(a.work)/r['run']/'counts_only'
    jf=Path(str(prefix)+'.qc.json');nf=Path(str(prefix)+'.coverage.npz')
    if not jf.exists() or not nf.exists():missing.append(r['run']);continue
    qc=json.loads(jf.read_text());assert qc['run']==r['run']
    if reference is None:reference=qc['map_reference_sha256'];counts=np.zeros(qc['cpgs'],dtype=np.uint16)
    assert reference==qc['map_reference_sha256'] and len(counts)==qc['cpgs']
    with np.load(nf) as d:counts+=np.unpackbits(d['eligible'])[:len(counts)]
result=dict(expected_donors=len(rows),completed_donors=len(rows)-len(missing),missing_runs=missing,
  eligible_pairs_at_20_donors=int(np.sum(counts>=20)) if counts is not None else 0,
  reference_sha256=reference,status='complete counts-only screen' if not missing else 'incomplete counts-only screen',
  caveat='Coverage precedes methylation-state conflict masking; final eligible count may be lower.')
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
