"""Join archived GEO donor metadata to verified ENA runs, without methylation data."""
import csv, hashlib, json, re
from pathlib import Path

base=Path("Results/Task68_Metadata")
samples=[]
for block in (base/"GSE165915.soft").read_text().split("^SAMPLE = ")[1:]:
    def field(name):
        return re.findall(r"^!Sample_"+name+r" = (.+)$",block,re.M)
    characteristics=dict(x.split(": ",1) for x in field("characteristics_ch1"))
    relations=field("relation")
    samples.append(dict(gsm=block.splitlines()[0],donor=field("title")[0],
      biosample=next(re.search(r"SAMN\d+",x).group() for x in relations if "BioSample:" in x),
      experiment=next(re.search(r"SRX\d+",x).group() for x in relations if "SRA:" in x),
      exposure=characteristics["tertile"],original_group=characteristics["sample group"],
      tissue=characteristics["tissue"],batch=characteristics["batch"]))
runs=list(csv.DictReader((base/"PRJNA698569_files.tsv").open(),delimiter="\t"))
assert len(samples)==52 and len({s['donor'] for s in samples})==52
assert len(runs)==52 and all(r['library_layout']=='PAIRED' for r in runs)
lookup={r['sample_accession']:r for r in runs}; out=[]
for s in sorted(samples,key=lambda s:s['gsm']):
    r=lookup[s['biosample']]; assert r['experiment_accession']==s['experiment']
    urls=r['fastq_ftp'].split(';'); md5=r['fastq_md5'].split(';'); sizes=r['fastq_bytes'].split(';')
    assert len(urls)==len(md5)==len(sizes)==2 and s['tissue']=='sperm'
    assert all(re.fullmatch('[0-9a-f]{32}',x) for x in md5)
    assert urls[0].endswith('_1.fastq.gz') and urls[1].endswith('_2.fastq.gz')
    out.append(dict(index=len(out)+1,**s,run=r['run_accession'],read_pairs=r['read_count'],
       bytes_total=sum(map(int,sizes)),r1_url='https://'+urls[0],r2_url='https://'+urls[1],
       r1_md5=md5[0],r2_md5=md5[1]))
assert {e:sum(r['exposure']==e for r in out) for e in ('First','Third')}=={'First':26,'Third':26}
for name,rows in [('Task68_SpermSamples.tsv',out),('Task68_SpermPilotSamples.tsv',
                  [min((r for r in out if r['exposure']==e),key=lambda r:r['bytes_total']) for e in ('First','Third')])]:
    with (Path('Results')/name).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
manifest=dict(status='metadata-only cohort choice; not an external-analysis lock',donors=52,
  compressed_fastq_bytes=sum(r['bytes_total'] for r in out),
  sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [base/'GSE165915.soft',base/'PRJNA698569_files.tsv']})
Path('Results/Task68_SelectionManifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:v for k,v in manifest.items() if k!='sources'},indent=2))
