"""Single-donor pilot: download, align, deduplicate and run counts-only feasibility.
This does not calculate the external covariance endpoint or claim a protocol lock.
"""
import argparse,csv,fcntl,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path

def sha256(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def fetch(url,path,expected):
    def md5(p):
        with p.open('rb') as f:return hashlib.file_digest(f,'md5').hexdigest()
    if path.exists():
        if md5(path)!=expected:raise RuntimeError('existing FASTQ checksum mismatch: '+str(path))
        return
    tmp=path.with_suffix(path.suffix+'.partial')
    if tmp.exists() and md5(tmp)==expected:
        tmp.rename(path);return
    subprocess.run(['curl','--fail','--location','--retry','3','--silent','--show-error',
                    '--continue-at','-','--output',str(tmp),url],check=True)
    if md5(tmp)!=expected:raise RuntimeError('download checksum mismatch: '+str(tmp))
    tmp.rename(path)

def run_stage(cmd,outputs,out,number):
    """Resume identical completed stages; failed stages never get a success marker."""
    marker=out/f'step{number}.done.json'
    if marker.exists():
        prior=json.loads(marker.read_text())
        if prior['command']!=cmd or any(not p.is_file() or p.stat().st_size!=prior['output_bytes'].get(p.name) for p in outputs):
            raise RuntimeError(f'completed stage {number} changed; inspect before reuse')
        print(f'Reusing completed stage {number}',flush=True);return
    start=time.monotonic();wrapped=cmd
    if Path('/usr/bin/time').exists():
        flag='-l' if sys.platform=='darwin' else '-v'
        wrapped=['/usr/bin/time',flag,'-o',str(out/f'step{number}.resources.txt'),*cmd]
    print(f'Starting stage {number}: {cmd[0]}',flush=True)
    with (out/f'step{number}.log').open('a') as log:
        subprocess.run(wrapped,stdout=log,stderr=subprocess.STDOUT,check=True)
    if any(not p.is_file() or p.stat().st_size==0 for p in outputs):
        raise RuntimeError(f'stage {number} returned success without all outputs')
    record=dict(command=cmd,elapsed_seconds=time.monotonic()-start,
      output_bytes={p.name:p.stat().st_size for p in outputs},
      directory_bytes_after=sum(p.stat().st_size for p in out.rglob('*') if p.is_file()))
    marker.write_text(json.dumps(record,indent=2)+'\n')
    print(f'Completed stage {number} in {record["elapsed_seconds"]:.1f} seconds',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--row',type=int,required=True)
    p.add_argument('--work',required=True);p.add_argument('--reference-dir',required=True);p.add_argument('--map-dir',required=True)
    p.add_argument('--threads',type=int,default=8);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    rows=list(csv.DictReader(open(a.manifest),delimiter='\t'))
    if not 1<=a.row<=len(rows):raise ValueError('row outside manifest')
    s=rows[a.row-1]
    root=Path(__file__).resolve().parents[1];out=Path(a.work).resolve()/s['run']
    ref=Path(a.reference_dir).resolve();mp=Path(a.map_dir).resolve()
    commands=[
      ['trim_galore','--paired','--quality','20','--length','30','--clip_r1','8','--clip_r2','20',
       '--three_prime_clip_r1','8','--three_prime_clip_r2','8','--cores','2','--output_dir',str(out),
       str(out/'R1.fastq.gz'),str(out/'R2.fastq.gz')],
      ['bismark','--genome',str(ref),'--parallel','1','-p',str(max(1,a.threads//2)),
       '--maxins','1000','--output_dir',str(out),'-1',str(out/'R1_val_1.fq.gz'),'-2',str(out/'R2_val_2.fq.gz')],
      ['deduplicate_bismark','--paired','--bam','--output_dir',str(out),str(out/'R1_val_1_bismark_bt2_pe.bam')],
      ['samtools','sort','-n','-@',str(max(1,a.threads//2)),'-m','1G','-o',str(out/'deduplicated.names.bam'),
       str(out/'R1_val_1_bismark_bt2_pe.deduplicated.bam')],
      [os.sys.executable,str(root/'Scripts/SharedReadNoise.CoverageQC.py'),'--bam',str(out/'deduplicated.names.bam'),
       '--map-dir',str(mp),'--run',s['run'],'--prefix',str(out/'counts_only')]]
    if a.dry_run:
        print(json.dumps(dict(sample=s,commands=commands,status='counts-only pilot; not external test'),indent=2));return
    out.mkdir(parents=True,exist_ok=True)
    lock=(out/'.pilot.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('another process is already running this donor')
    lock.seek(0);lock.truncate();lock.write(str(os.getpid())+'\n');lock.flush()
    for t in ['curl','trim_galore','bismark','deduplicate_bismark','samtools','bowtie2']:
        if not shutil.which(t):raise RuntimeError('missing tool: '+t)
    if not (mp/'manifest.json').exists():raise RuntimeError('missing reference CpG map')
    fasta=ref/'hg38.fa'
    if not fasta.exists():raise RuntimeError('expected reference-dir/hg38.fa')
    reference_sha=sha256(fasta)
    if json.loads((mp/'manifest.json').read_text())['reference_sha256']!=reference_sha:
        raise RuntimeError('CpG map does not match reference FASTA')
    for strand in ['CT','GA']:
        folder=ref/'Bisulfite_Genome'/(strand+'_conversion')
        if not list(folder.glob('BS_'+strand+'.1.bt2*')):raise RuntimeError('missing Bismark index: '+str(folder))
    provenance=dict(sample=s,reference_sha256=reference_sha,map_manifest_sha256=sha256(mp/'manifest.json'),
      code_sha256={p.name:sha256(p) for p in [Path(__file__),root/'Scripts/SharedReadNoise.CoverageQC.py']})
    if out.exists() and (out/'pilot_commands.json').exists():
        if json.loads((out/'pilot_commands.json').read_text())!=commands:
            raise RuntimeError('started pilot commands changed; inspect before reuse')
        if not (out/'pilot_provenance.json').is_file() or json.loads((out/'pilot_provenance.json').read_text())!=provenance:
            raise RuntimeError('started pilot input or counting code changed; inspect before reuse')
        if (out/'PILOT_COMPLETE').exists():
            if not (out/'counts_only.coverage.npz').is_file():raise RuntimeError('completed coverage output is missing')
            prior=json.loads((out/'counts_only.qc.json').read_text())
            if prior['map_reference_sha256']!=reference_sha:raise RuntimeError('completed reference changed')
            print('Already completed identical counts-only pilot:',s['run']);return
        print('Resuming identical incomplete pilot:',s['run'],flush=True)
    out.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(out).free<max(100_000_000_000,4*int(s['bytes_total'])):raise RuntimeError('insufficient free scratch space')
    (out/'pilot_commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    (out/'pilot_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    for mate in [1,2]:fetch(s[f'r{mate}_url'],out/f'R{mate}.fastq.gz',s[f'r{mate}_md5'])
    outputs=[['R1_val_1.fq.gz','R2_val_2.fq.gz'],['R1_val_1_bismark_bt2_pe.bam'],
      ['R1_val_1_bismark_bt2_pe.deduplicated.bam'],['deduplicated.names.bam'],
      ['counts_only.coverage.npz','counts_only.qc.json']]
    for i,cmd in enumerate(commands):run_stage(cmd,[out/n for n in outputs[i]],out,i+1)
    (out/'PILOT_COMPLETE').write_text('Counts-only pilot complete; external covariance analysis not run.\n')

if __name__=='__main__':main()
