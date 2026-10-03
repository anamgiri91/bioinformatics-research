"""Bounded, resumable ranged download of one pilot donor; verify ENA MD5s.
Uses four HTTPS connections, keeps completed chunks, and never runs alignments.
"""
import argparse,concurrent.futures,csv,hashlib,json,shutil,subprocess,time
from pathlib import Path

CHUNK=128*1024*1024

def md5(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'md5').hexdigest()

def fetch_chunk(job):
    url,total,start,end,dest=job;size=end-start+1
    if dest.exists() and dest.stat().st_size==size:return
    temp=dest.with_suffix('.partial')
    result=subprocess.run(['curl','--fail','--location','--retry','3','--silent','--show-error',
      '--max-time','1800','--range',f'{start}-{end}','--max-filesize',str(size),
      '--output',str(temp),'--write-out','%{http_code}',url],capture_output=True,text=True,check=True)
    if result.stdout.strip()!='206' or temp.stat().st_size!=size:
        raise RuntimeError('server did not return the requested byte range')
    temp.rename(dest)

def main():
    p=argparse.ArgumentParser();p.add_argument('--row',type=int,required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    samples=list(csv.DictReader((root/'Results/Task68_SpermPilotSamples.tsv').open(),delimiter='\t'))
    if a.row not in (1,2):raise ValueError('pilot row must be 1 or 2')
    s=samples[a.row-1]
    ena=next(r for r in csv.DictReader((root/'Results/Task68_Metadata/PRJNA698569_files.tsv').open(),delimiter='\t')
             if r['run_accession']==s['run'])
    sizes=list(map(int,ena['fastq_bytes'].split(';')))
    out=root/'Data/gse165915_sperm_wgbs'/s['run'];out.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(out).free<100_000_000_000:raise RuntimeError('less than 100 GB free')
    jobs=[];parts={};start_time=time.monotonic()
    for mate,total in enumerate(sizes,1):
        final=out/f'R{mate}.fastq.gz'
        if final.exists():
            if md5(final)!=s[f'r{mate}_md5']:raise RuntimeError('existing FASTQ checksum mismatch')
            continue
        d=out/f'R{mate}.download_parts';d.mkdir(exist_ok=True)
        parts[mate]=[]
        for start in range(0,total,CHUNK):
            end=min(total-1,start+CHUNK-1);dest=d/f'{start:012d}-{end:012d}.chunk'
            jobs.append((s[f'r{mate}_url'],total,start,end,dest));parts[mate].append(dest)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(fetch_chunk,j) for j in jobs]
        for n,f in enumerate(concurrent.futures.as_completed(futures),1):
            f.result()
            print(s['run'],n,'/',len(jobs),'chunks complete;',round(time.monotonic()-start_time),'seconds',flush=True)
    for mate,files in parts.items():
        temp=out/f'R{mate}.assembled.partial';final=out/f'R{mate}.fastq.gz'
        with temp.open('wb') as dst:
            for f in files:
                with f.open('rb') as src:shutil.copyfileobj(src,dst,8*1024*1024)
        if md5(temp)!=s[f'r{mate}_md5']:raise RuntimeError('assembled FASTQ checksum mismatch; chunks retained for inspection')
        temp.rename(final)
    result=dict(run=s['run'],elapsed_seconds=time.monotonic()-start_time,connections=4,
      bytes=sum(sizes),md5={f'R{m}.fastq.gz':md5(out/f'R{m}.fastq.gz') for m in (1,2)},
      note='Completed download chunks retained; earlier sequential partial file, if present, is unused.')
    (out/'DOWNLOAD_COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Verified both FASTQs:',s['run'],flush=True)

if __name__=='__main__':main()
