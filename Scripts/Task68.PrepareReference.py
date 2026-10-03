"""Prepare the pinned full UCSC hg38 reference locally, with provenance.
Download hg38.fa.gz.partial and md5sum.txt first; this script uses no network.
"""
import gzip,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'Data/gse165915_sperm_reference'
REF=BASE/'hg38_bismark';MAP=BASE/'hg38_cpg_map'
EXPECTED_MD5='1c9dcaddfa41027f17cd8f7a82c7293b'

def digest(p,algorithm='sha256'):
    with p.open('rb') as f:return hashlib.file_digest(f,algorithm).hexdigest()

def run(cmd,name):
    start=time.monotonic();print('Starting',name,flush=True)
    clock=['/usr/bin/time','-l' if sys.platform=='darwin' else '-v','-o',str(BASE/(name+'.resources.txt'))]
    with (BASE/(name+'.log')).open('a') as f:subprocess.run(clock+cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    print('Completed',name,'in',round(time.monotonic()-start,1),'seconds',flush=True)

def main():
    BASE.mkdir(parents=True,exist_ok=True);REF.mkdir(exist_ok=True)
    tool=ROOT/'Data/shared_fragment_tools/Bismark-0.24.2/bismark_genome_preparation'
    for t in ['samtools','bowtie2-build']:
        if not shutil.which(t):raise RuntimeError('missing tool '+t)
    listed=[x.split()[0] for x in (BASE/'md5sum.txt').read_text().splitlines() if x.split()[-1]=='hg38.fa.gz']
    if listed!=[EXPECTED_MD5]:raise RuntimeError('reference source checksum changed')
    archive=BASE/'hg38.fa.gz'
    candidate=archive if archive.exists() else BASE/'hg38.fa.gz.partial'
    if digest(candidate,'md5')!=EXPECTED_MD5:raise RuntimeError('reference download is incomplete or corrupt')
    if candidate!=archive:candidate.rename(archive)
    fasta=REF/'hg38.fa';done=BASE/'reference_fasta.json'
    if not done.exists():
        if fasta.exists():raise RuntimeError('unrecorded reference FASTA already exists')
        temp=REF/'hg38.fa.partial'
        with gzip.open(archive,'rb') as src,temp.open('wb') as dst:shutil.copyfileobj(src,dst,8*1024*1024)
        temp.rename(fasta)
        done.write_text(json.dumps(dict(source='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/initial/hg38.fa.gz',
          compressed_md5=EXPECTED_MD5,reference_sha256=digest(fasta)),indent=2)+'\n')
    info=json.loads(done.read_text())
    if digest(fasta)!=info['reference_sha256']:raise RuntimeError('reference FASTA changed')
    if not Path(str(fasta)+'.fai').exists():run(['samtools','faidx',str(fasta)],'faidx')
    if not (MAP/'manifest.json').exists():
        run([sys.executable,str(ROOT/'Scripts/SharedReadNoise.CoverageQC.py'),'--prepare-map',str(fasta),'--map-dir',str(MAP)],'cpg_map')
    if json.loads((MAP/'manifest.json').read_text())['reference_sha256']!=info['reference_sha256']:
        raise RuntimeError('CpG map reference mismatch')
    marker=BASE/'REFERENCE_COMPLETE.json'
    if not marker.exists():
        # Bismark creates CT and GA indices concurrently: one thread each
        # keeps the indexing load bounded on this 24-GiB machine.
        # Omit --parallel: Bismark accepts only >=2 for that option;
        # its default is one thread per indexer.
        run([str(tool),'--bowtie2',str(REF)],'bismark_genome_preparation')
        indices=[]
        for strand in ['CT','GA']:
            folder=REF/'Bisulfite_Genome'/(strand+'_conversion')
            for suffix in ['1','2','3','4','rev.1','rev.2']:
                found=list(folder.glob('BS_'+strand+'.'+suffix+'.bt2*'))
                if len(found)!=1 or found[0].stat().st_size==0:raise RuntimeError('incomplete index '+strand+' '+suffix)
                indices+=found
        info['index_sha256']={str(p.relative_to(BASE)):digest(p) for p in indices}
        info['preparation_script_sha256']=digest(tool)
        info['prepare_driver_sha256']=digest(Path(__file__))
        info['versions']={t:subprocess.check_output([t,'--version'],text=True,errors='replace').splitlines()[0]
                          for t in ['samtools','bowtie2-build']}
        marker.write_text(json.dumps(info,indent=2)+'\n')
    else:
        for name,expected in json.loads(marker.read_text())['index_sha256'].items():
            if digest(BASE/name)!=expected:raise RuntimeError('completed index changed: '+name)
    print('Full hg38 reference, CpG map and both Bismark indices are ready.',flush=True)

if __name__=='__main__':main()
