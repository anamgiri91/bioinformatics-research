"""Counts-only WGBS feasibility: whole-mate pairing, no methylation values retained.

Input BAM must be query-name sorted, deduplicated Bismark paired-end alignments.
CpG maps are generated from the exact reference FASTA by --prepare-map.
All output paths must be private: site-level coverage resolves one donor.
"""
import argparse, collections, gzip, hashlib, itertools, json, re
from pathlib import Path
import numpy as np
import pysam

def partition(run,qname,seed):
    key=f'{seed}\0{run}\0{qname}'.encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8],'big') % 3

def prepare_map(reference,dest):
    dest=Path(dest);dest.mkdir(parents=True,exist_ok=True)
    fasta=pysam.FastaFile(reference); entries=[];offset=0
    for chrom in ['chr%d'%i for i in range(1,23)]:
        seq=fasta.fetch(chrom).upper()
        pos=np.fromiter((m.start()+1 for m in re.finditer('CG',seq)),dtype='<u4')
        pos.tofile(dest/(chrom+'.u32'))
        entries.append(dict(chrom=chrom,offset=offset,cpgs=len(pos)));offset+=len(pos)
    with open(reference,'rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
    (dest/'manifest.json').write_text(json.dumps(dict(reference_sha256=sha,entries=entries),indent=2))

def indices(read,positions,min_quality):
    if not read.has_tag('XM') or not read.has_tag('XG'): raise ValueError('missing Bismark XM/XG tag')
    xg=read.get_tag('XG')
    if xg not in ('CT','GA'): raise ValueError('unsupported XG strand')
    xm=read.get_tag('XM');qual=read.query_qualities
    if len(xm)!=read.query_length or qual is None: raise ValueError('invalid XM length or absent qualities')
    q2r=dict(read.get_aligned_pairs(matches_only=True)); result=[]
    for q,c in enumerate(xm):
        # Deliberately ignore case: z and Z both mean coverage, not a state.
        if c not in 'zZ' or qual[q]<min_quality: continue
        rp=q2r.get(q)
        if rp is None: continue
        pos=rp+(1 if xg=='CT' else 0)
        i=int(np.searchsorted(positions,pos))
        if i>=len(positions) or positions[i]!=pos: raise ValueError('CpG call off reference map')
        result.append(i)
    return result

def coverage(bam_path,map_dir,run,prefix,seed=2026100268,min_quality=20):
    manifest=json.loads((Path(map_dir)/'manifest.json').read_text());entries=manifest['entries']
    maps={e['chrom']:np.fromfile(Path(map_dir)/(e['chrom']+'.u32'),dtype='<u4') for e in entries}
    offsets={e['chrom']:e['offset'] for e in entries};total=sum(e['cpgs'] for e in entries)
    N=np.zeros((3,total),dtype=np.uint32);rA=np.zeros(total,dtype=np.uint32)
    st=collections.Counter();inserts=collections.Counter()
    with pysam.AlignmentFile(bam_path,'rb') as bam:
        if bam.header.to_dict().get('HD',{}).get('SO')!='queryname':
            raise ValueError('require samtools sort -n and SO:queryname header')
        for name,group in itertools.groupby(bam,key=lambda r:r.query_name):
            reads=list(group);st['query_groups']+=1
            reads=[r for r in reads if not (r.is_secondary or r.is_supplementary)]
            if len(reads)!=2 or sum(r.is_read1 for r in reads)!=1 or sum(r.is_read2 for r in reads)!=1:
                st['incomplete_or_ambiguous_pair']+=1;continue
            if any(r.is_unmapped or r.is_qcfail or r.is_duplicate or not r.is_proper_pair for r in reads):
                st['filtered_pair']+=1;continue
            a,b=reads;chrom=a.reference_name
            if chrom!=b.reference_name or chrom not in maps: st['off_autosome_pair']+=1;continue
            calls=sorted(set(indices(a,maps[chrom],min_quality))|set(indices(b,maps[chrom],min_quality)))
            st['complete_autosomal_fragments']+=1
            inserts[min(abs(a.template_length)//50*50,2000)]+=1
            if not calls: st['no_usable_cpg']+=1;continue
            ix=np.asarray(calls,dtype=np.int64);global_ix=ix+offsets[chrom]
            part=partition(run,name,seed);N[part,global_ix]+=1
            st['fragments_'+str(part)]+=1
            if part==0 and len(ix)>1:
                adjacent=np.diff(ix)==1
                gap=np.diff(maps[chrom][ix].astype(np.int64))
                rA[global_ix[:-1][adjacent & (gap<=200)]]+=1
    valid=np.zeros(total,dtype=bool)
    for e in entries:
        p=maps[e['chrom']];o=e['offset'];left=np.arange(o,o+len(p)-1)
        valid[left]=np.diff(p)<=200
    eligible=np.zeros(total,dtype=bool);left=np.flatnonzero(valid)
    eligible[left]=np.all(N[:,left]>=2,axis=0)&np.all(N[:,left+1]>=2,axis=0)&((rA[left]==0)|(rA[left]>=2))
    prefix=Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(str(prefix)+'.coverage.npz',N=N,rA=rA,eligible=np.packbits(eligible))
    result=dict(run=run,seed=seed,min_quality=min_quality,counts=dict(st),
      map_reference_sha256=manifest['reference_sha256'],cpgs=total,
      covered_cpgs=int(np.count_nonzero(N.sum(axis=0))),eligible_pairs=int(eligible.sum()),
      insert_bins_50bp=dict(sorted(inserts.items())),
      scope='counts only; no beta/covariance/ranking calculated; no state-dependent consensus filter')
    Path(str(prefix)+'.qc.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepare-map',metavar='FASTA');p.add_argument('--map-dir',required=True)
    p.add_argument('--bam');p.add_argument('--run');p.add_argument('--prefix');args=p.parse_args()
    if args.prepare_map: prepare_map(args.prepare_map,args.map_dir)
    else:
        if not all([args.bam,args.run,args.prefix]):p.error('--bam, --run and --prefix required')
        coverage(args.bam,args.map_dir,args.run,args.prefix)
