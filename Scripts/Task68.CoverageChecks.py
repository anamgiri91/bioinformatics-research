"""Independent synthetic alignment checks for counts-only cluster feasibility."""
import importlib.util, tempfile, json
from pathlib import Path
import numpy as np
import pysam
spec=importlib.util.spec_from_file_location('qc',Path(__file__).with_name('SharedReadNoise.CoverageQC.py'))
qc=importlib.util.module_from_spec(spec);spec.loader.exec_module(qc)
with tempfile.TemporaryDirectory() as td:
    d=Path(td);mp=d/'map';mp.mkdir()
    np.array([3,9],dtype='<u4').tofile(mp/'chr1.u32')
    (mp/'manifest.json').write_text(json.dumps(dict(reference_sha256='synthetic',entries=[dict(chrom='chr1',offset=0,cpgs=2)])))
    header={'HD':{'VN':'1.6','SO':'queryname'},'SQ':[{'SN':'chr1','LN':100}]}
    names={p:[] for p in range(3)}
    for i in range(200):
        name=f'f{i:04d}';p=qc.partition('test',name,2026100268)
        if len(names[p])<2:names[p].append(name)
    assert all(len(v)==2 for v in names.values())
    def read(name,flag,xg='CT',soft=False):
        r=pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header));r.query_name=name;r.flag=flag
        r.reference_id=0;r.reference_start=2 if xg=='CT' else 3;r.mapping_quality=40
        seq='CGAAAACGAA';xm='Z.....z...'
        if soft:seq='AA'+seq;xm='..'+xm;r.cigartuples=[(4,2),(0,10)]
        else:r.cigartuples=[(0,10)]
        r.query_sequence=seq;r.query_qualities=pysam.qualitystring_to_array('I'*len(seq))
        r.next_reference_id=0;r.next_reference_start=r.reference_start;r.template_length=10
        r.set_tag('XM',xm);r.set_tag('XG',xg);return r
    r=read('cigar',99,soft=True);assert qc.indices(r,np.array([3,9]),20)==[0,1]
    r=read('bottom',99,xg='GA');assert qc.indices(r,np.array([3,9]),20)==[0,1]
    r.query_qualities=pysam.qualitystring_to_array('!'*10);assert qc.indices(r,np.array([3,9]),20)==[]
    r=read('state',99);a=qc.indices(r,np.array([3,9]),20);r.set_tag('XM',r.get_tag('XM').swapcase())
    assert qc.indices(r,np.array([3,9]),20)==a
    bam=d/'paired.bam'
    with pysam.AlignmentFile(bam,'wb',header=header) as h:
        for name in sorted(sum(names.values(),[])):
            h.write(read(name,99));h.write(read(name,147))
        h.write(read('orphan',99))
    qc.coverage(str(bam),str(mp),'test',str(d/'out'))
    c=np.load(str(d/'out.coverage.npz'));report=json.loads((d/'out.qc.json').read_text())
    assert np.array_equal(c['N'],np.full((3,2),2)) # mates count once, never twice
    assert list(c['rA'])==[2,0] and report['eligible_pairs']==1
    assert report['counts']['incomplete_or_ambiguous_pair']==1
    assert qc.partition('test','same',1)==qc.partition('test','same',1)
    def single_site(name,flag,second=False):
        r=read(name,flag);r.set_tag('XM','......z...' if second else 'Z.........');return r
    # The mates need not share any called CpG: the physical fragment still links them.
    bam=d/'nonoverlap.bam'
    with pysam.AlignmentFile(bam,'wb',header=header) as h:
        for name in sorted(sum(names.values(),[])):
            h.write(single_site(name,99));h.write(single_site(name,147,second=True))
    qc.coverage(str(bam),str(mp),'test',str(d/'nonoverlap'))
    with np.load(d/'nonoverlap.coverage.npz') as c:
        assert np.array_equal(c['N'],np.full((3,2),2)) and list(c['rA'])==[2,0]
        assert np.unpackbits(c['eligible'])[0]==1
    # Adequate site coverage does not make r_A=1 supported: replace one A
    # shared fragment with two different fragments, one covering each site.
    extras=[]
    for i in range(200,500):
        name=f'f{i:04d}'
        if qc.partition('test',name,2026100268)==0:extras.append(name)
        if len(extras)==2:break
    groups=[]
    for name in sum(names.values(),[]):
        if name!=names[0][1]:groups.append((name,[single_site(name,99),single_site(name,147,True)]))
    for second,name in enumerate(extras):
        groups.append((name,[single_site(name,99,bool(second)),single_site(name,147,bool(second))]))
    bam=d/'single_shared.bam'
    with pysam.AlignmentFile(bam,'wb',header=header) as h:
        for name,reads in sorted(groups):
            for r in reads:h.write(r)
    qc.coverage(str(bam),str(mp),'test',str(d/'single_shared'))
    with np.load(d/'single_shared.coverage.npz') as c:
        assert np.array_equal(c['N'],np.full((3,2),2)) and list(c['rA'])==[1,0]
        assert np.unpackbits(c['eligible'])[0]==0
print('PASS: mate union, non-overlapping mate linkage, shared-fragment counts, r_A=1 exclusion, strand coordinates, soft clipping, base quality, state-blind coverage, orphan reporting and deterministic partition')
