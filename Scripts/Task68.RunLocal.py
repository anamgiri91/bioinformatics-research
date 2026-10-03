"""Wait for this donor's local prerequisites, then run the first pilot to QC.
Bounded to 24 hours of prerequisite waiting. No SSH, SLURM or second-donor launch.
"""
import argparse,csv,datetime,json,os,subprocess,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def alive(pid):
    try:os.kill(pid,0);return True
    except ProcessLookupError:return False

def main():
    p=argparse.ArgumentParser();p.add_argument('--row',type=int,required=True)
    p.add_argument('--reference-pid',type=int,required=True);p.add_argument('--download-pid',type=int,required=True)
    a=p.parse_args()
    if a.row not in (1,2):raise ValueError('pilot row must be 1 or 2')
    s=list(csv.DictReader((ROOT/'Results/Task68_SpermPilotSamples.tsv').open(),delimiter='\t'))[a.row-1]
    out=ROOT/'Data/gse165915_sperm_wgbs'/s['run'];out.mkdir(parents=True,exist_ok=True)
    status=out/'LOCAL_PIPELINE_STATUS.json'
    def record(phase,**kwargs):
        status.write_text(json.dumps(dict(run=s['run'],pid=os.getpid(),phase=phase,
          updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kwargs),indent=2)+'\n')
        print(phase,kwargs,flush=True)
    try:
        needed=[(ROOT/'Data/gse165915_sperm_reference/REFERENCE_COMPLETE.json',a.reference_pid),
                (out/'DOWNLOAD_COMPLETE.json',a.download_pid)]
        deadline=time.monotonic()+24*3600
        while True:
            missing=[(f,pid) for f,pid in needed if not f.exists()]
            if not missing:break
            if any(not alive(pid) for f,pid in missing):raise RuntimeError('prerequisite process stopped without a completion record')
            if time.monotonic()>deadline:raise TimeoutError('local prerequisite wait exceeded 24 hours')
            record('waiting_for_local_prerequisites',missing=[f.name for f,pid in missing])
            time.sleep(30)
        record('running_counts_only_pilot')
        subprocess.run(['bash',str(ROOT/'Scripts/Task68.LocalPilot.sh'),str(a.row)],cwd=ROOT,check=True)
        qc=json.loads((out/'counts_only.qc.json').read_text())
        steps=[json.loads((out/f'step{k}.done.json').read_text()) for k in range(1,6)]
        summary=dict(run=s['run'],status='counts-only pilot completed; QC review required before next donor',
          counts=qc['counts'],covered_cpgs=qc['covered_cpgs'],eligible_pairs=qc['eligible_pairs'],
          insert_bins_50bp=qc['insert_bins_50bp'],reference_sha256=qc['map_reference_sha256'],
          stage_elapsed_seconds=[x['elapsed_seconds'] for x in steps],
          directory_bytes_after=[x['directory_bytes_after'] for x in steps],
          scope=qc['scope'],covariance_outcomes_computed=False)
        (ROOT/f'Results/Task68_LocalPilot_{s["run"]}.json').write_text(json.dumps(summary,indent=2)+'\n')
        record('completed_pending_qc_review',eligible_pairs=qc['eligible_pairs'])
    except Exception as e:
        record('failed',error=str(e));raise

if __name__=='__main__':main()
