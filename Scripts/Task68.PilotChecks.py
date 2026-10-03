"""Check restart integrity with real child commands; profiler is omitted in tests."""
import importlib.util,subprocess,sys,tempfile
from pathlib import Path
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('pilot',Path(__file__).with_name('Task68.SpermPilot.py'))
pilot=importlib.util.module_from_spec(spec);spec.loader.exec_module(pilot)
real_run=subprocess.run
def without_profiler(cmd,**kw):
    if cmd[0]=='/usr/bin/time':cmd=cmd[4:]
    return real_run(cmd,**kw)

with tempfile.TemporaryDirectory(prefix='task68_restart_') as name,patch.object(pilot.subprocess,'run',side_effect=without_profiler):
    d=Path(name);output=d/'product';counter=d/'invocations'
    cmd=[sys.executable,'-c',
         'from pathlib import Path;import sys;Path(sys.argv[1]).write_text("ok");p=Path(sys.argv[2]);p.write_text((p.read_text() if p.exists() else "")+"x")',str(output),str(counter)]
    pilot.run_stage(cmd,[output],d,1)
    pilot.run_stage(cmd,[output],d,1)
    assert counter.read_text()=='x','completed stage executed twice'
    output.write_text('changed product')
    try:pilot.run_stage(cmd,[output],d,1)
    except RuntimeError:pass
    else:raise AssertionError('changed product was silently reused')
    try:pilot.run_stage([sys.executable,'-c','raise SystemExit(7)'],[d/'absent'],d,2)
    except subprocess.CalledProcessError as e:assert e.returncode==7
    else:raise AssertionError('failed subprocess accepted')
    assert not (d/'step2.done.json').exists(),'failed stage received success marker'
    try:pilot.run_stage([sys.executable,'-c','pass'],[d/'absent'],d,3)
    except RuntimeError:pass
    else:raise AssertionError('missing output accepted')
    assert not (d/'step3.done.json').exists()
print('PASS: completed stage reused once, changed output rejected, failed/missing-output stages never marked complete; resource profiler excluded from this test.')
