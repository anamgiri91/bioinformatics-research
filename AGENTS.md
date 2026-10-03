# Research session continuity

Before continuing the shared-fragment CpG project, read
`SESSION_HANDOFF.md`, then `plan_shared_read_noise.md` and the linked task
record. The outlier-analysis `plan.md` is a different roadmap.

Current user constraints:

- Do not use the cluster, SSH or SLURM. Work locally.
- Leave the plasma cohort redo for later.
- The user requested shutdown at the end of this session. All Task 68
  jobs were stopped. Resume computation only when the user asks to continue;
  do not restart merely because this file is read.
- Preserve existing frozen analyses and unrelated working-tree changes.
  Individual-level data, alignments and coverage arrays stay under ignored
  `Data/`; aggregate reports and code can be versioned.

Update the handoff when the run state changes. Never reuse historical PIDs.
