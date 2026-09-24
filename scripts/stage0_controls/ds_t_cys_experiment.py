"""
Cys-constrained, receptor-aware D_s(T) re-run.

Differs from the original ds_t_experiment.py in TWO ways, both deliberate and
both chosen to match what Stage 2 actually runs in production:
  1. --fixed_positions_jsonl pins the two motif cysteines, so every designed
     sequence can form the disulfide. The original produced 0/2400 cyclizable.
  2. Receptor-aware (--pdb_path_chains L on the full RFdiffusion complex)
     rather than binder-only.

CONSEQUENCE, stated up front: because two variables changed, a difference
against the original published D_s table CANNOT be attributed to the cysteine
constraint alone. This run measures production's D_s; it is not a controlled
comparison against the original.

32 backbones x 7 temperatures x 300 sequences. No new RFdiffusion cost.
"""
import subprocess, os, sys, json

SHARD, NSHARD, GPU = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
D = "/scratch/drewdog/denovo_binder_100_pilot/stage_1_backbones/disulfide_100"
OUT = "/scratch/drewdog/denovo_binder_100_pilot/stage_0_1_benchmark/ds_t_experiment_cys"
os.makedirs(OUT, exist_ok=True)

BACKBONES = json.load(open("/home/drewdog/.claude/jobs/e2254ab9/tmp/ds_backbones.json"))
TEMPS = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
N_SEQ = 300

jobs = [(b, t) for b in BACKBONES for t in TEMPS][SHARD::NSHARD]
print(f"[shard {SHARD}] {len(jobs)} jobs on GPU {GPU}", flush=True)

env = dict(os.environ, CUDA_VISIBLE_DEVICES=GPU)
done = fail = 0
for i, (b, t) in enumerate(jobs, 1):
    folder = f"{OUT}/b{b}_T{t}"
    if os.path.exists(f"{folder}/seqs/out_{b}.fa"):
        done += 1
        continue
    os.makedirs(folder, exist_ok=True)
    cmd = [
        "/scratch/drewdog/rfdiffusion/env_rfd/bin/python", "/scratch/drewdog/ProteinMPNN/protein_mpnn_run.py",
        "--pdb_path", f"{D}/run/out/out_{b}.pdb",
        "--pdb_path_chains", "L",
        "--out_folder", folder,
        "--num_seq_per_target", str(N_SEQ),
        "--sampling_temp", str(t),
        "--batch_size", "50",
        "--fixed_positions_jsonl", f"{D}/mpnn_out_v2/fixed_out_{b}.jsonl",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=900, env=env)
        done += 1
        if i % 10 == 0:
            print(f"[shard {SHARD}] {i}/{len(jobs)} ok={done} fail={fail}", flush=True)
    except Exception as e:
        fail += 1
        print(f"[shard {SHARD}] FAILED b={b} T={t}: {str(e)[:200]}", flush=True)
print(f"[shard {SHARD}] DONE ok={done} fail={fail}", flush=True)
