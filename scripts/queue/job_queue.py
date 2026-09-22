"""
OXTR pipeline job queue - SQLite-backed, resumable, GPU/CPU-aware orchestration
for the scale-up run. Built for v2.0.0 to replace the ad hoc nohup/backgrounded
shell commands used throughout the 100-backbone pilot, which don't hold up at
~40,000-candidate scale (no crash recovery, no automatic resource assignment,
no centralized status).

Design: this queue orchestrates WHEN and WHERE work runs and tracks WHETHER it
succeeded - it does not reimplement pipeline science. Each job's "command" is a
shell command string (typically invoking one of this project's existing stage
scripts) that the worker pool executes and monitors.

Usage:
  python job_queue.py init                          # create the DB
  python job_queue.py enqueue <stage> <ids_file> <cmd_template> [--resources gpu|cpu]
  python job_queue.py worker --resource gpu:0        # run one worker bound to GPU 0
  python job_queue.py worker --resource cpu --slots 8
  python job_queue.py status
  python job_queue.py requeue-failed [--stage STAGE]
"""
import argparse
import sqlite3
import subprocess
import sys
import time
import os
from pathlib import Path
from datetime import datetime, timezone

DB_PATH = Path(__file__).parent / "queue.db"
LOG_DIR = Path(__file__).parent / "logs"


def get_conn():
    conn = sqlite3.connect(DB_PATH, timeout=30)
    conn.execute("PRAGMA journal_mode=WAL")  # allows concurrent workers safely
    return conn


def cmd_init(args):
    LOG_DIR.mkdir(exist_ok=True)
    conn = get_conn()
    conn.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            stage TEXT NOT NULL,
            candidate_id TEXT NOT NULL,
            command TEXT NOT NULL,
            resource_type TEXT NOT NULL DEFAULT 'cpu',  -- 'gpu' or 'cpu'
            status TEXT NOT NULL DEFAULT 'pending',      -- pending/running/done/failed
            assigned_resource TEXT,
            started_at TEXT,
            finished_at TEXT,
            exit_code INTEGER,
            log_path TEXT,
            attempt INTEGER NOT NULL DEFAULT 0,
            UNIQUE(stage, candidate_id)
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_status_resource ON jobs(status, resource_type)")
    conn.commit()
    print(f"Initialized queue DB at {DB_PATH}")


def cmd_enqueue(args):
    conn = get_conn()
    ids = [l.strip() for l in open(args.ids_file) if l.strip()]
    n = 0
    for cid in ids:
        command = args.cmd_template.replace("{id}", cid)
        try:
            conn.execute(
                "INSERT INTO jobs (stage, candidate_id, command, resource_type) VALUES (?,?,?,?)",
                (args.stage, cid, command, args.resources),
            )
            n += 1
        except sqlite3.IntegrityError:
            pass  # already enqueued (stage, candidate_id) - resumability
    conn.commit()
    print(f"Enqueued {n} new jobs for stage '{args.stage}' ({len(ids) - n} already existed)")


def _claim_job(conn, resource_type, resource_label, stage_filter=None):
    """Atomically claim one pending job for this resource type.

    stage_filter restricts claiming to one named stage - use it whenever the
    queue holds jobs from more than one stage at once (e.g. a large staged-
    but-not-yet-authorized run sitting alongside smaller ad hoc jobs).
    Without it, workers claim strictly by lowest job id regardless of stage,
    which silently launches whatever was enqueued first - a real footgun if
    a big job is staged (queued, not yet approved to run) and a worker for
    something else gets started before it's meant to.
    """
    conn.execute("BEGIN IMMEDIATE")
    if stage_filter:
        row = conn.execute(
            "SELECT id, stage, candidate_id, command FROM jobs "
            "WHERE status='pending' AND resource_type=? AND stage=? ORDER BY id LIMIT 1",
            (resource_type, stage_filter),
        ).fetchone()
    else:
        row = conn.execute(
            "SELECT id, stage, candidate_id, command FROM jobs "
            "WHERE status='pending' AND resource_type=? ORDER BY id LIMIT 1",
            (resource_type,),
        ).fetchone()
    if row is None:
        conn.execute("COMMIT")
        return None
    job_id = row[0]
    now = datetime.now(timezone.utc).isoformat()
    conn.execute(
        "UPDATE jobs SET status='running', assigned_resource=?, started_at=?, attempt=attempt+1 WHERE id=?",
        (resource_label, now, job_id),
    )
    conn.execute("COMMIT")
    return row


def cmd_worker(args):
    resource_type, _, label = args.resource.partition(":")
    label = args.resource
    LOG_DIR.mkdir(exist_ok=True)
    conn = get_conn()
    print(f"Worker started, bound to resource '{args.resource}' (type={resource_type})")
    idle_count = 0
    while True:
        row = _claim_job(conn, resource_type, label, stage_filter=args.stage)
        if row is None:
            idle_count += 1
            if idle_count >= args.idle_exit:
                print(f"No pending jobs for {idle_count} checks, exiting.")
                break
            time.sleep(args.poll_interval)
            continue
        idle_count = 0
        job_id, stage, candidate_id, command = row
        log_path = LOG_DIR / f"{stage}_{candidate_id}.log"
        env = os.environ.copy()
        if resource_type == "gpu":
            gpu_id = label.split(":")[1] if ":" in label else label
            env["CUDA_VISIBLE_DEVICES"] = gpu_id
        print(f"[{label}] Running {stage}/{candidate_id}")
        with open(log_path, "w") as logf:
            proc = subprocess.run(command, shell=True, env=env, stdout=logf, stderr=subprocess.STDOUT)
        now = datetime.now(timezone.utc).isoformat()
        status = "done" if proc.returncode == 0 else "failed"
        conn.execute(
            "UPDATE jobs SET status=?, finished_at=?, exit_code=?, log_path=? WHERE id=?",
            (status, now, proc.returncode, str(log_path), job_id),
        )
        conn.commit()
        print(f"[{label}] {stage}/{candidate_id}: {status} (exit {proc.returncode})")


def cmd_status(args):
    conn = get_conn()
    rows = conn.execute(
        "SELECT stage, status, COUNT(*) FROM jobs GROUP BY stage, status ORDER BY stage, status"
    ).fetchall()
    if not rows:
        print("No jobs in queue.")
        return
    stages = {}
    for stage, status, count in rows:
        stages.setdefault(stage, {})[status] = count
    print(f"{'Stage':<25} {'pending':>8} {'running':>8} {'done':>8} {'failed':>8}")
    for stage, counts in stages.items():
        print(f"{stage:<25} {counts.get('pending',0):>8} {counts.get('running',0):>8} "
              f"{counts.get('done',0):>8} {counts.get('failed',0):>8}")


def cmd_requeue_failed(args):
    conn = get_conn()
    q = "UPDATE jobs SET status='pending' WHERE status='failed'"
    params = ()
    if args.stage:
        q += " AND stage=?"
        params = (args.stage,)
    cur = conn.execute(q, params)
    conn.commit()
    print(f"Requeued {cur.rowcount} failed jobs" + (f" (stage={args.stage})" if args.stage else ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("init")

    p_enq = sub.add_parser("enqueue")
    p_enq.add_argument("stage")
    p_enq.add_argument("ids_file", help="text file, one candidate id per line")
    p_enq.add_argument("cmd_template", help="shell command with {id} placeholder")
    p_enq.add_argument("--resources", choices=["gpu", "cpu"], default="cpu")

    p_work = sub.add_parser("worker")
    p_work.add_argument("--resource", required=True, help="'gpu:0', 'gpu:1', or 'cpu'")
    p_work.add_argument("--stage", default=None,
                         help="restrict this worker to one stage - REQUIRED whenever the queue "
                              "might hold jobs from more than one stage, to avoid accidentally "
                              "launching a large staged-but-not-yet-authorized job by lowest id")
    p_work.add_argument("--poll-interval", type=float, default=5.0)
    p_work.add_argument("--idle-exit", type=int, default=6, help="exit after N consecutive empty polls")

    sub.add_parser("status")

    p_req = sub.add_parser("requeue-failed")
    p_req.add_argument("--stage", default=None)

    args = parser.parse_args()
    {
        "init": cmd_init,
        "enqueue": cmd_enqueue,
        "worker": cmd_worker,
        "status": cmd_status,
        "requeue-failed": cmd_requeue_failed,
    }[args.command](args)


if __name__ == "__main__":
    main()
