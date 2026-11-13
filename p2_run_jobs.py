"""Runs a docking manifest (resumable; skips jobs whose pose + json exist). Usage: PYTHONDONTWRITEBYTECODE=1 python p2_run_jobs.py manifest.csv [--workers N]
Workers use spawn; each job uses the cpu value in the manifest. Rebuilds docking_runs.csv at the end."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_lib as L
import multiprocessing as mp
if __name__ == '__main__':
    man = sys.argv[1]; workers = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 1
    log = L.Log('p2_run_jobs', f'manifest={man} workers={workers} cpu_override={sys.argv[sys.argv.index("--cpu")+1] if "--cpu" in sys.argv else None}')
    jobs = L.read_manifest(man)
    if '--cpu' in sys.argv:
        for j in jobs: j['cpu'] = int(sys.argv[sys.argv.index('--cpu') + 1])
    pending = [j for j in jobs if not (os.path.exists(L.pose_path(j)) and os.path.exists(L.meta_path(j)))]
    log.p(f'{len(jobs)} jobs in manifest, {len(jobs)-len(pending)} already done, {len(pending)} to run')
    t0 = time.time(); nerr = 0
    if pending:
        ctx = mp.get_context('spawn')
        with ctx.Pool(workers) as pool:
            for k, (st, jid, row) in enumerate(pool.imap_unordered(L.run_job_safe, pending), 1):
                if st == 'error': nerr += 1; log.p('ERROR', jid, row)
                elif k % 10 == 0 or k == len(pending): log.p(f'  {k}/{len(pending)} done, {time.time()-t0:.0f} s elapsed')
    n = L.rebuild_docking_runs_csv()
    log.p(f'finished {time.time()-t0:.0f} s, errors {nerr}, docking_runs.csv rows {n}')
