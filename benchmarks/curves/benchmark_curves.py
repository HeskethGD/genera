#!/usr/bin/env python3
"""Serial, process-isolated genera/Sage first-kind period benchmarks."""
import argparse
import datetime as dt
import hashlib
import fcntl
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import sys
import tempfile
import time

from curve_cases import catalog, SMOKE, QUICK
from validation import compare

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORKER = HERE/'period_worker.py'


def acquire_run_lock(path):
    """Hold a POSIX advisory lock to prevent competing copies of this runner."""
    handle = path.open('a+')
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        raise RuntimeError('Another curve benchmark is running from this folder')
    return handle


def source_hash():
    files = sorted((REPO/'src/genera/curves').rglob('*.py'))
    files += sorted(HERE.glob('*.py'))
    return hashlib.sha256(b''.join(str(p.relative_to(REPO)).encode()+p.read_bytes() for p in files)).hexdigest()


def git_state(path):
    def run(*args):
        p = subprocess.run(['git','-C',str(path),*args],capture_output=True,text=True)
        return p.stdout.strip() if p.returncode==0 else None
    return dict(commit=run('rev-parse','HEAD'),status=run('status','--short'))


def worker(request, sage, seconds):
    command = ([sage,'-python'] if request['engine']=='sage' else [sys.executable])+[str(WORKER)]
    env = dict(os.environ)
    env.pop('PYTHONPATH',None)
    env.pop('PYTHONHOME',None)
    tick = time.perf_counter()
    # Sage MUST start outside the development checkout. Use the same isolated
    # cwd for both workers, and leave Sage's normal HOME/cache configuration alone.
    with tempfile.TemporaryDirectory(prefix='curve-benchmark-') as cwd:
        try:
            process = subprocess.Popen(command,cwd=cwd,env=env,stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                       text=True,start_new_session=True)
        except OSError as exc:
            return dict(status='launch_error',error=str(exc))
        try:
            out,err = process.communicate(json.dumps(request),timeout=seconds)
        except KeyboardInterrupt:
            os.killpg(process.pid,signal.SIGKILL)
            process.communicate()
            for stream in (process.stdin,process.stdout,process.stderr):
                stream.close()
            raise
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL)
            out,err = process.communicate()
            for stream in (process.stdin,process.stdout,process.stderr):
                stream.close()
            return dict(status='timeout',wall_seconds=time.perf_counter()-tick,stderr=err[-4000:])
        try:
            result = json.loads(out)
        except ValueError:
            result = dict(status='protocol_error',error='Worker did not return JSON',stdout=out[-4000:])
        if process.returncode and result.get('status')=='ok':
            result.update(status='worker_exit_error',error=f'Worker exited {process.returncode}')
        result.update(wall_seconds=time.perf_counter()-tick,exit_code=process.returncode,stderr=err[-4000:])
        return result


def report(result):
    lines = ['# Curve-to-full-periods benchmark','',
             'Timing excludes imports, fixed warmup, serialization and validation. Fresh process per sample.',
             'Sage discovers its own basis. Bases are aligned exactly; no fitted factor-of-two correction.',
             'Target digits label the accuracy request. Both engines receive '+
             str(result.get('options',{}).get('working_extra_digits',0))+
             ' extra working digits, including coefficient conversion; their cost is timed.',
             '', '| Case | Digits | Implementation | Successes/trials | Median seconds | Accuracy |',
             '|---|---:|---|---:|---:|---|']
    groups = {}
    for row in result['samples']:
        label = ('hyperelliptic' if row.get('actual_engine')=='hyperelliptic'
                 else row['implementation'])
        groups.setdefault((row['case'],row['digits'],label),[]).append(row)
    for (name,digits,engine),rows in groups.items():
        successes = [r for r in rows if r['status']=='ok']
        times = [r['seconds'] for r in successes]
        checks = [r.get('accuracy',{}).get('status','not_checked') for r in rows]
        accuracy = 'passed' if checks and all(v=='passed' for v in checks) else ', '.join(sorted(set(checks)))
        median = f'{statistics.median(times):.4f}' if times else '—'
        lines.append(f'| {name} | {digits} | {engine} | {len(successes)}/{len(rows)} | {median} | {accuracy} |')
    lines += ['', '## Validated timing ratios','', '| Case | Digits | Sage / genera | Ratio |', '|---|---:|---|---:|']
    for (name,digits,engine),rows in groups.items():
        if engine=='sage':
            continue
        sage = groups.get((name,digits,'sage'),[])
        if not sage or any(r['status']!='ok' or r.get('accuracy',{}).get('status')!='passed'
                           for r in rows+sage):
            continue
        ratio = statistics.median(r['seconds'] for r in sage)/statistics.median(r['seconds'] for r in rows)
        lines.append(f'| {name} | {digits} | {engine} | {ratio:.3f} |')
    lines += ['', 'Ratio > 1 means genera was faster. Single trials are observations, not statistical estimates.',
              'Failures, timeouts and budget skips are retained in results.json and samples.jsonl.',
              'Accuracy tolerance is relative 10^(-digits+5), against a separately computed higher-precision Sage result.',
              'Passing these checks is numerical evidence, not a rigorous error bound.', '',
              f"Source unchanged during run: {result.get('source_unchanged', 'running')}", '']
    for row in result['samples']+result['references']:
        if row['status']!='ok':
            lines.append(f"- {row['case']} ({row['digits']} digits, {row['implementation']}): {row['status']} {row.get('error','')}")
        elif row.get('accuracy',{}).get('status') in ('failed','comparison_error'):
            lines.append(f"- {row['case']} / {row['implementation']}: {row['accuracy']}")
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list',action='store_true')
    parser.add_argument('--preset',choices=('smoke','quick','established','expanded','all'),default='smoke')
    parser.add_argument('--cases',nargs='+',help='Exact names; overrides preset')
    parser.add_argument('--digits',type=int,nargs='+',default=[20])
    parser.add_argument('--trials',type=int,default=1)
    parser.add_argument('--budget',type=float,default=180,help='Whole-run wall budget, seconds')
    parser.add_argument('--timeout',type=float,default=60,help='Worker wall limit including startup')
    parser.add_argument('--sage',default='/usr/local/bin/sage')
    parser.add_argument('--engines',nargs='+',choices=('genera','sage'),default=['genera','sage'])
    parser.add_argument('--reference-extra-digits',type=int,default=10)
    parser.add_argument('--working-extra-digits',type=int,default=15,
                        help='Extra working digits for BOTH engines, including input conversion')
    parser.add_argument('--output',type=Path,help='New results directory; refuses an existing one')
    args = parser.parse_args()
    cases = catalog()
    if args.list:
        for c in cases.values():
            print(f"{c['name']:34s} g={c['genus']} {c['group']:12s} {c['note']}")
        return 0
    if min(args.digits)<15 or args.trials<1 or args.budget<=0 or args.timeout<=0 or args.reference_extra_digits<5 or args.working_extra_digits<0:
        parser.error('digits >=15, trials >=1, positive budgets, reference-extra-digits >=5, working-extra-digits >=0 required')
    if args.cases:
        names = args.cases
    elif args.preset in ('smoke','quick'):
        names = SMOKE if args.preset=='smoke' else QUICK
    else:
        allowed = {'established'} if args.preset=='established' else {'established','expanded'}
        names = [n for n,c in cases.items() if args.preset=='all' or c['group'] in allowed]
    if any(n not in cases for n in names):
        parser.error('Unknown case: '+', '.join(n for n in names if n not in cases))
    try:
        run_lock = acquire_run_lock(HERE/'.run.lock')
    except RuntimeError as exc:
        parser.error(str(exc))
    output = args.output or HERE/'results'/dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(REPO/'src'))
    from mpmath import mp
    import genera  # Ensure genera is importable from development checkout
    started = time.perf_counter()
    original_hash = source_hash()
    result = dict(schema=2,started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  python=sys.version,platform=platform.platform(),repo=str(REPO),
                  git=git_state(REPO),examples_git=git_state(HERE.parent),source_hash=original_hash,
                  options={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                  catalog=[cases[n] for n in names],samples=[],references=[],complete=False)
    def save():
        (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        (output/'report.md').write_text(report(result))
    def execute(case, engine, digits, trial, reference=False):
        mp.dps = digits
        target_bits = mp.prec
        mp.dps = digits + args.working_extra_digits
        row = dict(case=case['name'],digits=digits,bits=target_bits,working_bits=mp.prec,
                   working_digits=digits+args.working_extra_digits,implementation=engine,
                   trial=trial)
        remaining = args.budget-(time.perf_counter()-started)
        if remaining<=0:
            row['status'] = 'skipped_budget'
        else:
            row.update(worker(dict(case=case,engine=engine,bits=mp.prec,repo=str(REPO)),args.sage,min(args.timeout,remaining)))
        if engine=='sage' and row.get('mpmath_path') and Path(row['mpmath_path']).is_relative_to(REPO):
            row.update(status='wrong_import',error='Sage imported development genera/mpmath')
        if case['expected_rejection'] and engine!='sage':
            if row['status']=='error' and row.get('error_type')=='ValueError':
                row['status']='expected_rejection'
            elif row['status']=='ok':
                row['status']='unexpected_acceptance'
        if row['status']=='ok':
            if row.get('genus')!=case['genus'] or row.get('validation') is False:
                row.update(status='validation_failed',error='Genus/internal validation failed')
            if engine!='sage' and case['group'] != 'negative':
                expected_marking = ('baker' if case['name'].startswith('hyper-')
                                    else 'geometric-polygon')
                if row.get('marking')!=expected_marking:
                    row.update(status='routing_mismatch',error=f'Expected marking {expected_marking}')
        (result['references'] if reference else result['samples']).append(row)
        with (output/'samples.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(reference=reference,**row))+'\n')
        print(f"{case['name']} {digits}dps {engine} {'reference' if reference else trial}: "
              f"{row['status']} {row.get('seconds','')} {row.get('error','')}",flush=True)
        save()
        return row
    save()
    for name in names:
        case = cases[name]
        for digits in args.digits:
            cells = []
            mp_engines = [e for e in args.engines if e!='sage']
            engines = mp_engines+(['sage'] if 'sage' in args.engines else [])
            for trial in range(args.trials):
                order = engines if trial%2==0 else list(reversed(engines))
                for engine in order:
                    cells.append(execute(case,engine,digits,trial))
            if 'sage' in args.engines and not case['expected_rejection']:
                ref = execute(case,'sage',digits+args.reference_extra_digits,'reference',True)
                for row in cells:
                    if row['status']!='ok' or ref['status']!='ok':
                        row['accuracy']={'status':'not_checked','reason':'sample or reference unsuccessful'}
                        continue
                    try:
                        check = compare(mp,row,ref,digits)
                        row['accuracy'] = dict(status='passed' if check['passed'] else 'failed',**check)
                    except Exception as exc:
                        row['accuracy'] = dict(status='comparison_error',error=str(exc))
            else:
                for row in cells:
                    row['accuracy']={'status':'not_checked','reason':'Sage disabled or negative case'}
            save()
    result['complete']=True
    result['elapsed_seconds']=time.perf_counter()-started
    result['source_unchanged']=source_hash()==original_hash
    if not result['source_unchanged']:
        for row in result['samples']:
            row['accuracy']={'status':'not_checked','reason':'sources changed during run'}
    save()
    print(f'Results: {output.resolve()}',flush=True)
    failures = [r for r in result['samples']+result['references']
                if r['status'] not in ('ok','expected_rejection') or
                r.get('accuracy',{}).get('status') in ('failed','comparison_error')]
    run_lock.close()
    return 1 if failures or not result['source_unchanged'] else 0


if __name__=='__main__':
    raise SystemExit(main())
