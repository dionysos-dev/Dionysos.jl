#!/usr/bin/env python3
"""Serve a read-only, local progress view of the single-query evaluation."""
import argparse
from datetime import datetime, timezone
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import subprocess
import time
from urllib.parse import urlsplit


def read_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return default


def modified_at(path):
    try:
        return path.stat().st_mtime
    except OSError:
        return None


def initial_inventory(output):
    """Recover the original missing-work total from the first saved inventory."""
    try:
        with (output/'run.log').open() as log:
            for line in log:
                if line.startswith('INVENTORY '):
                    return json.loads(line[len('INVENTORY '):])
    except (OSError, json.JSONDecodeError):
        return None
    return None


def runner_alive(background):
    pid = background.get('pid')
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        result = subprocess.run(['ps', '-p', str(pid), '-o', 'stat=,command='],
                                capture_output=True, text=True, timeout=3)
    except (OSError, subprocess.TimeoutExpired):
        return False
    fields = result.stdout.strip().split(maxsplit=1)
    return (result.returncode == 0 and len(fields) == 2 and 'Z' not in fields[0]
            and any(name in fields[1] for name in
                    ('CalibrationTable/single_query_baseline.py', 'CalibrationTable/resume_baseline.py')))


class ProgressReader:
    def __init__(self, output):
        self.output = output
        self.normalization = {}
        self.normalization_mtime = None
        self.original_inventory = None

    def current_case(self, formula):
        if not formula:
            return None
        path = self.output/'normalization.json'
        stamp = modified_at(path)
        if stamp != self.normalization_mtime:
            self.normalization = read_json(path, {})
            self.normalization_mtime = stamp
        normalized = self.normalization.get(formula)
        if normalized is None:
            return None
        job_id = hashlib.sha256(normalized.encode()).hexdigest()
        return self.output/'cases'/job_id

    def snapshot(self):
        now = time.time()
        summary = read_json(self.output/'summary.json', {})
        progress = read_json(self.output/'progress.json', {})
        inventory = read_json(self.output/'inventory.json', {})
        stats = summary.get('statistics', {})
        total = stats.get('n_prompts', inventory.get('test_prompts', 0))
        pending = stats.get('pending_prompts')
        batch_total = inventory.get('pending_prompts')
        resolved = total-pending if pending is not None else None
        batch_done = max(0, batch_total-pending) if batch_total is not None and pending is not None else None
        if self.original_inventory is None:
            self.original_inventory = initial_inventory(self.output)
        original_total = (self.original_inventory or {}).get('pending_prompts')
        original_done = original_total-pending if original_total is not None and pending is not None else None
        completed_before_batch = original_total-batch_total if original_total is not None and batch_total is not None else None
        stage = progress.get('stage', 'preparing')
        complete = bool(summary.get('complete'))
        alive = runner_alive(read_json(self.output/'background-process.json', {}))
        if alive and stage == 'stopped' and not complete:
            progress = {}
            stage = 'preparing'
        if complete:
            state = 'complete'
        elif not alive:
            state = 'stopped'
        elif stage == 'cooldown':
            state = 'resting'
        else:
            state = 'running'
        case = self.current_case(progress.get('formula'))
        progress_mtime = modified_at(self.output/'progress.json')
        worker_path = case/'worker-stage.json' if case else None
        worker_mtime = modified_at(worker_path) if worker_path else None
        substage = None
        if progress.get('synthesis_process_active') and worker_mtime is not None and worker_mtime >= (progress_mtime or 0):
            substage = read_json(worker_path, {}).get('stage')
        return dict(generated_at=now, state=state, complete=complete, runner_alive=alive,
                    total=total, resolved=resolved, pending=pending,
                    batch_total=batch_total, batch_done=batch_done,
                    original_missing_total=original_total, original_missing_done=original_done,
                    completed_before_batch=completed_before_batch,
                    current_index=progress.get('pending_index'), stage=stage, substage=substage,
                    stage_started_at=worker_mtime if substage else progress_mtime,
                    last_saved_at=modified_at(self.output/'summary.json'),
                    prompt_index=progress.get('prompt_index'), record_id=progress.get('record_id'),
                    formula=progress.get('formula'), reason=progress.get('reason'),
                    seconds_remaining=progress.get('seconds_remaining'),
                    error=progress.get('interruption_or_error'),
                    final_statistics=stats if complete else None)


def handler_for(reader, template):
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            path = urlsplit(self.path).path
            if path == '/':
                content, mime = template.read_bytes(), 'text/html; charset=utf-8'
            elif path == '/api/progress':
                content = json.dumps(reader.snapshot(), ensure_ascii=False).encode()
                mime = 'application/json; charset=utf-8'
            elif path == '/favicon.ico':
                self.send_response(204)
                self.end_headers()
                return
            else:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header('Content-Type', mime)
            self.send_header('Content-Length', str(len(content)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; connect-src 'self'; base-uri 'none'")
            self.end_headers()
            try:
                self.wfile.write(content)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def log_message(self, format, *args):
            pass

    return Handler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--port', type=int, default=8769)
    parser.add_argument('--snapshot', action='store_true')
    args = parser.parse_args()
    if not 0 <= args.port <= 65535:
        parser.error('Port must be between 0 and 65535')
    root = args.root.resolve()
    output = (args.output or root/'Results/deepseek-single-query-baseline-cal300-v1').resolve()
    if not output.is_dir():
        parser.error('The evaluation output directory does not exist')
    reader = ProgressReader(output)
    if args.snapshot:
        print(json.dumps(reader.snapshot(), indent=2))
        return
    template = Path(__file__).with_name('progress_dashboard.html')
    server = ThreadingHTTPServer(('127.0.0.1', args.port), handler_for(reader, template))
    server.daemon_threads = True
    port = server.server_address[1]
    url = f'http://127.0.0.1:{port}/'
    metadata = dict(pid=os.getpid(), url=url,
                    started_at=datetime.now(timezone.utc).isoformat(), read_only=True)
    destination = output/'progress-dashboard.json'
    temporary = destination.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(metadata, indent=2)+'\n')
    temporary.replace(destination)
    print(url, flush=True)
    try:
        server.serve_forever(poll_interval=1)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == '__main__':
    main()
