#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Separate, six-job source-only build against existing read-only Linux dependencies
(the recorded build_linux_diagnostic.py of coupon-accuracy-assessment-20260913/experiments,
re-pointed at the decision-62(4) source freeze source-freeze-perf62 -> linux-build-perf62)."""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

ROOT = Path('/data/home/simlap/coupon_accuracy_assessment_20260913')
OLD = Path('/data/home/simlap/palace_builds/palace')
DEPS = OLD/'build_c8g_main'
FROZEN = Path('/data/home/simlap/transmon_library_speed_20260909/bin/palace-fast')
FROZEN_SHA = 'e1d14c903c6b0b41ce65221a57c569ab18c9c804e9b0ac0a6e1dfcffd8cf6cbd'
SOURCE = ROOT/'source-freeze-perf62'
BUILD = ROOT/'linux-build-perf62'

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1048576), b''):
            h.update(data)
    return h.hexdigest()

def save(name, data):
    with (BUILD/name).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')

assert sha(FROZEN) == FROZEN_SHA
assert not BUILD.exists()
expected = json.loads((SOURCE/'manifest.json').read_text())
assert all(sha(SOURCE/p) == h for p,h in expected['Files'].items())
BUILD.mkdir()
commands = json.loads((DEPS/'palace-build/compile_commands.json').read_text())
commands = [c for c in commands if c['output'].startswith(('CMakeFiles/libpalace.dir/', 'CMakeFiles/palace.dir/'))]
assert len(commands) == 102
rewritten = []
for record in commands:
    cmd = shlex.split(record['command'])
    source = SOURCE/'palace'/Path(record['file']).relative_to(OLD/'palace')
    assert str(source.relative_to(SOURCE)) in expected['Files']
    cmd[cmd.index('-c')+1] = str(source)
    cmd = [arg.replace('-I'+str(OLD/'palace'), '-I'+str(SOURCE/'palace'))
           .replace('-I'+str(DEPS/'palace-build/generated'), '-I'+str(SOURCE/'generated')) for arg in cmd]
    obj = BUILD/record['output']
    obj.parent.mkdir(parents=True, exist_ok=True)
    cmd[cmd.index('-o')+1] = str(obj)
    rewritten.append({'Command':cmd, 'Object':str(obj), 'Source':str(source)})
link = shlex.split((DEPS/'palace-build/CMakeFiles/palace.dir/link.txt').read_text())
for i,arg in enumerate(link):
    if arg == 'libpalace.a' or arg.startswith('CMakeFiles/palace.dir/'):
        link[i] = str(BUILD/arg)
link[link.index('-o')+1] = str(BUILD/'palace-diagnostic.bin')
# Inventory contents, not just paths, of dependencies and configured build inputs.
inputs = {DEPS/'palace-build/compile_commands.json', DEPS/'palace-build/CMakeFiles/palace.dir/link.txt', FROZEN,
          DEPS/'palace-build/palace-arm64.bin', Path(__file__).resolve(), Path('/opt/gcc/bin/g++'),
          Path('/opt/openmpi/bin/mpirun')}
inputs.update(p for p in (DEPS/'include').rglob('*') if p.is_file())
inputs.update(p for p in (DEPS/'lib').glob('*') if p.is_file())
inputs.update(Path(p) for p in link if p.startswith('/') and Path(p).is_file())
ldd_before = subprocess.check_output(['ldd', str(FROZEN)], text=True)
for token in ldd_before.split():
    if token.startswith('/') and Path(token).is_file():
        inputs.add(Path(token))
inputs.update(Path(p) for p in ['/opt/gcc/libexec/gcc/aarch64-linux-gnu/14.3.0/cc1plus'] if Path(p).is_file())
inputs = sorted(inputs)
before = {str(p):sha(p) for p in inputs}
save('provenance-before.json', {'InputSHA256':before, 'SourceManifest':expected,
     'CompilerVersion':subprocess.check_output(['/opt/gcc/bin/g++','--version'], text=True),
     'MPI':subprocess.check_output(['/opt/openmpi/bin/mpirun','--version'], text=True),
     'Environment':{k:os.environ.get(k) for k in ['PATH','LD_LIBRARY_PATH','LOADEDMODULES']},
     'FrozenLdd':ldd_before, 'Jobs':6, 'CompileCommands':rewritten, 'LinkCommand':link})
status = 'incomplete'
try:
    def compile_one(record):
        print(shlex.join(record['Command']), flush=True)
        subprocess.run(record['Command'], cwd=BUILD, check=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(compile_one, rewritten))
    objects = [c['Object'] for c in rewritten if '/libpalace.dir/' in c['Object']]
    subprocess.run(['/usr/bin/ar','qc',str(BUILD/'libpalace.a'),*objects], check=True)
    subprocess.run(['/usr/bin/ranlib',str(BUILD/'libpalace.a')], check=True)
    print(shlex.join(link), flush=True)
    subprocess.run(link, cwd=BUILD, check=True)
    exe = BUILD/'palace-diagnostic.bin'
    digest = sha(exe)
    binary = ROOT/('palace-archive-estimate-'+digest+'.bin')
    with binary.open('xb') as stream:
        stream.write(exe.read_bytes())
    binary.chmod(0o555)
    ldd = subprocess.check_output(['ldd',str(binary)], text=True)
    dynamic = {token:sha(token) for token in ldd.split() if token.startswith('/') and Path(token).is_file()}
    save('binary.json', {'Path':str(binary), 'SHA256':digest, 'Ldd':ldd,
                        'DynamicDependenciesSHA256':dynamic})
    status = 'complete'
finally:
    changed = [str(p) for p in inputs if not p.exists() or sha(p) != before[str(p)]]
    changed_source = [p for p,h in expected['Files'].items() if sha(SOURCE/p) != h]
    save('provenance-after.json', {'Status':status,'ChangedInputs':changed,'ChangedSources':changed_source})
    assert not changed and not changed_source
