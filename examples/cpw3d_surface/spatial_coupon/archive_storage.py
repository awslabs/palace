# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Conservative cross-worker reservations for shared response-archive storage."""
import fcntl
import json
import os
import shutil
from pathlib import Path


def save(path,value):
    temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def reserve(root,key,bytes_required,margin=100*2**30,free_bytes=None):
    root=Path(root)
    if bytes_required<=0:raise ValueError('Invalid storage reservation')
    with (root/'archive-storage.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=root/'archive-storage.json'
        ledger=json.loads(path.read_text()) if path.exists() else {}
        if key in ledger:raise ValueError('Existing archive reservation requires explicit recovery')
        remaining=sum(r['Bytes'] for r in ledger.values() if not r.get('Sealed',False))
        # Until sealed, count the entire reservation, even if it has already
        # partly consumed disk. This is conservative, never optimistic.
        available=shutil.disk_usage(root).free if free_bytes is None else free_bytes
        if available-remaining-bytes_required<margin:
            raise RuntimeError(f'Insufficient unreserved storage: free={available}, outstanding={remaining}, requested={bytes_required}, margin={margin}')
        ledger[key]={'Bytes':int(bytes_required),'Sealed':False}
        save(path,ledger)
        return ledger[key]


def seal(root,key,actual_bytes):
    root=Path(root)
    with (root/'archive-storage.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=root/'archive-storage.json';ledger=json.loads(path.read_text())
        if actual_bytes>ledger[key]['Bytes']:
            raise RuntimeError('Archive exceeded its reserved bound')
        ledger[key].update(Sealed=True,ActualBytes=int(actual_bytes))
        # Existing files remain allocated and are reflected in filesystem free
        # space. Sealed reservations no longer promise future archive writes.
        save(path,ledger)


def raw_tetrahedral_h1_dofs(counts,order):
    """Raw topology diagnostic only; cracking/refinement can change this count."""
    if order<1:raise ValueError('Invalid H1 order')
    p=order
    return (counts['Vertices']+(p-1)*counts['Edges']+
            (p-1)*(p-2)//2*counts['Faces']+
            (p-1)*(p-2)*(p-3)//6*counts['Elements'])
