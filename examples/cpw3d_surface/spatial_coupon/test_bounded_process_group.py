# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import importlib.util
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

spec=importlib.util.spec_from_file_location('bounded',Path(__file__).with_name('run_bounded_mesher.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@unittest.skipUnless(os.name=='posix','POSIX process-group guard')
class ProcessGroupTest(unittest.TestCase):
    def test_term_ignoring_child_does_not_survive_launcher(self):
        with tempfile.TemporaryDirectory() as tmp:
            pidfile=Path(tmp)/'child.pid'
            child=f'import os,signal,time,pathlib; signal.signal(signal.SIGTERM,signal.SIG_IGN); pathlib.Path({str(pidfile)!r}).write_text(str(os.getpid())); time.sleep(60)'
            parent=f'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c",{child!r}]); time.sleep(60)'
            process=subprocess.Popen([sys.executable,'-c',parent],start_new_session=True)
            try:
                deadline=time.monotonic()+3
                while not pidfile.exists() and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue(pidfile.exists())
                pid=int(pidfile.read_text())
                module.stop(process)
                deadline=time.monotonic()+3
                state=''
                while time.monotonic()<deadline:
                    state=subprocess.run(['ps','-o','stat=','-p',str(pid)],capture_output=True,text=True).stdout.strip()
                    if not state or state.startswith('Z'):break
                    time.sleep(.02)
                self.assertTrue(not state or state.startswith('Z'),state)
            finally:
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                process.wait(timeout=5)


if __name__=='__main__':unittest.main()
