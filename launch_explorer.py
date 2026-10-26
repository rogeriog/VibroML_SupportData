#!/usr/bin/env python3
"""
Root-level launcher for the VibroML GA Phonon Results Explorer.
Forwards all arguments to whitespace_data/setup_and_launch.py.
"""

import runpy
import sys
from pathlib import Path

target = Path(__file__).resolve().parent / "whitespace_data" / "setup_and_launch.py"
sys.path.insert(0, str(target.parent))
runpy.run_path(str(target), run_name="__main__")
