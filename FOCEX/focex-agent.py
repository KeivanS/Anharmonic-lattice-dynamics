#!/usr/bin/env python3
"""Moved to FOCEX/AGENT/focex-agent.py (2026-10-23). This forwarder keeps the run scripts of
existing projects working; use AGENT/focex-agent.py directly."""
import os
import runpy
import sys

target = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'AGENT', 'focex-agent.py')
sys.argv[0] = target
runpy.run_path(target, run_name='__main__')
