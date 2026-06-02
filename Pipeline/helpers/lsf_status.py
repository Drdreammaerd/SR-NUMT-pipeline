#!/usr/bin/env python3
import sys
import subprocess
import re

job_id = sys.argv[1]

# Extract job ID if it contains extra text (Snakemake sometimes passes 'Job <12345> is submitted...')
match = re.search(r'(\d+)', job_id)
if match:
    job_id = match.group(1)

try:
    res = subprocess.run(["bjobs", "-o", "stat", "-noheader", job_id], 
                         capture_output=True, text=True)
    
    if res.returncode != 0:
        # If job is not found, LSF returns non-zero. Assume it was killed and purged.
        print("failed")
        sys.exit(0)
        
    status = res.stdout.strip()
    
    if status == "DONE":
        print("success")
    elif status in ["EXIT", "ZOMBI", "UNKWN"]:
        print("failed")
    elif status in ["PEND", "RUN", "SUSP", "USUSP", "SSUSP", "PSUSP"]:
        print("running")
    else:
        print("running")
except Exception:
    print("failed")
