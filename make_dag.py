#!/usr/env python 
import numpy as np

#run make_snr_bin.py for specified attenuations (att_unique) and channels (chs)

att_unique  = np.arange(0, 32, 0.5)
chs = [0,1,2,3]

dag = ""
for att in att_unique:
    for ch in chs:
        job_name = f"{att}_{ch}_sim"
        cmd = f"JOB {job_name} att.sub\n"
        cmd += f"VARS {job_name} job_name=\"{job_name}\""
        cmd += f" cmd=\"'python3 /data/condor_builds/users/avijai/RNO/tutorials-rnog/get_daqstatus/make_snr_bin.py "
        cmd += f"--file /data/i3store/users/avijai/rnog_tutorials/station23/run3400/combined.root "
        cmd += f"--att {att} "
        cmd += f"--ch {ch} "
        cmd += "'"
        cmd += f"\"\n"
        dag += cmd

open("att.dag", 'w').write(dag)





