# Trigger-Efficiency-vs-SNR
Calculating the trigger efficiency vs SNR for the 2/4 hit multiplicity trigger

1. Run make_dag.py to generate att.dag
2. Run att.dag to obtain an array of snr values for each attenuation (utilized make_snr_bin.py)
3. Run plot_snr_att.py to obtain a plot of snr vs the attenuation 
4. Run plot_trig_snr.py to obtain a plot of the trigger efficiency vs SNR. 
5. Can modify make_dag.py for select attenuations/channels
