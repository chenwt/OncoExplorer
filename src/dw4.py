#dw3
from multiprocessing import Pool
import subprocess
import glob,os
import random
#deepwalk --output result2/brca_km40_ss_wl1_nw10.embeddings --number-walks 10 --walk-length 1 --workers 64 --representation-size 40 --input data/graph_brca.mat --format mat --matfile-variable-name network

# function for multi-thread processing.
def f(wl, nw):
    subprocess.call("deepwalk --output result2/brca_km40_ss_wl"+str(wl)+"_nw"+str(nw)+".embeddings --number-walks "+str(nw)+" --walk-length "+str(wl)+" --workers 64 --representation-size 40 --input data/graph_brca.mat --format mat --matfile-variable-name network", shell=True)

for wl in [40]:
	for nw in [10,20,50,100,200]:
		f(wl, nw)
# EOF.