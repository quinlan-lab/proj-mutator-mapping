import pandas as pd 

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-bxd"

df = pd.read_csv(PROJDIR + "/data/S1_WS220_CeMEEv2_markerSet1.csv.gz")

print (df)