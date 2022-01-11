import pandas as pd
import os

SRR_IDs = pd.read_csv("resources/accession_info/Isakova_smRNA_atlas_SRR_IDs.csv")["SRRID"]
num_threads = 30
outdir="/workdir/dwm269/totalRNA/data/fastq/Isakova_smRNAseq"
tmpdir="/workdir/dwm269/totalRNA/tmp"

print(SRR_IDs)

for SRR in SRR_IDs:

    print("Downloading .fastq's for " + str(SRR))
    os.system("prefetch --max-size 999999999999 " + str(SRR))
    os.system("parallel-fastq-dump --sra-id " + str(SRR) + " --threads " + str(num_threads) + " --outdir " + outdir + " --tmpdir " + tmpdir + " --split-files --gzip")
    # shell("mv "+outdir+"")
