
# Source: https://github.com/BUStools/BUS_notebooks_python/blob/master/dataset-notebooks/10x_neuron_1k_v2chem_python/10x_neuron_1k_v2chem.ipynb
import os
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, collections

# read in arguments
outdir = str(sys.argv[1])
t2g_path = str(sys.argv[2])
bustext = str(sys.argv[3])


# Conversion script, almost exactly the same as in the link above
tr2g = {}
trlist = []
with open(t2g_path) as f:
    for line in f:
        l = line.split()
        tr2g[l[0]] = l[1]
        trlist.append(l[0])

genes = list(set(tr2g[t] for t in tr2g))

# load equivalence classes
ecs = {}
with open(outdir+'/matrix.ec') as f:
    for line in f:
        l = line.split()
        ec = int(l[0])
        trs = [int(x) for x in l[1].split(',')]
        ecs[ec] = trs

def ec2g(ec):
    if ec in ecs:
        return list(set(tr2g[trlist[t]] for t in ecs[ec]))
    else:
        return []

cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
pbar=None
pumi=None
with open(bustext) as f:
    gs = set()
    for line in f:
        l = line.split()
        barcode,umi,ec,count = line.split()
        ec = int(ec)

        if barcode == pbar:
            # same barcode
            if umi == pumi:
                # same UMI, let's update with intersection of genelist
                gl = ec2g(ec)
                gs.intersection_update(gl)
            else:
                # new UMI, process the previous gene set
                for g in gs:
                    cell_gene[barcode][g] += 1.0/len(gs)
                # record new umi, reset gene set
                pumi = umi
                gs = set(ec2g(ec))
        else:
            # work with previous gene list
            for g in gs:
                cell_gene[pbar][g] += 1.0/len(gs)

            if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                del cell_gene[pbar]

            pbar = barcode
            pumi = umi

            gs = set(ec2g(ec))

    for g in gs:
        cell_gene[pbar][g] += 1.0/len(gs)

    if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
        del cell_gene[pbar]

barcode_hist = collections.defaultdict(int)
for barcode in cell_gene:
    cg = cell_gene[barcode]
    s = len([cg[g] for g in cg])
    barcode_hist[barcode] += s

#Output a gene count histogram
# bcv = [x for b,x in barcode_hist.items()]
# plt.switch_backend('agg')
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.hist(bcv,bins=100)
# ax.set_title("Histogram")
# plt.xlabel("number of genes detected")
# plt.ylabel("number of barcodes")
# fig.savefig(outdir+'gene_hist.png')

outfile = outdir+'/counts_unfiltered/matrix.mtx'

gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
barcodes_to_use = [b for b,x in barcode_hist.items()]

num_entries = 0
for barcode in barcodes_to_use:
    num_entries += len([x for x in cell_gene[barcode].values() if x>0])

with open(outfile, 'w') as of:
    of.write('%%MatrixMarket matrix coordinate real general\n%\n')
    #number of genes
    of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), round(num_entries)))
    bcid = 0
    for barcode in barcodes_to_use:
        bcid += 1
        cg = cell_gene[barcode]
        gl = [(gene_to_id[g],cg[g]) for g in cg if cg[g] > 0]
        gl.sort()
        for x in gl:
            of.write("%d %d %f\n"%(x[0],bcid,x[1]))

gene_names = {}
with open(t2g_path) as f:
    f.readline()
    for line in f:
        t,ens,g,gn,chrom,start,stop,strand = line.split()
        gene_names[g] = gn

id_to_genes = dict((i,g) for (g,i) in gene_to_id.items())
gl = []
for i in range(1,len(genes)+1):
    g = id_to_genes[i]
    gid = g
#    gid = g[:g.find('.')]
    if gid in gene_names:
        gn = gene_names[gid]
    else:
        gn = ''
    gl.append((g,gn))

with open(outdir+'/counts_unfiltered/genes.tsv','w') as of:
    for g,gn in gl:
        of.write("%s\t%s\n"%(g,gn))

with open(outdir+'/counts_unfiltered/barcodes.tsv','w') as of:
    of.write('\n'.join(x + '' for x in barcodes_to_use))
    of.write('\n')
