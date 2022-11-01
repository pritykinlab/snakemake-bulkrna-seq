#!/bin/bash

sbcmd="sbatch --mem={cluster.mem} --cpus-per-task={cluster.cpus}"
sbcmd+=" --time={cluster.time} --nodes={cluster.n} --qos={cluster.qos}"

snakemake -j 18 --rerun-incomplete --keep-going \
    --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all
