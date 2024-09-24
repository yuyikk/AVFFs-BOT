#!/bin/zsh

# Total events=MULT*pow(10, POW)

export POW="10"
export MULT="1"
export OUTDIR="./output/0915-0916.Virtual"
export INPUT="./input/0915-0916.Virtual.root"


# export POW="7"
# export MULT="300"

# export OUTDIR="./output/0915.Virtual.Vacuum"
# export INPUT="./input/0915.Virtual.Vacuum.root"

# export POW="7"
# export MULT="394"

# export OUTDIR="./output/0918.Vacuum.Bfield1T"
# export INPUT="./input/0918.vir.vacuum.root"


mkdir -p $OUTDIR
root -l -b -q drawVirtual.cc