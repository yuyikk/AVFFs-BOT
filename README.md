# AVFFs-BOT
The Beam-on-target simulations on neutron arm.
## How to build
1. change to the directory where `src` and `include` located.
2. `mkdir -p build`
3. `cd build`
4. `cmake ..`
5. `make`

    You will get a executable file `doSimu` in the current directory.

## How to run
1. `cd build`
2. `./doSimu number_of_events runnumber`. 

    The `doSimu` accepts two arguements, the first is the number of events and the second is the runnumber that you can set freely.
    
    You will get a root file `runnumber.BOT.root` in the current directory, in which a tree is stored and also a histogram storing number of generated enents

## Read tree
Some root macros of reading simulation trees can be found in `readTree`.
### How to run
1. `compile processVirtual.cc`
2. `./processVritual input_root_file output_root_file`

To draw plots, please modify the `draw.zsh`. There are some variables defined in this script.
```bash
# POW and MULT are used to calculate the number of generated events 
# nevent = MULT * pow(10, POW)
export POW="10" 
export MULT="1"
## directory where plots are stored
export OUTDIR="./output/0915-0916.Virtual" 
## input root files obtained with `processVirtual`
export INPUT="./input/0915-0916.Virtual.root" 
```
## Some words about the codes
The size of detectors, for example, are defined in `include/Constant.hh`

Two virtual detectors are placed at Magnet and TOF with copy numbers 10001 and 10002.