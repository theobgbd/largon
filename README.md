# largon
A simple Monte-Carlo/Metropolis simulation code for a Lennard Jones fluid.

## Compilation :
`./compile.sh`

The script will produce both a serial and a parallel executable. 

## Execution :
`./largon.x < input` will execute the serial code

`./largon-p.x < input` will execute the parallel version

## Input file sample
```
outputfile statsfile
trajfile traj.xyz
nstep 50
natoms 10
temperature 300
```
The program stores the internal energy in `outputfile` and the sampled configurations in `trajfile`, as a `.xyz` format.
