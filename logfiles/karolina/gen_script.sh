#!/usr/bin/bash -l

allTasks=(1 2 4 8 12 16 32 64 128)
allGPUS=(1 2 4 8 12 16 32 64 128)
allNodes=(1 1 1 1 2 2 4 8 16)
idim=(2449 3464 4899 6928 9798 13856 19596 27713)
script=$1
outscript=${script%.sh}
i=0

for batchnp in ${allTasks[@]};
do
 cat $1 | sed "s/@NNODES@/${allNodes[$i]}/g
 s/@NGPUS@/${allGPUS[$i]}/g
 s/@SIZE@/${idim[$i]}/g
 s/@NTASKS@/${allTasks[$i]}/g" > ${outscript}-$batchnp.sh
 let i=$i+1
done
