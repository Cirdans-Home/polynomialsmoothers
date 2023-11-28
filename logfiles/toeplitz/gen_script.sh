#!/bin/bash -l

allTasks=(1 2 4 8 12 16 32 64 128 216)
idim=(500 1000 2000 4000 8000 16000 32000 64000 10800)
script=$1
outscript=${script%.sh}
i=0
nthread=1
theta=30
eps=100

for batchnp in ${allTasks[@]};
do
 cat $1 | sed "s/@NTASK@/$batchnp/g
 s/@NTHREAD@/$nthread/g
 s/@SIZE@/${idim[$i]}/g
 s/@THETA@/$theta/g
 s/@EPS@/$eps/g" > ${outscript}-$batchnp.sh
 let i=$i+1
done
