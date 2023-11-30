#!/bin/bash -l

allTasks=(1 2 4 8 12 16 32 64 128 216)
idim=(500 707 1000 1414 2000 2828 4000 5657 7348)
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
