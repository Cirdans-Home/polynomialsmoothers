#!/usr/bin/bash -l

allTasks=(1 2 4 8 16 32 64)
allGPUS=(1 2 4 8 16 32 64)
allNodes=(1 1 1 2 4 8 16)
idim=(2000 2840 4017 5681 8034 11362 16069)
script=$1
outscript=${script%.sh}
i=0

for batchnp in ${allTasks[@]};
do
 echo "Building script with ${batchnp} tasks and ${allGPUS[${i}]} gpus size is ${idim[$i]}"
 if [[ ${allNodes[$i]} -ge 64  ]]
 then
	queue='#SBATCH --qos=boost_qos_bprod'	
 else
	queue=''
 fi
 if [[ ${allGPUS[$i]} -le 4 ]]
 then
	gpunumber=${allGPUS[${i}]}
 else
	gpunumber=4
 fi
 cat $1 | sed "s/@NNODES@/${allNodes[$i]}/g
 s/@NGPUS@/${gpunumber}/g
 s/@SIZE@/${idim[$i]}/g
 s/@NTASKS@/${allTasks[$i]}/g 
 s/@MORETHANONENODE@/${queue}/g
 s/@DEGITER@/${2}/g" > ${outscript}-$batchnp.sh
 let i=$i+1
done
