#!/usr/bin/bash -l

allTasks=(1 2 4 8 16 32 64 128 256 512 1024)
allGPUS=(1 2 4 8 16 32 64 128 256 512 1024)
allNodes=(1 1 1 2 4 8 16 32 64 128 512)
idim=(2449 3464 4898 6928 9797 13856 19595 27712 39191 55425 78383)
script=$1
outscript=${script%.sh}
i=0

for batchnp in ${allTasks[@]};
do
 echo "Building script with ${batchnp} tasks and ${allGPUS[${i}]} gpus size is ${idim[$i]}"
 if [[ ${allNodes[$i]} -ge 2  ]]
 then
	queue='#SBATCH --qos=a100multi'	
 else
	queue=''
 fi
 if [[ ${allGPUS[$i]} -le 8 ]]
 then
	gpunumber=${allGPUS[${i}]}
 else
	gpunumber=8
 fi
 cat $1 | sed "s/@NNODES@/${allNodes[$i]}/g
 s/@NGPUS@/${gpunumber}/g
 s/@SIZE@/${idim[$i]}/g
 s/@NTASKS@/${allTasks[$i]}/g 
 s/@MORETHANONENODE@/${queue}/g
 s/@DEGITER@/${2}/g" > ${outscript}-$batchnp.sh
 let i=$i+1
done
