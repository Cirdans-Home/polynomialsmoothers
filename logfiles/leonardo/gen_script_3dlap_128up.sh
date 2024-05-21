#!/usr/bin/bash -l

allTasks=(128 256 512 1024)
allGPUS=(128 256 512 1024)
allNodes=(32 64 128 256)
idim=(917 1155 1455 1833)
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
