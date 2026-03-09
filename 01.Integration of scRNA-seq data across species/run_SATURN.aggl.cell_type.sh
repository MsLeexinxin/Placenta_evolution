#!/bin/bash

size=800

mkdir -p ${input_data}_max_split_size_mb_${size}_mac${nmac}_hv0_aggl
cd ${input_data}_max_split_size_mb_${size}_mac${nmac}_hv0_aggl
mkdir -p multiple_seeds_results
ln -fs ../${input_data} .

echo "path,species,embedding_path" >${input_data}.csv
for i in $PWD/${input_data}/*.h5ad
do
    dir=`dirname $i`
    sp1=`basename $i|cut -f 1 -d '_'|cut -c 1-4`
    sp2=`basename $i|cut -f 1 -d '.'`
    if [ -e $dir/$sp1.*.pt ]
    then
        emb=`ls -1 $dir/$sp1.*.pt`
        echo -e "$i,$sp2,$emb"
    fi
done >>${input_data}.csv

source /share/home/zhanglab/user/lixin/miniconda3/etc/profile.d/conda.sh
conda activate python39


############################################################################################################################
echo $SLURMD_NODENAME

echo "----------------------------------------------------------------------------------------"
path=`which nvidia-smi`
echo
ls -lh $path
echo
nvidia-smi

while [ -z $device ]
do
    device=`nvidia-smi|grep 'Default'|awk '{gsub(/%/,"",$13);if($13<10 && !a){print FNR-1;a++}}'`
    sleep 10s
done

export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:$size

python3 -u /share/home/zhanglab/user/lixin/software/SATURN/train-saturn.ncores.aggl.py \
    --in_data=./${input_data}.csv \
    --device_num=$device \
    --pretrain_batch_size 10240 \
    --batch_size 10240 \
    --in_label_col=cell_type \
    --ref_label_col=cell_type \
    --work_dir=./multiple_seeds_results/ \
    --num_macrogenes=${nmac} \
    --pretrain \
    --model_dim=256 \
    --polling_freq=201 \
    --epochs=50 \
    --hv_genes=0 \
    --hv_span=1 \
    --pretrain_epochs=200 \
    --pe_sim_penalty=1.0 \
    --l1_penalty=0 \
    --centroid_score_func=default \
    --seed=0 \
    --org=${input_data}_l1_0_pe_1.0_ESM2_macrogenes_${nmac}_hv_genes_0_centroid_score_func_default_batch_label_split \
    --embedding_model=ESM2 >s_mac${nmac}.log 2>&1

input=`ls -1 ./multiple_seeds_results/saturn_results/*seed_0.h5ad`
prefix=`echo $input|sed "s/_org_${input_data}/ /g"|awk '{print $1}'|sed 's/test256_data_//g'`
ln -fs $input $prefix.h5ad
