#!/bin/bash
##############################################################################
# author      : Ljj
# description : RNAseq pipeline step2 count-> FPKM matrix-> TPM matrix
##############################################################################


#reads average length 作为 -l 的输入,结果如下:
#zcat WT-1h-1_1.fq.gz | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'
#139.759 平均的read长度都是在50左右

echo "请质控以及过滤后的reads的长度 ,输入文件样子 bash 03_extract_exp.sh 150 pair/single" ;
input_dir=./05_quantify
out_dir=./result/05_exp
echo "定量完成,记得生成质量控制报告!!!" ;

if [ $2 == "single" ];then
	#rpkm
	#提取gene对应的count文件
	python2 ./00_script/prepDE.py -l $1 -i $input_dir -g $out_dir/gene_ref_count.csv -t $out_dir/transcript_ref_count.csv
	
	#提取gene对应的rpkm文件
	python2 ./00_script/getFPKM.py -l $1 -i $input_dir -g $out_dir/gene_ref_rpkm.csv -t $out_dir/transcript_ref_rpkm.csv

	#提取gene对应的TPM文件
	python2 ./00_script/getTPM.py -l $1 -i $input_dir -g $out_dir/gene_ref_TPM.csv -t $out_dir/transcript_ref_TPM.csv

elif [ $2 == "pair" ];then
	#提取gene对应的count文件
	python2 ./00_script/prepDE.py -l $1 -i $input_dir -g $out_dir/gene_ref_count.csv -t $out_dir/transcript_ref_count.csv 

	#提取gene对应的FPKM文件
	python2 ./00_script/getFPKM.py -l $1 -i $input_dir -g $out_dir/gene_ref_FPKM.csv -t $out_dir/transcript_ref_FPKM.csv

	#提取gene对应的TPM文件
	python2 ./00_script/getTPM.py -l $1 -i $input_dir -g $out_dir/gene_ref_TPM.csv -t $out_dir/transcript_ref_TPM.csv
fi
