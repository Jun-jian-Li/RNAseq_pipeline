#!/bin/bash

#将结果拷贝进差异分析文件夹

#FPKM
if [ $1 == "single" ];then

    #count
	cp result/05_exp/gene_ref_count.csv 06_diff_rpkm/01_Raw_Data/ref_count/gene_ref_count.csv
	cp result/05_exp/transcript_ref_count.csv 06_diff_rpkm/01_Raw_Data/ref_count/transcript_ref_count.csv

	#rpkm
	cp result/05_exp/gene_ref_rpkm.csv 06_diff_rpkm/01_Raw_Data/ref_rpkm/gene_ref_rpkm.csv
	cp result/05_exp/transcript_ref_rpkm.csv 06_diff_rpkm/01_Raw_Data/ref_rpkm/transcript_ref_rpkm.csv

     #TPM
	cp result/05_exp/gene_ref_TPM.csv 06_diff_rpkm/01_Raw_Data/ref_TPM/gene_ref_TPM.csv
	cp result/05_exp/transcript_ref_TPM.csv 06_diff_rpkm/01_Raw_Data/ref_TPM/transcript_ref_TPM.csv

elif [ $1 == "pair" ];then
    #count
	cp result/05_exp/gene_ref_count.csv 06_diff_fpkm/01_Raw_Data/ref_count/gene_ref_count.csv
	cp result/05_exp/transcript_ref_count.csv 06_diff_fpkm/01_Raw_Data/ref_count/transcript_ref_count.csv


	cp result/05_exp/gene_ref_FPKM.csv 06_diff_fpkm/01_Raw_Data/ref_fpkm/gene_ref_FPKM.csv
	cp result/05_exp/transcript_ref_FPKM.csv 06_diff_fpkm/01_Raw_Data/ref_fpkm/transcript_ref_FPKM.csv

	#TPM
	cp result/05_exp/gene_ref_TPM.csv 06_diff_fpkm/01_Raw_Data/ref_TPM/gene_ref_TPM.csv
	cp result/05_exp/transcript_ref_TPM.csv 06_diff_fpkm/01_Raw_Data/ref_TPM/transcript_ref_TPM.csv
fi



#删除运行完的文件节约计算机储存资源
rm -rf 01_raw_data/*
rm -rf 02_trim_data/*
rm -rf 03_cut_base/*
rm -rf 04_sort_bam/*
rm -rf result/01_raw_data_qc/*
rm -rf result/02_clean_data_qc/*
rm -rf result/03_cut_data_qc/*
rm -rf result/04_bam_qc/*
rm -rf result/05_exp/*

echo "文件移除完成,记得进行下游分析,希望会有好的结果!!!"