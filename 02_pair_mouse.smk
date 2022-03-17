configfile: "env/config.yaml"
import os
import re

samples = []
reads = []
for item in os.listdir(os.path.abspath(".")+"/"+config["dir"]["raw_data_dir"]+"/"):
     #输出指定后缀类型的文件
     if(item.endswith('gz')):
         m=re.findall(r'(.*)_(.*\d).*q.gz',item)       #findall匹配返回一个list
         sample_name_list = [x[0] for x in m]          #提取样本命的
         sample_name_list="".join(sample_name_list)    #list 转 str
         samples.append(sample_name_list)              #样本加入list
         samples=list(set(samples))                    #去除重复样本名字
         samples = [i for i in samples if(len(str(i))!=0)]#去除空
     
         #提取Read信息列表
         read_name_list = [x[1] for x in m]            #提取read名
         read_name_list="".join(read_name_list)        #list 转 str
         reads.append(read_name_list)                  #样本加入list
         reads=list(set(reads))                        #去除重复样本名字
         reads = [i for i in reads if(len(str(i))!=0)] #去除空


print (samples) #获取的样本的名字
print (reads)   #获取的样本的测序read方向

rule all:
     input:
        raw_data_fastqc = expand(os.path.join( config["dir"]["raw_data_qc_result_dir"],"{sample}_{R}_fastqc.html" ), sample=samples, R=reads),
        clean = expand(os.path.join( config["dir"]["clean_data_dir"], "{sample}_{trim_read}.fq.gz"),sample=samples,trim_read=['R1_val_1','R2_val_2']),
        clean_data_fastqc1=expand(os.path.join( config["dir"]["clean_data_qc_result_dir"],"{sample}_R1_val_1_fastqc.html"),sample=samples),
        clean_data_fastqc2=expand(os.path.join( config["dir"]["clean_data_qc_result_dir"],"{sample}_R2_val_2_fastqc.html"),sample=samples),
        cut_base_data1 = expand(os.path.join(config["dir"]["cut_base_dir"],"{sample}_cut_R1.fq.gz"),sample=samples),
        cut_base_data2 = expand(os.path.join(config["dir"]["cut_base_dir"],"{sample}_cut_R2.fq.gz"),sample=samples),
        cut_data_fastqc1=expand(os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_cut_R1_fastqc.html"),sample=samples),
        cut_data_fastqc2=expand(os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_cut_R2_fastqc.html"),sample=samples),
        sort_bam_file = expand(os.path.join( config["dir"]["sort_bam_dir"],"{sample}.sort.bam" ), sample=samples),
        sort_bam_bai = expand(os.path.join( config["dir"]["sort_bam_dir"],"{sample}.sort.bam.bai" ), sample=samples),
        quantify = expand(os.path.join( config["dir"]["quantify_dir"],"{sample}/{sample}.gtf"), sample=samples)


# 01 raw data quality control
rule raw_fastqc:
    input:
         raw_data = os.path.join(config["dir"]["raw_data_dir"],"{sample}_{R}.fq.gz")
    output:
         raw_data_fastqc = os.path.join(config["dir"]["raw_data_qc_result_dir"],"{sample}_{R}_fastqc.html")
    params:
         raw_data_fastqc_dir = os.path.join(config["dir"]["raw_data_qc_result_dir"])
    message: """----------------- Quality check of raw data with Fastqc.--------------"""
    shell:
         "fastqc {input.raw_data} -o {params.raw_data_fastqc_dir} "

# 02 filter low-quality sequences
rule trim:
    input:
         raw_data1 = os.path.join(config["dir"]["raw_data_dir"],"{sample}_R1.fq.gz"),
         raw_data2 =os.path.join(config["dir"]["raw_data_dir"],"{sample}_R2.fq.gz")
    output:
         clean1 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_R1_val_1.fq.gz"),
         clean2 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_R2_val_2.fq.gz")
    message: """-------------------------- Trim. ------------------------------------"""
    params:
         quality = config["para"]["trim_quality"],   #碱基最小质量分数
         stringency = config["para"]["trim_stringency"], #设定可以忍受的前后adapter重叠的碱基数,默认为1
         length = config["para"]["trim_length"],    #设定输出reads长度阈值,小于设定值会被抛弃
         output_dir = config["dir"]["clean_data_dir"]
    shell:
         "trim_galore --gzip --phred33 -q {params.quality} -j 16\
         --stringency {params.stringency} --length {params.length} --paired {input[0]} {input[1]} -o {params.output_dir}"

# 03 clean data quality control
rule trim_fastqc:
     input:
          clean_data1 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_R1_val_1.fq.gz"),
          clean_data2 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_R2_val_2.fq.gz")
     output:
          clean_data_fastqc1=os.path.join(config["dir"]["clean_data_qc_result_dir"],"{sample}_R1_val_1_fastqc.html"),
          clean_data_fastqc2=os.path.join(config["dir"]["clean_data_qc_result_dir"],"{sample}_R2_val_2_fastqc.html")
     params:
          os.path.join(config["dir"]["clean_data_qc_result_dir"])
     message: """------------ Quality check of clean data with Fastqc.---------------"""
     shell:
          "fastqc {input.clean_data1} -o {params[0]} && fastqc {input.clean_data2} -o {params[0]}"

# 04 cut base 
rule cut_base:
     input:
          clean_data1 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_R1_val_1.fq.gz"),
          clean_data2 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_R2_val_2.fq.gz")
     output:
          cut_base_data1 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_cut_R1.fq.gz"),
          cut_base_data2 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_cut_R2.fq.gz")
     message: """----------------------- Cut base ---------------------------------------."""
     params:
         left_base = config["para"]["cut_left_base"],
         right_base = config["para"]["cut_right_base"]
     shell:
          "cutadapt -u {params.left_base} -U {params.right_base} -o {output[0]} -p {output[1]} {input[0]} {input[1]}"

#05 cut fastqc
rule cut_fastqc:
     input:
          cut_base_data1 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_cut_R1.fq.gz"),
          cut_base_data2 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_cut_R2.fq.gz")
     output:
          cut_data_fastqc1=os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_cut_R1_fastqc.html"),
          cut_data_fastqc2=os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_cut_R2_fastqc.html")
     params:
          os.path.join(config["dir"]["cut_data_qc_result_dir"])
     message: """------------ Cut base data with Fastqc.---------------"""
     shell:
          "fastqc {input.cut_base_data1} -o {params[0]} && fastqc {input.cut_base_data2} -o {params[0]} "

# 06 align to reference and sort bam file
rule align_sort_bam:
     input:
         cut_base_data1 = os.path.join(config["dir"]["cut_base_dir"], "{sample}_cut_R1.fq.gz"),
         cut_base_data2 = os.path.join(config["dir"]["cut_base_dir"], "{sample}_cut_R2.fq.gz")
     output:
         sort_bam_file = os.path.join(config["dir"]["sort_bam_dir"], "{sample}.sort.bam")
     log:
         align_log = os.path.join(config["dir"]["bam_qc_dir"], "{sample}.sort.bam.log")
     message: """---------------------- Alignment ---------------------------------------."""
     params:
         reference = config["reference"]["mouse_reference"]
     shell:
         "hisat2 -p 12 --summary-file {log.align_log} -x {params.reference} -1 {input.cut_base_data1} -2 {input.cut_base_data2} \
         |samtools sort -O BAM -@ 16 -o {output.sort_bam_file}"

# 07 construct sort bam file index
rule sort_bam_index:
    input:
        sort_bam = os.path.join(config["dir"]["sort_bam_dir"],"{sample}.sort.bam")
    output:
        sort_bam_bai = os.path.join(config["dir"]["sort_bam_dir"],"{sample}.sort.bam.bai")
    shell:
         "samtools index {input.sort_bam} {output.sort_bam_bai}"

# 08 quantify
rule quantify:
     input:
         sort_bam = os.path.join(config["dir"]["sort_bam_dir"],"{sample}.sort.bam")
     output:
         quantify = os.path.join(config["dir"]["quantify_dir"],"{sample}/{sample}.gtf")
     message: """------------- Quantify-------------."""
     params:
          annotation = config["annotation"]["mouse_gtf"]
     shell:
         "stringtie {input.sort_bam} -p 8 -e -B -G {params.annotation} -o {output.quantify}"


print("上游分析已经完成,记得去对样本表达定量成表达谱")