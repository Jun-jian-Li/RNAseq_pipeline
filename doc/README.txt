#-------------------------------------------------------------------------
#流程说明:RNAseq数据的处理及差异分析流程
#作者:Junjian Li
#描述:本流程基于fastqc trim-galore cutadapter(剪切前端质量不好的碱基) 
				 hisat2 stringtie edgeR流程进行RNAseq数据的处理及差异分析
#-------------------------------------------------------------------------

#文件整体流程的分析说明及使用
01_环境搭建:搭建linux系统的高通量测序环境
    1:基于系统Centos7.9稳定版
	2:基于环境conda环境
		2.1:下载miniconda 清华镜像源:https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/
		     下载最新的 Miniconda https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
			 bash Miniconda3-latest-Linux-x86_64.sh
			 此时注意,在个人主目录下安装conda,尽量选择安装的位置在/home/用户/software 的目录下,方便管理
             
			 随后在software的环境下将miniconda 中的bin加入环境变量中 ~/.bahsrc
			 expert PATH='/home/zs/software/miniconda3/bin:$PATH'
			 激活环境变量
			 source ~/.bashrc
			 
		2.2:#添加国内软件下载镜像
			 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
			 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
			 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
			 conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/

			 conda config --set show_channel_urls yes            #显示已经安装的频道
			 conda config --get channels                             #查看安装的频道
		
		2.3:创建conda的分析环境,注意写明python需要的环境
			 conda activate -n rnaseq python=3.9
			
			 下载mamba的国内软件管理源
			 conda install mamba 
			 mamba install hisat2 sra-tools trim-galore stringtie samtools bwa snakemake -Y 
			
			 此时激活conda的环境 
		     conda activate rnaseq
				
	3:基于snakemake环境搭建
			 整体rnaseq项目是以snakemake工具搭建的流程化的项目,我们需要对项目创建相应的文件系统
			 env:conda的环境以及文件系统的环境 
				conda env export > env/environments.yaml 
				conda env create -f py36.yaml
			 doc:存放说明文档
			 result:结果文件目录
			 02*.smk:整个项目的代码块
			 03_extract_exp: 调用提取表达结果的shell脚本
			 00_script:存放提取表达的python脚本
02_常规流程:基于已有的基因注释信息对基因,转录本进行组装定量分析
			 分析物种:现目前存放了小鼠以及人类
			 分析流程:fastqc->trim_galore(fastqc)->cutadapter(剪切前端质量不好的碱基+fastqc)->hisat2 ->samtools (sort index)->stringtie ->Extract_exp
			 reference目录下:参考基因组选取ensembl的GRCm38 (构建了snp trans索引的版本)
			 annotation目录下:参考基因组注释文件选择 GRCm38 
			 
			运行操作:
			1:激活环境
			 1.1:conda环境:conda activate rnaseq
			 
		    2:运行代码
			 2.1: 01 理下载以及自己测序文件格式,运行规则:bash 01_sample_name_format.sh public/self single/pair
			 2.2: 02 snakemake -s snakemake.py -j 16 (-n :试运行一下) 文件整体流程的分析 选择对应物种以及测序类型进行计算
			 2.3: 03 获取表达谱以及count矩阵 : bash 03_extract_exp.sh + read的长度  
			 2.4: 04 获取质量控制文件:分别对原始,过滤,剪切,比对之后的文件进行质量控制,生成质控报告
			 2.5: 05 对分析之后的文件删除,释放储存空间



