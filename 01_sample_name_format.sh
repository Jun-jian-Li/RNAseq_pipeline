#!/usr/bash
echo "#----该文件是整理公共数据下载下的数据文件以及自测数据的文件命名----#";
echo "----请按照以下格式输入 bash file.sh public/self single/pair----"
echo "参数1:public测序或者是self测序-->$1";
echo "参数2:single测序或者是pair测序-->$2";

if [ $1 == "public" ];then    #这个后面一定有个空格
echo "#----你选择的数据是public数据库数据----#";
    if [ $2 == "single" ];then
    echo "#----你选择的数据是single数据----#";
    #对原始文件进行压缩
        cd ./01_raw_data
        ls *.fastq|while read id
        do
        pigz $id -p 16
        done
        #返回路径
        cd ..
    #对压缩后的文件进行重命名
    #单端测序 修改fastq的名字
        cd ./01_raw_data
        ls |grep '.gz'|while read id
        do
        newname1=$(echo $id |sed 's/ast//g')
        mv $id $newname1
        done

    elif [ $2 == "pair" ];then
    echo "#----你选择的数据是pair数据----#";
    #对原始文件进行压缩
        cd ./01_raw_data
        ls *.fastq|while read id
        do
        pigz $id -p 16
        done
        #返回路径
        cd ..
    #对压缩后的文件进行重命名
    #单端测序 修改fastq的名字
        cd ./01_raw_data
        #替换fastq文件名字 ,双端测序的数据
        ls |grep '.gz'|while read id
        do
        newname=$(echo $id |sed 's/_/_R/g')
        newname1=$(echo $newname |sed 's/ast//g')
        mv $id $newname1
        done
    fi
fi
    #----------------------------------------------------------------
    #自测数据集
if [ $1 == "self" ];then    #这个后面一定有个空格
    echo "#----你选择的数据是self数据库数据----#";
        if [ $2 == "single" ];then
        echo "#----你选择的数据是single数据----#";
        #单端测序 修改fastq的名字
            cd ./01_raw_data
            ls |grep '.gz'|while read id
            do
            newname1=$(echo $id |sed 's/ast//g')
            mv $id $newname1
            done
        #双端测序
        elif [ $2 == "pair" ];then
        echo "#----你选择的数据是pair数据----#";
        #对压缩后的文件进行重命名
        #单端测序 修改fastq的名字
            cd ./01_raw_data
            #替换fastq文件名字 ,双端测序的数据
            ls |grep '.gz'|while read id
            do
            newname=$(echo $id |sed 's/_/_R/g')
            newname1=$(echo $newname |sed 's/ast//g')
            mv $id $newname1
            done
        fi
fi



echo "请选择对应测序类型以及物种进行运行脚本:如物种:人类,测序类型:单端 -> 运行02_snakemake_human_single.py文件"
