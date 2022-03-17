
#创建AS的分析conda环境

# conda create -n AS python=2
# conda install -y rmats
# conda install -y rmats2sashimiplot #rmats的可视化软件
# RNASeq-MATS.py -h # 测试安装成功

#这里需要设置自己的样本分组
echo Treat1.bam,Treat2.bam > 07_alternative_splicing/group1.txt
echo Control1.bam,Control2.bam > 07_alternative_splicing/group2.txt

# -t paired 测序类型
# --readLength 140 reads的长度

RNASeq-MATS.py --b1 07_alternative_splicing/group1.txt \
    --b2 07_alternative_splicing/group2.txt \
    --gtf /home/ljj/Genome/mouse_annotation/Mus_musculus.GRCm38.99.chr.gtf \
    --od 07_alternative_splicing/ \
    -t paired \
    --nthread 10 \
    --cstat 0.0001 \
    --readLength 140 \
    --tmp ../tmp