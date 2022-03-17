
#对质量文件进行生成质控报告
multiqc result/01_raw_data_qc/*.zip -o ./result/01_raw_data_qc
multiqc result/02_clean_data_qc/*.zip -o ./result/02_clean_data_qc
multiqc result/03_cut_data_qc/*.zip -o ./result/03_cut_data_qc
multiqc result/04_bam_qc/*.log -o ./result/04_bam_qc
echo "质控报告完成,记得将生成文件下载，最后移除文件释放服务器储存空间!!!" ;
echo "运行05_copy_rm_files.sh"
