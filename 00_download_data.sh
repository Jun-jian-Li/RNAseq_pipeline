cd 01_raw_data

grep "SRR*" SRR_Acc_List.txt|while read id;
do
        kingfisher get -r $id -m aws-http ena-ftp prefetch --download-threads 16
done

cd ../
