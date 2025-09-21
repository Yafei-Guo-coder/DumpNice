# 172.24.116.24
cat sample.txt | while read sample dir lib ref row col
do
cd ${dir}
/work_home/chenchuan/soft/saw-8.0.2/bin/saw convert gef2gem --gef ${sample}.raw.gef --gem ${sample}.gem --bin-size 1
sed '1,8d' ${sample}.gem |awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6}' | gzip > ${sample}.gem.gz
mkdir ${sample}
gem_ext part -c "#932d22" -b -s ${sample} -p /home/chenchuan/yafei/spatio/${dir}/${sample}_DAPI_regist_ps.tif -o ${sample} -f /home/chenchuan/yafei/spatio/${dir}/${sample}_DAPI_regist_ps.tif /home/chenchuan/yafei/spatio/${dir}/${sample}.gem.gz ${lib}
cd ${sample}
gem_ext merge ${sample}.info.txt -s ${sample} ${row} ${col}
#sed 's///' /home/chenchuan/yafei/spatio/stomics_bin_slice_cc.R
Rscript /home/chenchuan/yafei/spatio/stomics_bin_slice_cc.R ${sample} /home/chenchuan/yafei/spatio/${dir}/${sample}/${sample}.merged.gem.gz 50 0 30 0.5
done

# 转格式
/work_home/chenchuan/soft/saw-8.0.2/bin/saw convert gef2gem --gef D05507C4.raw.gef --gem D05507C4.gem --bin-size 1
sed '1,8d' D05507C4.gem | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6}' | gzip > D05507C4.gem.gz

mkdir D05507C4
gem_ext part -c "#cc1c67" -b -s D05507C4 -p D05507C4_DAPI_regist_ps.tif -o D05507C4 -f D05507C4_DAPI_regist_ps.tif D05507C4.gem.gz 6

cd D05507C4
gem_ext merge D05507C4.info.txt -s D05507C4 2 3

STO_mtx2pic D05936C1.gem.gz 1 bin1 -u -m 1

Rscript /data/chenchuan/4.Soybean/2.Script/plot.spatial.cluster.R /home/guoyafei1/01_soybean_spatial/input/merge/cd_07/cd_07_bin50/cd_07_bin50_bin50_cluster.txt /home/guoyafei1/01_soybean_spatial/input/merge/cd_07/cd_07_bin50/cd_07_bin50_bin50_col.txt 50 /home/guoyafei1/01_soybean_spatial/input/merge/cd_07/cd_07_bin50/cd_07_bin50

# R 有个包出错
export LD_LIBRARY_PATH=/home/guoyafei1/software/anaconda3/envs/r_env_icuv58/lib:$LD_LIBRARY_PATH
gem_ext part -c "#932d22" -b -s D05936C1 -p /home/guoyafei1/yafei/spatio/D05936C1-CD_Nod_21_dpi/D05936C1_DAPI_regist_ps.tif -o D05936C1 -f /home/guoyafei1/yafei/spatio/D05936C1-CD_Nod_21_dpi/D05936C1_DAPI_regist_ps.tif /home/guoyafei1/yafei/spatio/D05936C1-CD_Nod_21_dpi/D05936C1.gem.gz 6

STO_mtx2pic D05936C1.gem.gz 1 bin1 -u -m 1

# 合并文件
# 大豆7天
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.mask_with_mark.tif
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.mask_with_mark.tif
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.mask_with_mark.tif
cat /home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.info.txt /home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.info.txt /home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.info.txt
# info.txt-1
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.6.gem.gz	0	0	1	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.6.part.tif	4	3	0	-2
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.4.gem.gz	0	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.4.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.5.gem.gz	1	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.5.part.tif	0	2	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.1.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.1.part.tif	0	2	0	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.2.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.2.part.tif	0	1	-1	-2
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.3.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.3.part.tif	0	0	-1	0
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.4.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.4.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.5.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.5.part.tif	0	0	-1	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.6.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.6.part.tif	0	1	-1	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.1.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.1.part.tif	2	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.2.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.2.part.tif	0	1	0	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.3.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.3.part.tif	0	1	-1	-2
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.8.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.8.part.tif	0	1	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.7.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.7.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.5.part.tif	1	0	-2	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.6.part.tif	0	1	-2	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.3.part.tif	2	1	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.1.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.2.part.tif	1	0	-2	-1
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.4.part.tif	1	0	-1	-1
# info.txt-2
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.1.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.1.part.tif	0	2	0	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.2.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.2.part.tif	0	1	-1	-2
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.3.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.3.part.tif	0	0	-1	0
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.6.gem.gz	0	0	1	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.6.part.tif	4	3	0	-2
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.4.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.4.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.5.gem.gz	1	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507F2/D05507F2/D05507F2.5.part.tif	0	2	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.1.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.1.part.tif	2	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.2.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.2.part.tif	0	1	0	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.3.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.3.part.tif	0	1	-1	-2
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.4.gem.gz	1	1	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.4.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.5.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.5.part.tif	0	0	-1	-1
/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.6.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/D05507C6/D05507C6/D05507C6.6.part.tif	0	1	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.3.part.tif	2	1	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.1.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.2.part.tif	1	0	-2	-1
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.4.part.tif	1	0	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.8.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.8.part.tif	0	1	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.7.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.7.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.5.part.tif	1	0	-2	0
/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507E5/D05507E5/D05507E5.6.part.tif	0	1	-2	0

# 大豆10天
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.mask_with_mark.tif
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.mask_with_mark.tif
cat /home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.info.txt /home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.info.txt
# info.txt-1
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.6.part.tif	1	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.4.part.tif	0	2	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.5.part.tif	0	0	-3	-1
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.3.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.3.part.tif	1	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.1.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.1.part.tif	0	2	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.2.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.2.part.tif	0	2	0	-3
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.4.part.tif	0	0	-3	-2
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.6.part.tif	0	0	0	-2
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.5.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.1.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.1.part.tif	0	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.2.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.2.part.tif	0	1	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.3.part.tif	0	1	0	0
# info.txt-2
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.2.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.2.part.tif	1	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.1.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.1.part.tif	0	2	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.3.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.3.part.tif	0	2	0	-3
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.5.part.tif	1	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.4.part.tif	0	2	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507A3/D05507A3/D05507A3.6.part.tif	0	0	-3	-1
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.1.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.1.part.tif	0	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.3.gem.gz	1	1	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.3.part.tif	0	1	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.2.part.tif	0	1	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.4.part.tif	0	0	-3	-2
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.6.part.tif	0	0	0	-2
/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507F6/D05507F6/D05507F6.5.part.tif	0	0	0	0

# 大豆21天
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.mask_with_mark.tif
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.mask_with_mark.tif
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.mask_with_mark.tif
cat /home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.info.txt /home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.info.txt /home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.info.txt
# info.txt-1
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.6.part.tif	2	2	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.5.part.tif	0	0	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.4.part.tif	1	3	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.3.part.tif	3	2	-3	0
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.2.part.tif	1	5	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.1.part.tif	0	0	-6	-2
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.5.part.tif	0	2	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.4.part.tif	1	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.6.part.tif	0	0	-4	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.3.part.tif	1	3	0	-2
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.2.part.tif	1	1	-2	-1
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.1.part.tif	1	12	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.3.part.tif	1	0	-5	-2
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.4.part.tif	0	1	-2	-1
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.1.part.tif	15	20	0	-2
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.2.part.tif	0	15	0	0
# info.txt-2
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.2.part.tif	1	5	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.3.part.tif	3	2	-3	0
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.1.part.tif	0	0	-6	-2
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.6.part.tif	2	2	0	0
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.5.part.tif	0	0	-1	-1
/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05507C4/D05507C4/D05507C4.4.part.tif	1	3	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.2.part.tif	1	1	-2	-1
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.3.part.tif	1	3	0	-2
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.1.part.tif	1	12	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.5.part.tif	0	2	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.4.part.tif	1	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936D2/D05936D2/D05936D2.6.part.tif	0	0	-4	0
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.1.part.tif	15	20	0	-2
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.2.part.tif	0	15	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.3.part.tif	1	0	-5	-2
/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936C1/D05936C1/D05936C1.4.part.tif	0	1	-2	-1

# 苜蓿7天
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.mask_with_mark.tif
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.mask_with_mark.tif
cat /home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.info.txt /home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.info.txt
# info.txt-1
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.5.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.5.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.4.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.4.part.tif	0	0	-1	0
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.6.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.6.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.2.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.2.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.3.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.3.part.tif	0	0	0	0
/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.1.gem.gz	0	0	0	/home/xiaomanyi/1_soybean_spatial/input/C05912C6_A17/C05912C6/C05912C6.1.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.5.part.tif	1	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.6.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.4.part.tif	0	1	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.7.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.7.part.tif	0	1	-2	-1
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.1.part.tif	0	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.2.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.2.part.tif	0	0	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.3.part.tif	0	0	0	-1
# info.txt-2
/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.2.gem.gz	0	1	0	/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.2.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.1.gem.gz	0	1	0	/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.1.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.4.part.tif	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.3.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.3.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/C05912C6/C05912C6/C05912C6.5.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.1.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.1.part.tif	0	-1	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.2.gem.gz	0	1	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.2.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.3.gem.gz	0	1	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.3.part.tif	0	-1
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.5.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.5.part.tif	1	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.6.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.6.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.4.gem.gz	0	0	0	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.4.part.tif	0	0
/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.7.gem.gz	0	0	1	/home/guoyafei1/01_soybean_spatial/input/D05936A3/D05936A3/D05936A3.7.part.tif	0	-2	-1

/home/guoyafei1/01_soybean_spatial/input/merge/cd_07/cd_07.merged.ssDNA.tif
/home/guoyafei1/01_soybean_spatial/input/merge/cd_07/cd_07_bin50/cd_07_bin50_bin50_cluster.pdf
/home/guoyafei1/01_soybean_spatial/input/merge/cd_10/cd_10.merged.ssDNA.tif
/home/guoyafei1/01_soybean_spatial/input/merge/cd_10/cd_10_bin50/cd_10_bin50_bin50_cluster.pdf
/home/guoyafei1/01_soybean_spatial/input/merge/cd_21/cd_21.merged.ssDNA.tif
/home/guoyafei1/01_soybean_spatial/input/merge/cd_21/cd_21_bin50/cd_21_bin50_bin50_cluster.pdf
/home/guoyafei1/01_soybean_spatial/input/merge/mt_07/mt_07.merged.ssDNA.tif
/home/guoyafei1/01_soybean_spatial/input/merge/mt_07/mt_07_bin50/mt_07_bin50_bin50_cluster.pdf

LD_PRELOAD=/usr/lib64/libstdc++.so.6 LD_LIBRARY_PATH=/home/guoyafei1/software/anaconda3/lib:$LD_LIBRARY_PATH R








