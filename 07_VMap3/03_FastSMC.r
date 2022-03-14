#66:/data1/home/yafei/004_Vmap3/FastSMC/input
awk '{if(NR>2) print $2,$2,$3,$4;else print $0}' Alineage001_Free_threshing_phased_imputed.sample