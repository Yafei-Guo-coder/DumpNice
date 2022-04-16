#Raxml-ng
#yafei@204:/data2/yafei/004_Vmap3/Tree
raxml -f a -m GTRGAMMA --JC69 -p 12345 -x 12345 -# 100 -s ABsubgenome_01.phy -n ABsubgenome_01.raxml -o out -T 100.
raxml-ng --bootstrap --msa Alineage001_5000.phy --model GTR+G --prefix Alineage --seed 333 --threads 80 --bs-trees 100