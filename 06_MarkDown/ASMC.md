# ASMC
> P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. Nature Genetics, 2018. 
### Summary (TL;DR)
ASMC是一种有效估计基因组`pairwise coalescence times`的方法。它可以使用SNP芯片或全基因组测序(WGS)数据运行。
### 输入文件格式
* HAP文件
是一个没有标题行的文本文件，有 2N+5 或 2N 列，其中 N 是样本数。 在前一种情况下，前五列是：
  1. 染色体代码
  2. 变体 ID
  3. 碱基对坐标 等位基因 0（通常是次要的，导入时使用 'ref-first' 视为 REF）
  4. 等位基因 1（通常是主要的，导入时使用 'ref-last' 视为 REF） 
接下来是第一个样本的一对 0/1 值的单倍型列，然后是第二个样本的一对单倍型列
* Sample文件
  gen 或 .bgen 基因型剂量文件或 .haps 分阶段参考面板随附的样本信息文件。 加载了--data/--sample，并在某些情况下由--export 生成。

默认情况下，由 --export 发出的 .sample 空格分隔文件有两个标题行，然后每个样本一行有 4 个以上的字段：


* map文件


