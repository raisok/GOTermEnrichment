**这是一个做富集分析的python包**
* 需要用到的python环境包有
    * pandas,numpy,re,logging,os,sys,statsmodels.sandbox.stats.multicomp
* 需要用到的文件
    * GO数据库里面的obo文件
    * 差异表达基因的ID文件
    * 基因的GO注释文件
* 用法
    * ``python GOTermEnrichment.py --oboFile gene_ontology.1_2.obo --selectionSetFiles gene.lst --gene2GoFile hg19.annot --outdir ./``
