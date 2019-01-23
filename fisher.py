#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: fisher.py
@time: 2019/01/08
'''

from scipy.stats import fisher_exact, hypergeom
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import re
import logging
import time

'''
    超几何分布的方法
    GO example:
    extracellular region    61 in 715 DEGs, 195 in 4564 all GO genes

                        差异      非差异
    属于某一个GO term    100     1000
    不属于某一个GO term   500     10000

    表示总共有11600个基因,其中600个位差异基因，属于某一个GO注释目录的有1100，其中属于差异基因的有100个
    x=100,m=1100,n=10500,k=600
    x,m,n,k

    x = vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
    m	the number of white balls in the urn.
    n	the number of black balls in the urn.
    k	the number of balls drawn from the urn.


    1- phyper(100-1, 1100, 10500, 600)

    x,m+n,m,k
    stats.hypergeom.sf(99,11600,1100,600)

    extracellular region    61 in 715 DEGs, 195 in 4564 all GO genes
    x=60,m+n=4564,m=195 ,k=715
    输入依次为61个差异基因属于某一个GO term，总共有4564个基因有GO注释，195个属于某一个GO term，我们的DEG有715个


    fisher检验的方法
                path    nonpath
interest          a         b   num_interestgene
non interest      c         d   num_allgene_num_interestgene
                num_sp   num_allgene-num_sp

a   为在通路里面的感兴趣的基因   61
b   为不在通路里面的感兴趣的基因  715-61
c   为在通路里面不属于感兴趣基因  195-61
d   为不在通路里面不感兴趣的基因  4564-195

'''


def read_fg(file):
    dic_fg={}
    with open(file,'r') as fr:
        for line in fr.readlines():
            line = line.strip()
            dic_fg[line]=1
    return dic_fg

def read_bg(file):
    dic_bg={}
    allgene=[]
    with open(file,'r') as fr:
        for line in fr.readlines():
            line = line.strip()
            goid=line.split('\t')[1]
            geneid=line.split('\t')[0]
            allgene.append(geneid)
            if goid in dic_bg:
                dic_bg[goid].append(geneid)
            else:
                dic_bg[goid]=[geneid]
    allgene_num=len(set(allgene))
    Gene2GO,GOnum=GO2GeneToGene2GO(dic_bg)
    logging.warning("GeneNum:\t"+str(allgene_num))
    logging.warning("GOnum:\t"+str(GOnum))
    return dic_bg,allgene_num

def GO2GeneToGene2GO(dic_bg):
    Gene2GO={}
    GOnum=len(dic_bg.keys())
    for iv,ik in dic_bg.items():
        for gene in ik:
            if gene in Gene2GO:
                Gene2GO[gene].append(iv)
            else:
                Gene2GO[gene]=[iv]
    return Gene2GO,GOnum

def Gene2GoToGo2Gene(Gene2Go):
    Go2Gene={}
    for Gene,Gos in Gene2Go.items():
        for go in Gos:
            if go in Go2Gene:
                Go2Gene[go].append(Gene)
            else:
                Go2Gene[go]=[Gene]
    return Go2Gene

def parseAllGO():
    parsedGene2Go={}
    parsedGo2Gene={}

    return parsedGene2Go,parsedGo2Gene

def get_num(dic_bg,dic_fg):
    fw=open('go.lst','w')
    for iv,ik in dic_bg.items():
        fg_count=0
        for gl in ik:
            bg_num=len(ik)
            if gl in dic_fg:
                fg_count +=1
        print(iv+"\t"+str(fg_count)+"\t"+str(bg_num))
        fw.write(iv+"\t"+str(fg_count)+"\t"+str(bg_num)+"\n")
    fw.close()

def count_gene(dic_fg,dic_bg):
    alldegnum=[]
    dicgo={}
    deg_in_go={}
    for goid,genel in dic_bg.items():
        allgene_in_term_num=len(genel)
        for gene in genel:
            if gene in dic_fg:
                if goid in deg_in_go:
                    deg_in_go[goid].append(gene)
                else:
                    deg_in_go[goid]=[]
                    deg_in_go[goid].append(gene)
                alldegnum.append(gene)
                deg_in_term_num=len(deg_in_go[goid])
                dicgo[goid]=[deg_in_term_num,allgene_in_term_num]
    alldeg_in_go_num=len(set(alldegnum))
    return alldeg_in_go_num,deg_in_go,dicgo

def getParents(id,parents,results):
    if id not in parents.keys():
        return
    for iv in parents[id]:
        if iv not in parents.keys():
            results[id] = 1
        else:
            getParents(iv,parents,results)

def read_go_obo(file):
    logging.warning(time.asctime()+"\tread file:"+file)
    id_regex01 = re.compile('^id:\s+(GO:.*)')
    name_regex02 = re.compile('^name:\s+(.*)')
    namespace_regex03 = re.compile('^namespace:\s+(.*)')
    altid_regex05 = re.compile('^alt_id:\s+(GO:.*)')
    isa_regex06 = re.compile('^is_a:\s+(GO:[^\s]+)')
    relation_regex07 = re.compile('^relationship:.*(GO:[\d]+)')
    parents={}
    alias={}
    gos={}
    obsoletes={}
    id=""
    name=""
    with open (file,'r') as fr:
        for line in fr.readlines():
            line = line.strip()
            if line.startswith("[Typedef]"):
                break
            if not line.startswith("[Term]"):
                if line.startswith("id"):
                    id=re.findall(id_regex01, line)[0]
                elif line.startswith("name:"):
                    name=re.findall(name_regex02, line)[0]
                elif line.startswith("namespace"):
                    namespace=re.findall(namespace_regex03, line)[0]
                    gos[id] = [name, namespace]
                elif line.startswith("is_obsolete"):
                    obsoletes[id] = 1
                elif line.startswith("alt_id"):
                    if id in alias:
                        alias[id].append(re.findall(altid_regex05, line)[0])
                    else:
                        alias[id] = [re.findall(altid_regex05, line)[0]]
                elif line.startswith("is_a"):
                    if id in parents:
                        parents[id].append(re.findall(isa_regex06, line)[0])
                    else:
                        parents[id] = [re.findall(isa_regex06, line)[0]]
                elif line.startswith("relationship"):
                    if id in parents:
                        parents[id].append(re.findall(relation_regex07, line)[0])
                    else:
                        parents[id] = [re.findall(relation_regex07, line)[0]]
    return parents,alias,gos,obsoletes

def format_go_result(dicgo,deg_in_go,alldeg_in_go_num,allgene,gos):
    dic_hy={}
    for go,num in dicgo.items():
        deg_in_term_num=num[0]
        a=deg_in_term_num-1
        allgene_in_term_num=num[1]
        if go in deg_in_go and go in gos:
            pvalue = hypergeom.sf(a, allgene, allgene_in_term_num, alldeg_in_go_num)
            # print(go+"\t"+gos[go][0]+"\t"+str(deg_in_term_num)+"of "+str(alldeg_in_go_num)+" in the list\t"+str(allgene_in_term_num)+" of "+str(allgene)+" in the genome\t"+str(pvalue)+"\t"+";".join(deg_in_go[go]))
            dic_hy[go]=[deg_in_term_num, allgene, allgene_in_term_num, alldeg_in_go_num, pvalue,gos[go][0],";".join(deg_in_go[go])]
    return dic_hy

def read_enrichfile_hypergeom(file):
    '''
    这里用来存GO注释的结果，包含了差异表达基因的GO注释和所有背景基因的GO注释
    hypergeom.sf(deg_in_special_go-1,allgene_has_go,allgene_in_special_go,alldeg_has_go)
        超几何分布的方法
    GO example:
    extracellular region    61 in 715 DEGs, 195 in 4564 all GO genes
    这里作为例子的输入依次为61个差异基因属于某一个GO term，总共有4564个基因有GO注释，195个属于某一个GO term，我们的DEG有715个

                        差异      非差异
    属于某一个GO term    100     1000
    不属于某一个GO term   500     10000

    表示总共有11600个基因,其中600个位差异基因，属于某一个GO注释目录的有1100，其中属于差异基因的有100个
    x=100,m=1100,n=10500,k=600

    x,m+n,m,k
    stats.hypergeom.sf(99,11600,1100,600)

    '''

    dic_hy = {}
    with open(file, 'r') as go:
        for line in go.readlines():
            line = line.strip()
            goid = line.split('\t')[0]
            deg_in_special_go = int(line.split('\t')[2].split(' ')[0])
            a = deg_in_special_go - 1
            alldeg_has_go = int(line.split('\t')[2].split(' ')[2])
            allgene_has_go = int(line.split('\t')[3].split(' ')[2])
            allgene_in_special_go = int(line.split('\t')[3].split(' ')[0])
            qvalue = line.split('\t')[4]
            pvalue = hypergeom.sf(a, allgene_has_go, allgene_in_special_go, alldeg_has_go)
            dic_hy[goid] = [deg_in_special_go, allgene_has_go, allgene_in_special_go, alldeg_has_go, pvalue, qvalue]
    return dic_hy


def read_enrichfile_fisher(file):
    dic_go = {}
    with open(file, 'r') as go:
        for line in go.readlines():
            line = line.strip()
            goid = line.split('\t')[0]
            genlis = line.split('\t')[2]
            a = genlis.split(' ')[0]
            b = int(genlis.split(' ')[2]) - int(a)
            allgenelis = line.split('\t')[3]
            c = int(allgenelis.split(' ')[0]) - int(a)
            d = int(allgenelis.split(' ')[2]) - int(allgenelis.split(' ')[0])
            qvalue = line.split('\t')[4]
            pvalue = fisher_exact([[a, b], [c, d]])[1]
            dic_go[goid] = [a, b, c, d, pvalue, qvalue]
    return dic_go


def p_adjust_bh(p):
    '''
        这里使用BH的方法进行p值校正，华大的p值校正是使用bonferroni的方法
        将一系列的p值按照从大到小排序，然后利用下述公式计算每个p值所对应的FDR值
        p是这一次检验的pvalue，n是检验的次数，i是排序后的位置ID（如最大的P值的i值肯定为n，第二大则是n-1，依次至最小为1
        主要公式为p * (n/i)
        将计算出来的FDR值赋予给排序后的p值，如果某一个p值所对应的FDR值大于前一位p值（排序的前一位）所对应的FDR值，
        则放弃公式计算出来的FDR值，选用与它前一位相同的值。因此会产生连续相同FDR值的现象；反之则保留计算的FDR值
        将FDR值按照最初始的p值的顺序进行重新排序，返回结果
    '''
    p = np.asfarray(p)  # 将列表元素的类型变为浮点型
    by_descend = p.argsort()[::-1]  # 将元素从大到小进行排序，并将索引存起来
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)  # 利用位置除以p值列表的数目，来构建BH检验的公式
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def p_adjust_bonferroni(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    qvalue = multipletests(p, method='bonferroni')
    return qvalue[1][by_orig], qvalue[0][by_orig]


def get_output(dic_hy):
    data = pd.DataFrame.from_dict(dic_hy, orient="columns", dtype=None)
    dataout = pd.DataFrame.transpose(data)
    dataout.columns = ["HitsGenesCountsInSelectedSet", "AllGenesCountsInBackground", "HitsGenesCountsInBackground", "AllGenesCountsInSelectedSet", "pValue", "GO_Name","GenesOfSelectedSetInGOterm"]
    datasort = dataout.sort_values(by="pValue", axis=0)
    pValue = datasort["pValue"]
    fdr = p_adjust_bh(pValue)
    fdr = pd.DataFrame(fdr, index=datasort.index)
    datasort.insert(5, "q-value(BH method)", fdr)
    datasort.to_csv('goenrich.xls', encoding='utf-8', index=True, sep="\t",index_label="GO_ID")


if __name__ == '__main__':
    file = 'HBRR-VS-UHRR.DEseq2_Method_C.xls'
    # dic_go=read_enrichfile_fisher(file)

    dic_hy = read_enrichfile_hypergeom(file)
    get_output(dic_hy)

    # obo = 'gene_ontology.1_2.obo'
    # deglis='fg.lst'
    # dic_fg=read_fg(deglis)
    # c_annot='hg19.c.annot'
    # dic_bg, allgene_num=read_bg(c_annot)
    #
    # get_num(dic_bg,dic_fg)
    # parents, alias, gos, obsoletes=read_go_obo(obo)
    #
    # alldeg_in_go_num, deg_in_go, dicgo = count_gene(dic_fg,dic_bg)
    # # print(deg_in_go)
    # dic_hy=format_go_result(dicgo,deg_in_go,alldeg_in_go_num,allgene_num,gos)
    #
    # get_output(dic_hy)
