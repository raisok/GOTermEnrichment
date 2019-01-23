#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Enricher.py
@time: 2019/01/21
'''

from scipy.stats import fisher_exact, hypergeom


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


class Enricher(object):
    def __init__(self):
        self.enrich=0.0
        self.pValue=0.0
        self.allBalls=0
        self.allSelectBalls=0
        self.allWhiteBalls=0
        self.allSelectWhiteBalls=0

    def getEnrichment(self):
        return self.enrich

    def getPValue(self):
        return self.pValue

    def getAllBalls(self):
        return self.allBalls

    def setAllBalls(self,allBalls):
        self.allBalls=allBalls
        return self

    def getAllSelectBalls(self):
        return self.allSelectBalls

    def setAllSelectBalls(self,allSelectBalls):
        self.allSelectBalls=allSelectBalls
        return self

    def getAllWhiteBalls(self):
        return self.allWhiteBalls

    def setAllWhiteBalls(self,allWhiteBalls):
        self.allWhiteBalls = allWhiteBalls
        return self

    def getAllSelectWhiteBalls(self):
        return self.allSelectWhiteBalls

    def setAllSelectWhiteBalls(self,allSelectWhiteBalls):
        self.allSelectWhiteBalls=allSelectWhiteBalls
        return self

    def doHyperGeometric(self):
        self.HyperGeometric(self.allSelectWhiteBalls,self.allBalls,self.allWhiteBalls,self.allSelectBalls)

    def HyperGeometric(self,deg_in_special_go,allgene_has_go,allgene_in_special_go,alldeg_has_go):
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
        a = deg_in_special_go - 1
        self.pValue = float (hypergeom.sf(a, allgene_has_go, allgene_in_special_go, alldeg_has_go))
        self.enrich = float ((deg_in_special_go/alldeg_has_go)/(allgene_in_special_go/allgene_has_go))