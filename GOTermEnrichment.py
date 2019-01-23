#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GOTermEnrichment.py
@time: 2019/01/21
'''

import argparse
import sys
import logging
import time
import os
from GOLevelGetter import GOLevelGetter
from GObasicOboDatabase import GObasicOboDatabase
from Enricher import Enricher
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import logging



class GOTermEnrichment(object):

    def __init__(self):
        logging.warning("Using null creator... Remenber to Set DB file...")
        self.IntraColSeperator="\t"
        self.InterSeperator=","
        self.GoDb=None
        self.allGenes=set()
        self.selectGenes=set()
        self.Go2Genes={}
        self.parsedGo2Genes={}
        self.Gene2GO={}
        self.parsedGene2Go={}

    def getGoDb(self):
        return self.GoDb

    def setGoDb(self,GoDb):
        self.GoDb=GoDb

    def getAllGenes(self):
        return self.allGenes

    def setAllGenes(self,allGenes):
        self.allGenes=allGenes

    def getSelectGenes(self):
        return self.selectGenes

    def getGo2Genes(self):
        return self.Go2Genes

    def setGo2Genes(self,Go2Genes):
        self.Go2Genes=Go2Genes
        self.Go2GeneToGene2Go(Go2Genes,self.Gene2GO)

    def getParsedGo2Genes(self):
        return self.parseGo2Genes

    def setParsedGo2Genes(self,parsedGo2Genes):
        self.parseGo2Genes=parsedGo2Genes
        self.Go2GeneToGene2Go(parsedGo2Genes,self.parsedGene2Go)

    def getGene2Go(self):
        return self.Gene2GO

    def setGene2Go(self,Gene2Go):
        self.Gene2GO=Gene2Go

    def getGene2Go(self):
        return self.Gene2GO

    def setGene2Go(self,Gene2Go):
        self.Gene2GO=Gene2Go
        self.Gene2GoToGo2Gene(Gene2Go,self.Go2Genes)

    def getParsedGene2Go(self):
        return self.getParsedGene2Go

    def setParsedGene2Go(self,parsedGene2Go):
        self.parsedGene2Go=parsedGene2Go
        self.Gene2GoToGo2Gene(parsedGene2Go,self.parsedGo2Genes)

    def GOTermEnrichment(self,GoDb):
        self.GoDb=GoDb
        self.Go2Genes={}
        self.allGenes=set()
        return self

    def setSelectGenes(self,selectGenesFile):
        curSelectGenes=set()
        with open(selectGenesFile,'r') as fr:
            for genes in fr.readlines():
                genes=genes.strip()
                curSelectGenes.add(genes)
        self.selectGenes=curSelectGenes
        return self

    def Gene2GoToGo2Gene(self,Gene2Go,Go2Gene):
        for GeneId in Gene2Go.keys():
            for GoId in Gene2Go[GeneId]:
                if GoId in Go2Gene:
                    Go2Gene[GoId].add(GeneId)
                else:
                    newGenes=set()
                    newGenes.add(GeneId)
                    Go2Gene[GoId]=newGenes

    def Go2GeneToGene2Go(self,Go2Gene,Gene2Go):
        for GoId in Go2Gene.keys():
            for geneId in Go2Gene[GoId]:
                if geneId in Gene2Go:
                    Gene2Go[geneId].add(GoId)
                else:
                    newGos=set()
                    newGos.add(GoId)
                    Gene2Go[geneId]=newGos

    def readGene2GoFile(self,fileForEnrich):
        logging.warning("read file:\t"+fileForEnrich)
        with open(fileForEnrich,'r') as fr:
            for line in fr.readlines():
                line=line.strip()
                inlineInfo=line.split("\t")
                if len(inlineInfo) >= 2:
                    geneIds=inlineInfo[0]
                    goId=inlineInfo[1]
                    if geneIds not in self.Gene2GO:
                        nexGos=set()
                        nexGos.add(goId)
                        self.Gene2GO[geneIds]=nexGos
                    else:
                        self.Gene2GO[geneIds].add(goId)
        self.Gene2GoToGo2Gene(self.Gene2GO,self.Go2Genes)
        logging.warning("GeneNums:"+str(len(self.Gene2GO)))
        logging.warning("GoNums:"+str(len(self.Go2Genes)))
        return self

    def parseAllGO(self):
        if not self.Gene2GO:
            pass
        for geneId in self.Gene2GO.keys():
            GoList=[]
            if len(self.Gene2GO[geneId]) == 0:
                logging.warning("GeneId: "+geneId)
            GoList.extend(self.Gene2GO[geneId])
            if geneId in self.parsedGene2Go:
                logging.warning("This is a big bug")
            else:
                Gos=set()
                for gi in self.GoDb.getAncestor_list(GoList):
                    Gos.add(gi)
                self.parsedGene2Go[geneId]=Gos

        self.Gene2GoToGo2Gene(self.parsedGene2Go,self.parsedGo2Genes)
        logging.warning("Parsed GeneNums:"+str(len(self.parsedGene2Go)))
        logging.warning("Parsed GoNum:"+str(len(self.parsedGo2Genes)))
        return self

    def writeOutParsedGOFile(self,outFilePath):
        if self.parsedGo2Genes:
            pass
        with open(outFilePath,'w') as fw:
            for GoId in self.parsedGo2Genes.keys():
                if GoId in self.GoDb.getGoHashDb():
                    sb=[]
                    sb.append(GoId)
                    sb.append("\t")
                    sb.append(self.GoDb.getGoHashDb()[GoId].getNameSpace())
                    sb.append("\t")
                    sb.append(self.GoDb.getGoHashDb()[GoId].getName())
                    sb.append("\t")
                    for geneId in self.parsedGo2Genes[GoId]:
                        sb.append(geneId)
                        sb.append(",")
                    sb.append("\n")
                    fw.write("".join(sb))
                else:
                    logging.warning("If u are not enriching for goslim, Skip GOID["+GoId+"] Because it is not included in given obo file and u may need to update the go-basic.obo file\n")

    def splitTo3GoPart(self,Go2Genes,nameSpaceShortName):
        if nameSpaceShortName == "CC":
            nameSpace="cellular_component"
        elif nameSpaceShortName == "BP":
            nameSpace="biological_process"
        elif nameSpaceShortName == "MF":
            nameSpace="molecular_function"
        elif nameSpaceShortName == "cellular_component":
            nameSpace="cellular_component"
        elif nameSpaceShortName == "biological_process":
            nameSpace="biological_process"
        elif nameSpaceShortName == "molecular_function":
            nameSpace="molecular_function"
        else:
            logging.warning("Error! unknown NameSpace...")
            return GOTermEnrichment()
        returnGOTermEnrichMent=GOTermEnrichment().GOTermEnrichment(self.GoDb)
        GoHashDb=self.GoDb.getGoHashDb()
        curGo2Genes={}
        for GoId in Go2Genes.keys():
            if not GoId in GoHashDb:
                logging.warning("If u are not enriching for goslim, Skip GOID[" + GoId + "] Because it is not included in given obo file and u may need to update the go-basic.obo file\n")
            elif GoHashDb[GoId].getNameSpace() == nameSpace:
                curGo2Genes[GoId]=Go2Genes[GoId]
        returnGOTermEnrichMent.setGo2Genes(curGo2Genes)
        curGenes2Go={}
        self.Go2GeneToGene2Go(curGo2Genes,curGenes2Go)
        returnGOTermEnrichMent.setGene2Go(curGenes2Go)
        logging.warning(nameSpace+":"+str(len(curGenes2Go.keys())))
        return returnGOTermEnrichMent

    def prepareForEnrichMent(self,goBasicOboFile,Gene2GoFile):
        logging.warning("Parsing GO obo database File....")
        GOdb=GObasicOboDatabase()
        GOdb.setDatabaseFile(goBasicOboFile)
        logging.warning("Making GO Hash database....")
        goHashDb=GOdb.makeDb().getDbHash()
        glg=GOLevelGetter(goHashDb)
        self.setGoDb(glg)
        logging.warning("Parsing BackGround File....")
        self.readGene2GoFile(Gene2GoFile)
        logging.warning("Parsing All GO term....")
        self.parseAllGO()
        parsedAllGoFilePath=Gene2GoFile+"_ParsedAllGO.xls"
        logging.warning("Write All Parsed GO Map as "+parsedAllGoFilePath)
        self.writeOutParsedGOFile(parsedAllGoFilePath)
        logging.warning("Preparation Finished....")

    def AutoEnrichMent(self,selectionGenesFile,outDir):
        logging.warning("Begin Auto EnrichMents....")
        logging.warning("Reading Selection Set DataSet....")
        self.setSelectGenes(selectionGenesFile)
        ThreeTermPart=["MF","CC","BP"]
        for curTP in ThreeTermPart:
            logging.warning("Enrishing in "+curTP+"......")
            mixResult = self.doEnrichMent(self.getSelectGenes(),curTP)
            if outDir:
                TPOutFile=outDir+"/"+os.path.basename(selectionGenesFile)+"_"+curTP+"_EnrichResult.xls"
            else:
                TPOutFile=selectionGenesFile+"_"+curTP+"_EnrichResult.xls"
            with open(TPOutFile,'w') as fw:
                fw.write("GO_Name\tGO_ID\tGO_Level\tP_value\tEnrichmentScore\tHitsGenesCountsInSelectedSet\tHitsGenesCountsInBackground\tAllGenesCountsInSelectedSet\tAllGenesCountsInBackground\tGenesOfSelectedSetInGOterm\n")
                goHashDb = self.GoDb.getGoHashDb()
                for GoId in mixResult.keys():
                    gse = mixResult[GoId]
                    sb=[]
                    for geneId in self.getSelectGenes():
                        if geneId in self.parsedGene2Go:
                            if GoId in self.parsedGene2Go[geneId]:
                                sb.append(geneId)
                                sb.append(",")
                    if len(sb) >0:
                        fw.write(goHashDb[GoId].getName()+"\t")
                        fw.write(GoId+"\t")
                        fw.write(str(self.GoDb.getLevel(GoId))+"\t")
                        fw.write(str(gse.getPValue())+"\t")
                        fw.write(str(gse.getEnrichment())+"\t")
                        fw.write(str(gse.getAllSelectWhiteBalls()) + "\t")
                        fw.write(str(gse.getAllWhiteBalls()) + "\t")
                        fw.write(str(gse.getAllSelectBalls()) + "\t")
                        fw.write(str(gse.getAllBalls()) + "\t")
                        fw.write(str("".join(sb)) + "\n")
            data = pd.read_csv(TPOutFile,header=0,sep="\t",index_col=0,encoding="utf-8")
            # logging.warning(data.columns)
            # data.columns = ["GO_ID","GO_Level","P_value","EnrichmetSocre","HitsGenesCountsInSelectedSet", "HitsGenesCountsInBackground","AllGenesCountsInSelectedSet", "AllGenesCountsInBackground", "GenesOfSelectedSetInGOterm"]
            datasort = data.sort_values(by="P_value", axis=0)
            pValue = datasort["P_value"]
            try:
                fdr = self.p_adjust_bonferroni(pValue)
                fdr = pd.DataFrame(fdr, index=datasort.index)
                datasort.insert(9, "q-value(bonferroni method)", fdr)
                datasort.to_csv(TPOutFile.replace(".xls","") + '.sorted.padjust.xls', encoding='utf-8', index=True, sep="\t",
                                index_label="GO_Name")
            except ZeroDivisionError:
                logging.warning("This GO Enrichment Result: "+TPOutFile+"_"+curTP+"_EnrichResult.xls is not complete,you can't get qvalue file!")
                continue


    def AutoEnrich(self,selectionGenesFile):
        self.AutoEnrichMent(selectionGenesFile,"")

    def p_adjust_bh(self,p):
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

    def p_adjust_bonferroni(self,p):
        p = np.asfarray(p)
        by_descend = p.argsort()[::-1]
        # by_orig = by_descend.argsort()
        qvalue = multipletests(p, method='bonferroni')
        # return qvalue[1][by_orig], qvalue[0][by_orig]
        return qvalue[1]

    def doEnrichMent(self,selectedGenes,mode):
        curEnrich={}
        curGOTermEnrichMent=GOTermEnrichment().GOTermEnrichment(self.GoDb)
        curGene2Go={}
        for geneId in selectedGenes:
            if geneId in self.parsedGene2Go:
                curGoSet=set(self.parsedGene2Go[geneId])
                curGene2Go[geneId]=curGoSet
            else:
                logging.warning("Skip Selected Gene With Out GO Annotation:" + geneId)
        curGOTermEnrichMent.setGene2Go(curGene2Go)
        curGOTermEnrichMent.parseAllGO()
        if mode == "Mix":
            allBalls = len(self.parsedGene2Go)
            allSelectBalls = len(curGOTermEnrichMent.parsedGene2Go)
            for GoId in curGOTermEnrichMent.parsedGo2Genes:
                allWhiteBalls = len(self.parsedGo2Genes[GoId])
                allSelectWhiteBalls = len(curGOTermEnrichMent.parsedGo2Genes[GoId])
                gse=Enricher()
                gse.setAllBalls(int(allBalls)).setAllSelectBalls(int(allSelectBalls)).setAllWhiteBalls(int (allWhiteBalls)).setAllSelectWhiteBalls(int(allSelectWhiteBalls)).doHyperGeometric()
                curEnrich[GoId]=gse
        else:
            allWithMode=self.splitTo3GoPart(self.parsedGo2Genes,mode)
            selectedWithMode=curGOTermEnrichMent.splitTo3GoPart(curGOTermEnrichMent.parsedGo2Genes,mode)
            logging.warning("Parsed All GO in mode")
            allWithMode.parseAllGO()
            logging.warning("Parsed All Selected GO in mode")
            selectedWithMode.parseAllGO()
            allBalls=len(allWithMode.parsedGene2Go)
            allSelectBalls=len(selectedWithMode.parsedGene2Go)
            for GoId in selectedWithMode.parsedGo2Genes:
                allWhiteBalls=len(allWithMode.parsedGo2Genes[GoId])
                allSelectWhiteBalls=len(selectedWithMode.parsedGo2Genes[GoId])
                gse=Enricher()
                gse.setAllBalls(int(allBalls)).setAllSelectBalls(int(allSelectBalls)).setAllWhiteBalls(
                    int(allWhiteBalls)).setAllSelectWhiteBalls(int(allSelectWhiteBalls)).doHyperGeometric()
                curEnrich[GoId]=gse
        return curEnrich

    def getOriLevelAt(self,num):
        oriGo2GeneOfLevel={}
        for GoNum in self.Go2Genes.keys():
            if self.GoDb.getLevel[GoNum] == num:
                oriGo2GeneOfLevel[GoNum]=self.Go2Genes[GoNum]
        return oriGo2GeneOfLevel

    def getParsedLevelAt(self,num):
        parsedGo2GeneOfLevel={}
        for GoNum in self.parsedGo2Genes.keys():
            if self.GoDb.getLevel(GoNum) == num:
                parsedGo2GeneOfLevel[GoNum]=self.parsedGo2Genes[GoNum]
        return parsedGo2GeneOfLevel

    def test(self):
        logging.warning("start go enrich......")
        GOdb = GObasicOboDatabase()
        GOdb.setDatabaseFile("gene_ontology.1_2.obo")
        goHashDb = GOdb.makeDb().getDbHash()
        glg = GOLevelGetter(goHashDb)
        curEnrich = GOTermEnrichment().GOTermEnrichment(glg)
        curEnrich.readGene2GoFile("hg19.c.annot")
        curEnrich.parseAllGO()
        curEnrich.writeOutParsedGOFile("test.xls")
        curEnrich.setSelectGenes("fg.lst")
        mixResult = curEnrich.doEnrichMent(curEnrich.getSelectGenes(), "Mix")
        with open("test2.xls", 'w') as fw:
            for GoId in mixResult.keys():
                gse = mixResult[GoId]
                logging.warning(GoId + "\t" + str(gse.getPValue()) + "\t" + str(gse.getEnrichment()) + "\t" + str(
                    gse.getAllSelectWhiteBalls()) + "\t" + str(gse.getAllWhiteBalls()) + "\t" + str(
                    gse.getAllSelectBalls()) + "\t" + str(gse.getAllBalls()))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="do GOenrichment Slim from gene ontology in obo format")
    parser.add_argument('--oboFile', dest='oboFile', type=str,
                        help="gene ontology in obo format,you could download it from http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo")
    parser.add_argument('--selectionSetFiles', dest='selectFile', type=str,
                        help="DEG List you want to do GOenrichment, you just need to set the Gene name in the list.")
    parser.add_argument('--gene2GoFile',dest='backFile',type=str,
                        help="the Gene GO annotation of a species, format is like this \nGene_Id\tGO_Id\nGene01\tGO:000001\n")
    parser.add_argument('--outdir',dest='outdir',default="",help="set outdir for the GOEnrichment result.")
    localeArg = parser.parse_args()

    if localeArg.oboFile is not None and localeArg.selectFile is not None and localeArg.backFile is not None and localeArg.outdir == "" :
        logging.warning(time.asctime() + "\tStart GoTermEnrichment!")
        file=localeArg.oboFile
        fg=localeArg.selectFile
        bg=localeArg.backFile
        enrich=GOTermEnrichment()
        enrich.prepareForEnrichMent(file,bg)
        enrich.AutoEnrich(fg)
        logging.warning(time.asctime() +"\tGoTermEnrichment Finish!")
    elif localeArg.oboFile is not None and localeArg.selectFile is not None and localeArg.backFile is not None and localeArg.outdir != "":
        logging.warning(time.asctime() + "\tStart GoTermEnrichment!")
        file=localeArg.oboFile
        fg=localeArg.selectFile
        bg=localeArg.backFile
        od=localeArg.outdir
        enrich=GOTermEnrichment()
        enrich.prepareForEnrichMent(file,bg)
        enrich.AutoEnrichMent(fg,od)
        logging.warning(time.asctime() +"\tGoTermEnrichment Finish!")
    else:
        print('Use python --help')
        sys.exit()