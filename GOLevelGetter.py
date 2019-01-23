#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GOLevelGetter.py
@time: 2019/01/17
'''

import logging
import os
from GObasicOboDatabase import GObasicOboDatabase


class GOLevelGetter(object):

    def __init__(self,GoHashDb):
        self.GoHashDb=GoHashDb

    def setGoHashDb(self,GoHashDb):
        self.GoHashDb=GoHashDb

    def getGoHashDb(self):
        return self.GoHashDb

    def getLevelList(self,GoNum):
        curGo=self.GoHashDb[GoNum]
        allDis=[]
        self.CountDistance(curGo,allDis,2)
        return allDis

    def getLevel(self,GoNum):
        Max = 0
        LevelList=self.getLevelList(GoNum)
        for curLevel in LevelList:
            if Max < curLevel:
                Max=curLevel
        return Max


    def getGOinLevel(self,LevelNum):
        goLevelList=[]
        for CurGO in self.GoHashDb.keys():
            if self.getLevel(CurGO) == LevelNum:
                goLevelList.append(CurGO)
        return goLevelList



    def CountDistance(self,inGo,allDis,Counts):
        inGoParents = self.getDirectIs_aParent(inGo.getId())
        # logging.warning("Check\t"+str(inGo.getId())+"\t"+str(len(inGoParents)))
        for curGo in inGoParents:
            curCounts=Counts
            if curGo is None:
                os._exit(1)
            if not curGo.isIs_obsolete():
                if len(self.getDirectIs_aParent(curGo.getId())) == 0:
                    allDis.append(curCounts)
                else:
                    curCounts+=1
                    self.CountDistance(curGo,allDis,curCounts)

    def getDirectIs_aParent(self,GoNum):
        GOparent=[]
        curGo = self.GoHashDb[GoNum]
        if len(curGo.getIs_a()) != 0:
            for addGO in curGo.getIs_a():
                GOparent.append(self.GoHashDb[addGO])
        return GOparent

    def getDirectChildren(self,GoNum):
        GOchildren=[]
        for curTerm in self.GoHashDb.values():
            if GoNum in curTerm.getPart_of():
                GOchildren.append(curTerm)
            elif GoNum in curTerm.getIs_a():
                GOchildren.append(curTerm)
        return GOchildren


    def getDirectParent(self,GoNum):
        curGo = self.GoHashDb[GoNum]
        GOparents=[]
        if len(curGo.getRegulates()) !=0:
            for addGO in curGo.getRegulates():
                GOparents.append(self.GoHashDb[addGO])
        if len(curGo.getNegatively_regulates()) !=0:
            for addGO in curGo.getNegatively_regulates():
                GOparents.append(self.GoHashDb[addGO])
        if len(curGo.getPositively_regulates()) !=0:
            for addGO in curGo.getPositively_regulates():
                GOparents.append(self.GoHashDb[addGO])
        if len(curGo.getPart_of()) !=0:
            for addGO in curGo.getPart_of():
                GOparents.append(self.GoHashDb[addGO])
        if len(curGo.getIs_a()) !=0:
            for addGO in curGo.getIs_a():
                GOparents.append(self.GoHashDb[addGO])
        return GOparents

    def parseAncestors(self,inGo,allAncestor):
        inGoParents=self.getDirectParent(inGo.getId())
        for curGO in inGoParents:
            if not curGO.isIs_obsolete():
                if len(self.getDirectParent(curGO.getId())) != 0:
                    allAncestor.append(curGO.getId())
                    self.parseAncestors(curGO,allAncestor)

    def getAncestorRedundance(self,GoNum):
        ancestors=[]
        if GoNum in self.GoHashDb.keys():
            self.parseAncestors(self.GoHashDb[GoNum],ancestors)
        else:
            logging.warning("[Error]: Sorry GO ID "+GoNum+" is not exist in this DB...\nRemember we just need GO id in Integer Format like 988 represents GO:0000988")
        return ancestors

    def getAncestor(self,GoNum):
        uniqAncestors=[]
        uniqAncestors.append(GoNum)
        if GoNum in self.GoHashDb.keys():
            curGO = self.GoHashDb[GoNum]
            if not curGO.isIs_obsolete():
                if curGO.getNameSpace() == "biological_process":
                    uniqAncestors.append("GO:0008150")
                elif curGO.getNameSpace() == "cellular_component":
                    uniqAncestors.append("GO:0005575")
                elif curGO.getNameSpace() == "molecular_function":
                    uniqAncestors.append("GO:0003674")
            ancestors=self.getAncestorRedundance(GoNum)
            uniq=set()
            for i in ancestors:
                if i not in uniq:
                    uniq.add(i)
                    uniqAncestors.append(i)
        else:
            for curGONum in self.GoHashDb.keys():
                if GoNum in self.GoHashDb[curGONum].getAlt_ids():
                    return self.getAncestor(curGONum)
            logging.warning("[Hard Warning]:"+str(GoNum)+" is not found in go-basis.obo file, u may need to update go-basis.obo file")
        return uniqAncestors


    def getAncestor_list(self,GoNumList):
        allAncestors=[]
        for GoNum in GoNumList:
            allAncestors.extend(self.getAncestor(GoNum))
        return allAncestors


    def getGoListAtLevel(self):
        GoListAtLevel={}
        for go in self.GoHashDb.values():
            if not go.isIs_obsolete():
                curLevel=self.getLevel(go.getId())
                if curLevel in GoListAtLevel:
                    GoListAtLevel[curLevel].append(go.getId())
                else:
                    newGoList=[]
                    newGoList.append(go.getId())
                    GoListAtLevel[curLevel]=newGoList
        return GoListAtLevel

    def writeGOtermWithLevel(self,outFile):
        with open(outFile,'w') as fw:
            for go in self.GoHashDb.values():
                if not go.isIs_obsolete():
                    fw.write(go.getId()+"\t"+go.getName()+"\t"+go.getNameSpace()+"\t"+str(self.getLevel(go.getId()))+"\n")

if __name__ == '__main__':
    GOdb=GObasicOboDatabase()
    GOdb.setDatabaseFile("./data/gene_ontology.1_2.obo")
    goHashDb=GOdb.makeDb().getDbHash()
    glg = GOLevelGetter(goHashDb)
    logging.warning(glg.getLevel("GO:0044422"))
    logging.warning(glg.getAncestor("GO:0044422"))
    # logging.warning("+++++++")
    logging.warning(glg.getLevel("GO:0005581"))