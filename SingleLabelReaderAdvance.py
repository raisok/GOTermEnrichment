#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import re
import logging

'''
@author: yueyao
@file: SingleLabelReaderAdvance.py
@time: 2019/01/17
'''

class SingleLabelReaderAdvance(object):

    def __init__(self):
        self.tagetFile=""
        self.targetFilePath=""
        self.headerArr=[]
        self.returnArrayList=[]
        self.isSkipHeader=None
        self.inline=None
        self.endPattern=""
        self.endPattenString=""
        self.labelPattenString=""
        self.labelPatten=""
        self.br=None
        self.flag=True

    def setHeaderArr(self,headerArr):
        self.headerArr=headerArr

    def getHeaderArr(self):
        return self.headerArr

    def setIsSkipHeader(self,isSkipHeader):
        self.isSkipHeader=isSkipHeader

    def isIsSkipHeader(self):
        return self.isSkipHeader

    def setEndPattern(self,endPattern):
        self.endPattern=endPattern
        return self

    def getEndPattern(self):
        return self.endPattern

    def getEndPattenString(self):
        return self.endPattenString

    def setEndPattenString(self,endPattenString):
        self.endPattenString=endPattenString
        self.endPattern=re.compile(self.endPattenString)
        return self

    def getTagetFile(self):
        return self.tagetFile

    def getLabelPatten(self):
        return self.labelPatten

    def setLabelPatten(self,labelPatten):
        self.labelPatten=labelPatten
        return self

    def getLabelPattenString(self):
        return self.labelPattenString

    def setLabelPattenString(self,labelPattenString):
        self.labelPattenString=labelPattenString
        self.labelPatten=re.compile(labelPattenString)
        return self

    def hasNext(self):
        return self.flag

    def getNext(self):
        self.returnArrayList.clear()
        self.returnArrayList.append(self.inline)
        for self.inline in self.br:
            self.inline=self.inline.strip()
            if self.endPattern != None:
                if re.match(self.endPattern,self.inline):
                    logging.warning("END: meet Pattern :"+self.endPattern)
                    self.flag = False
                    self.br.close()
                    return self.returnArrayList
            if re.findall(self.labelPatten,self.inline):
                return self.returnArrayList
            self.returnArrayList.append(self.inline)
        return self.returnArrayList

    def setTargetFile(self,targetFile):
        F=open(targetFile,'r')
        if self.isSkipHeader:
            for self.inline in F:
                self.inline=self.inline.strip()
                if re.findall(self.labelPatten,self.inline):
                    logging.warning("Skip: until meet Pattern "+str(self.labelPatten))
                    break
                self.headerArr.append(self.inline)
        self.inline=F.readline().strip()
        self.br=F
        return self

if __name__ == '__main__':
    slr = SingleLabelReaderAdvance()
    slr.setIsSkipHeader(True)
    slr.setLabelPatten("^\[Term\]")
    slr.setEndPattern("^\[Typedef\]")
    slr.setTargetFile("gene_ontology.1_2.obo")
    while slr.hasNext():
        curArr = slr.getNext()
        cursize=len(curArr)
        logging.warning("curArr size\t"+str(cursize))
