#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GOoboTerm.py
@time: 2019/01/18
'''

class GOoboTerm(object):

    def __init__(self):
        self.id=None
        self.nameSpace=None
        self.name=None
        self.def_=None
        self.comment=None
        self.is_obsolete=None
        self.considers=[]
        self.replaced_by=[]
        self.is_a=[]
        self.alt_ids=[]
        self.subsets=[]
        self.synonym=[]
        self.xref=[]
        self.negatively_regulates=[]
        self.part_of=[]
        self.positively_regulates=[]
        self.regulates=[]

    def setId(self,id):
        self.id=id

    def getId(self):
        return self.id

    def setName(self,name):
        self.name=name

    def getName(self):
        return self.name

    def setNameSpace(self,nameSpace):
        self.nameSpace=nameSpace

    def getNameSpace(self):
        return self.nameSpace

    def setDef(self,def_):
        self.def_=def_

    def getDef(self):
        return self.def_

    def setComment(self,comment):
        self.comment=comment

    def getComment(self):
        return self.comment

    def setIs_obsolete(self,is_obsolete):
        self.is_obsolete=is_obsolete

    def isIs_obsolete(self):
        return self.is_obsolete

    def getConsiders(self):
        return self.considers

    def setConsiders(self,considers):
        self.considers=considers

    def setRepalced_by(self,replaced_by):
        self.replaced_by=replaced_by

    def getRepalced_by(self):
        return self.replaced_by

    def setIs_a(self,is_a):
        self.is_a=is_a

    def getIs_a(self):
        return self.is_a

    def setAlt_ids(self,alt_ids):
        self.alt_ids=alt_ids

    def getAlt_ids(self):
        return self.alt_ids

    def setSubsets(self,subsets):
        self.subsets=subsets

    def getSubsets(self):
        return self.subsets

    def setSynonym(self,synonym):
        self.synonym=synonym

    def getSynonym(self):
        return self.synonym

    def setXref(self,xref):
        self.xref=xref

    def getXref(self):
        return self.xref

    def setNegatively_regulates(self,negatively_regulates):
        self.negatively_regulates=negatively_regulates

    def getNegatively_regulates(self):
        return self.negatively_regulates

    def setPart_of(self,part_of):
        self.part_of=part_of

    def getPart_of(self):
        return self.part_of

    def setPositively_regulates(self,positively_regulates):
        self.positively_regulates=positively_regulates

    def getPositively_regulates(self):
        return self.positively_regulates

    def setRegulates(self,regulates):
        self.regulates=regulates

    def getRegulates(self):
        return self.regulates

    def toString(self):
        return "GoTerm{id=" + str(self.id) + ", name=" + str(self.name) + ", nameSpace=" + str(self.nameSpace) + ", def=" + str(self.def_) + ", comment=" + str(self.comment) + ", is_obsolete=" + str(self.is_obsolete) + ", considers=" + str(self.considers) + ", repalced_by=" + str(self.replaced_by) + ", is_a=" + str(self.is_a) + ", alt_ids=" + str(self.alt_ids) + ", subsets=" + str(self.subsets) + ", synonym=" + str(self.synonym) + ", xref=" + str(self.xref) + ", negatively_regulates=" + str(self.negatively_regulates) + ", part_of=" + str(self.part_of) + ", positively_regulates=" + str(self.positively_regulates) + ", regulates=" + str(self.regulates) + "}"

if __name__=='__main__':
    obo=GOoboTerm()
    print(obo.toString())

