#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GObasicOboDatabase.py
@time: 2019/01/18
'''

from GObasicOboReader import GObasicOboReader

class GObasicOboDatabase(object):

    def __init__(self):
        self.databaseFile=""
        self.DbHash={}

    def getDatabaseFile(self):
        return self.databaseFile

    def setDatabaseFile(self,databaseFile):
        self.databaseFile = databaseFile
        return self

    def getDbHash(self):
        return self.DbHash

    def setDbHash(self,DbHash):
        self.DbHash=DbHash

    def GObasicOboDatabase(self):
        self.DbHash={}

    def makeDb(self):
        gor = GObasicOboReader()
        gor.setTargetFile(self.databaseFile)
        while gor.hasNext():
            got = gor.getNextTerm()
            self.DbHash[got.getId()]=got
        return self

if __name__ == '__main__':
    GOdb=GObasicOboDatabase()
    GOdb.setDatabaseFile("gene_ontology.1_2.obo")
    goHashDb=GOdb.makeDb().getDbHash()
    for iv,ik in goHashDb.items():
        print(iv)