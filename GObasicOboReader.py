#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
from  GOoboTerm import GOoboTerm
from SingleLabelReaderAdvance import SingleLabelReaderAdvance
import re
import logging

""" 
@author:yueyao 
@file: GObasicOboReader.py 
@time: 2018/06/08 
"""

class GObasicOboReader(SingleLabelReaderAdvance):

    def __init__(self):
        super().__init__()
        self.setEndPattern("^\[Typedef\]")
        self.setIsSkipHeader(True)
        self.setLabelPatten("^\[Term\]")
        self.info=[]
        self.id_regex01 = re.compile('^id:\s+(GO:.*)')
        self.name_regex02 = re.compile('^name:\s+(.*)')
        self.namespace_regex03 = re.compile('^namespace:\s+(.*)')
        self.altid_regex05 = re.compile('^alt_id:\s+(GO:.*)')
        self.isa_regex06 = re.compile('^is_a:\s+(GO:[^\s]+)')
        self.relation_regex07 = re.compile('^relationship:.*(GO:[\d]+)')
        self.part_of=re.compile("part_of")
        self.regulates=re.compile("regulates")
        self.positively_regulates=re.compile("positively_resulates")
        self.negatively_regulates=re.compile("negatively_regulates")

    def makeGOoboTermFromFile(self,info):
        goOboTerm = GOoboTerm()
        for line in info:
            line=line.strip()
            if line.startswith("[Typedef]"):
                break
            if line.startswith("subsetdef") or line.startswith("synonymtypedef") or line.startswith("default"):
                continue
            if not line.startswith("[Term]"):
                if line.startswith("id"):
                    goOboTerm.setId(line.split(": ")[1])
                elif line.startswith("name:"):
                    goOboTerm.setName(line.split(": ")[1])
                elif line.startswith("namespace"):
                    goOboTerm.setNameSpace(line.split(": ")[1])
                elif line.startswith("def:"):
                    goOboTerm.setDef(line.split(": ")[1])
                elif line.startswith("comment"):
                    goOboTerm.setComment(line.split(": ")[1])
                elif line.startswith("is_obsolete"):
                    goOboTerm.setIs_obsolete(True)
                elif line.startswith("consider"):
                    goOboTerm.getConsiders().append(line.split(": ")[1])
                elif line.startswith("replaced_by"):
                    goOboTerm.getRepalced_by().append(line.split(": ")[1])
                elif line.startswith("is_a"):
                    goOboTerm.getIs_a().append(re.findall(self.isa_regex06, line)[0])
                elif line.startswith("alt_id"):
                    goOboTerm.getAlt_ids().append(line.split(": ")[1])
                elif line.startswith("subset:"):
                    goOboTerm.getSubsets().append(line.split(": ")[1])
                elif line.startswith("synonym:"):
                    goOboTerm.getSynonym().append(line.split(": ")[1])
                elif line.startswith("xref"):
                    goOboTerm.getXref().append(line.split(": ")[1])
                elif line.startswith("relationship"):
                    relationshipInfo=line.split(' ')
                    if re.search(self.part_of,relationshipInfo[1]):
                        goOboTerm.getPart_of().append(relationshipInfo[2])
                    elif re.search(self.regulates,relationshipInfo[1]):
                        goOboTerm.getRegulates().append(relationshipInfo[2])
                    elif re.search(self.positively_regulates,relationshipInfo[1]):
                        goOboTerm.getPositively_regulates().append(relationshipInfo[2])
                    elif re.search(self.negatively_regulates,relationshipInfo[1]):
                        goOboTerm.getNegatively_regulates().append(relationshipInfo[2])
                    else:
                        logging.warning("[Unknow Relationship]: "+line)
        return goOboTerm

    def getNextTerm(self):
        curRecord=self.getNext()
        return self.makeGOoboTermFromFile(curRecord)

if __name__ == '__main__':
    info="gene_ontology.1_2.obo"
    obo=GObasicOboReader()
    GOterm =obo.makeGOoboTermFromFile()
    print(GOterm.toString())