#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GOTermEnrichmentGOslim.py
@time: 2019/01/22
'''

import logging
import time
import argparse
import sys
from GOTermEnrichment import GOTermEnrichment


class GOTermEnrichmentGOslim(GOTermEnrichment):
    def __init__(self):
        super().__init__()
        self.goSlimSelected=""

    def getGoSlimSelected(self):
        return self.goSlimSelected

    def setGoSlimSelected(self,goSlimSelected):
        self.goSlimSelected=goSlimSelected

    def parseAllGO(self):
        if not self.Gene2GO:
            pass
        for geneId in self.Gene2GO.keys():
            GoList=[]
            if len(self.Gene2GO[geneId]) == 0:
                logging.warning("GeneId: "+geneId)
            GoList.extend(self.Gene2GO[geneId])
            if geneId in self.parsedGene2Go:
                logging.warning("如果运行到这里，应该是出了bug!")
            else:
                Gos=set()
                for gi in self.GoDb.getAncestor_list(GoList):
                    Gos.add(gi)
                slimGos=set()
                for CurGo in Gos:
                    if self.goSlimSelected in self.GoDb.getGoHashDb()[CurGo].getSubset():
                        slimGos.add(CurGo)
                if len(slimGos) >0:
                    self.parsedGene2Go[geneId]=slimGos
                else:
                    logging.warning("Gene Without Slim Go annotation:"+geneId)
        self.Gene2GoToGo2Gene(self.parsedGene2Go,self.parsedGo2Genes)
        logging.warning("Parsed GeneNums:"+str(len(self.parsedGene2Go)))
        logging.warning("Parsed GoNum:"+str(len(self.parsedGo2Genes)))
        return self

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="do GOenrichment Slim from gene ontology in obo format")
    parser.add_argument('--oboFile', dest='oboFile', type=str,
                        help="gene ontology in obo format,you could download it from http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo")
    parser.add_argument('--selectionSetFiles', dest='selectFile', type=str,
                        help="DEG List you want to do GOenrichment, you just need to set the Gene name in the list.")
    parser.add_argument('--gene2GoFile',dest='backFile',type=str,
                        help="the Gene GO annotation of a species, format is like this \nGene_Id\tGO_Id\nGene01\tGO:000001\n")
    parser.add_argument('--goslim',dest='goslim',default="",help="set slim set to do go term enrichment, Options are as followed:\n\t\t'goslim_aspergillus'\twhich means Aspergillus GO slim\n\t\t'goslim_candida'\twhich means Candida GO slim\n\t\t'goslim_chembl'\twhich means ChEMBL protein targets summary\n\t\t'goslim_generic'\twhich means Generic GO slim\n\t\t'goslim_goa\twhich' means GOA and proteome slim\n\t\t'goslim_metagenomics'\twhich means Metagenomics GO slim\n\t\t'goslim_pir'\twhich' means PIR GO slim\n\t\t'goslim_plant'\twhich means Plant GO slim\n\t\t'goslim_pombe'\twhich means Fission yeast GO slim\n\t\t'goslim_synapse'\twhich means synapse GO slim\n\t\t'goslim_virus'\twhich means Viral GO slim\n\t\t'goslim_yeast'\twhich means Yeast GO slim")
    parser.add_argument('--outdir',dest='outdir',help="set outdir for the GOEnrichment result.")
    localeArg = parser.parse_args()

    if localeArg.oboFile is not None and localeArg.goslim is not None and localeArg.selectFile is not None and localeArg.backFile is not None and localeArg.outdir is None:
        logging.warning(time.asctime() + "\tStart GoTermEnrichment Slim!")
        file=localeArg.oboFile
        selectFile=localeArg.selectFile
        backfile=localeArg.backFile
        goslim=localeArg.goslim
        curGOEnrich = GOTermEnrichmentGOslim()
        curGOEnrich.setGoSlimSelected(goslim)
        curGOEnrich.prepareForEnrichMent(file, backfile)
        curGOEnrich.AutoEnrich(selectFile)
        logging.warning(time.asctime() +"\tGoTermEnrichment Slim Finish!")
    elif localeArg.oboFile is not None and localeArg.selectFile is not None and localeArg.backFile is not None and localeArg.outdir is not None:
        logging.warning(time.asctime() + "\tStart GoTermEnrichment Slim!")
        file=localeArg.oboFile
        selectFile=localeArg.selectFile
        backfile=localeArg.backFile
        goslim=localeArg.goslim
        od=localeArg.outdir
        curGOEnrich = GOTermEnrichmentGOslim()
        curGOEnrich.setGoSlimSelected(goslim)
        curGOEnrich.prepareForEnrichMent(file, backfile)
        curGOEnrich.AutoEnrichMent(selectFile,od)
        logging.warning(time.asctime() +"\tGoTermEnrichment Slim Finish!")
    else:
        print('Use python --help')
        sys.exit()

    logging.warning(time.asctime()+"")
    goslim=localeArg.goslim
    obofile=localeArg.oboFile
    backfile=localeArg.backFile
    selectFile=localeArg.selectFile
    curGOEnrich=GOTermEnrichmentGOslim()
    curGOEnrich.setGoSlimSelected(goslim)
    curGOEnrich.prepareForEnrichMent(localeArg,backfile)
    curGOEnrich.AutoEnrichMent(selectFile)

    logging.warning(time.asctime()+"")