#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: dealGOObo.py
@time: 2019/01/11
'''

import argparse
import re
import sys
import logging
import time

def delimited(file, delimiter = '\n', bufsize = 4096):
    buf = ''
    while True:
        newbuf = file.read(bufsize)
        if not newbuf:
            yield buf
            return
        buf += newbuf
        lines = buf.split(delimiter)
        for line in lines[:-1]:
            yield line
        buf = lines[-1]

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

def write_alias(alias,obsoletes,prefix):
    logging.warning(time.asctime()+"\twrite file:"+prefix+".alias")
    with open (prefix+'.alias','w') as fw:
        for id,mm in alias.items():
            if id in obsoletes:
                continue
            #id alias id
            fw.write(id+'\t'+"\t".join(mm)+"\n")

def write_class(gos,obsoletes,parents,prefix):
    logging.warning(time.asctime()+"\twrite file:" + prefix + ".class")
    with open (prefix+'.class','w') as fw:
        for id,im in gos.items():
            if id in obsoletes:
                continue
            temps={}
            getParents(id,parents,temps)
            for iv,ik in temps.items():
                #first class second class term id term
                #namespase  parentsidname  goid    goname
                fw.write("\t".join([gos[id][1],gos[iv][0],id,gos[id][0]])+"\n")

def get_parent(parents,prefix):
    logging.warning(time.asctime()+"\twrite file:" + prefix + ".parents")
    with open(prefix+".parents",'w') as fw:
        for iv,ik in parents.items():
            fw.write(iv+"\t"+"\t".join(ik)+"\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="get go class and alias from gene ontology in obo format")
    parser.add_argument('--go', dest='obo', type=str,
                        help="gene ontology in obo format,you could download it from http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo")
    parser.add_argument('--prefix', dest='prefix', type=str,default="go",
                        help="directory with prefix for output files, default is \"go\" ")
    localeArg = parser.parse_args()

    if localeArg.obo is not None and localeArg.prefix is not None:
        file=localeArg.obo
        prefix=localeArg.prefix
        parents,alias,gos,obsoletes=read_go_obo(file)
        write_alias(alias,obsoletes,prefix)
        write_class(gos,obsoletes,parents,prefix)
        logging.warning(time.asctime() +"\tdone!")
        get_parent(parents,prefix)
    else:
        print('Use python --help')
        sys.exit()