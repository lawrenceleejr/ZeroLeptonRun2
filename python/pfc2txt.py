#!/bin/env python

import os
import xml.etree.ElementTree as ET


PFCFILE='PoolFileCatalog.xml'
fout = open('pfc.txt','w')
if os.path.exists(PFCFILE):
    doc = ET.parse(PFCFILE)
    for File in doc.getroot():
        logicalname = ''
        physicalname = ''
        for item in File:
            if item.tag == 'logical' and len(item)>0:
                if item[0].attrib.has_key('name'): logicalname = item[0].attrib['name']
            if item.tag == 'physical':
                for info in item:
                    if info.attrib.has_key('name'): 
                        physicalname = os.path.basename(info.attrib['name'])
                        if ':' in physicalname:
                            physicalname = physicalname.split(':')[1]
        name = logicalname
        if name == '': name = physicalname
        if name != '':
            fout.write("%s %s\n" % (name,File.attrib['ID']))
fout.close()

