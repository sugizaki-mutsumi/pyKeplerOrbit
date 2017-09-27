#!/usr/bin/env python

import os

parname_list = [
    "name",
    "ra",
    "dec",
    "binary",
    "pbinary",
    "pbdot",
    "epoch_type",
    "binaryepoch",
    "axsini",
    "periapse",
    "apsidalrate",
    "eccentricity",
    "egress",
    "ingress",
    ]

parlist_def = [
    "dummyname", 0.0, 0.0,
    "N", 0.0, 0.0,
    "P", 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0
    ]


### from text file
def readfile_OrbitalElements(fname) :
    #### default parlist
    parlist = parlist_def

    for line in open(fname).readlines() :
        if line[0]=="#" : continue
        if len(line.split())<3 : continue
        for ipar in range(len(parname_list)) :
            pname = parname_list[ipar]
            if line.split()[0]==pname :
                if ipar in [0, 3, 6] :
                    parlist[ipar] = line.split(",")[1].strip()
                else :
                    parlist[ipar] = float(line.split(",")[1])
                break
            
    return parlist



if __name__=="__main__" :


    parfname_list = [
        "orbpar/gx304m1_ms.param",
        "orbpar/groj1008_ms.param",
        "orbpar/4u0115_ms2.param",
        ]

    for fname in parfname_list :
        parlist = readfile_OrbitalElements(fname)
        if parlist[0]!="dummyname" :
            #print fname
            print parlist

            
            
