#!/usr/bin/env bash

SNORKTEST=/well/bsg/microbial/analysis/stophcv/20170602-Snork-PHE/data/snorktest
rawdir=/well/bsg/microbial/analysis/stophcv/20170602-Snork-PHE/data/snorktest/raw
datsnorkdir=/well/bsg/microbial/analysis/stophcv/20170602-Snork-PHE/data/snorktest/dat/SnorkHCVT

for bam in `ls ${rawdir}/SnorkHCV*.bam` ; do

    dataid=`basename ${bam} | sed "s,.bam,,g"`
    dataiddir=${datsnorkdir}/${dataid}

    if [ ! -d ${dataiddir} ] ; then
        mkdir -p ${dataiddir}
    fi
    scriptpath=${dataiddir}/${dataid}_snork_splitpop.sh
    shlogpath=${dataiddir}/${dataid}_snork_splitpop.sh.log
    sentinelpath=${dataiddir}/${dataid}.targetpop.stats.txt
    if [ ! -s ${sentinelpath} ] ; then
        echo "Running ${scriptpath}"
        ${scriptpath} &> ${shlogpath}
    else
        echo "Skipping ${scriptpath}"
    fi
done
