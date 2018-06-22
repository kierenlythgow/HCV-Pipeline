#!/usr/bin/env bash

SNORKTEST=/well/bsg/microbial/analysis/stophcv/20170602-Snork-PHE/data/snorktest
rawdir=/well/bsg/microbial/analysis/stophcv/20170602-Snork-PHE/data/snorktest/raw
datsnorkdir=/well/bsg/microbial/analysis/stophcv/20170602-Snork-PHE/data/snorktest/dat/SnorkHCVT

for bam in `ls ${rawdir}/SnorkHCVT*.bam` ; do

    dataid=`basename ${bam} | sed "s,.bam,,g"`
    dataiddir=${datsnorkdir}/${dataid}

    if [ ! -d ${dataiddir} ] ; then
        mkdir -p ${dataiddir}
    fi
    scriptpath=${dataiddir}/${dataid}_snork_splitpop.sh
    echo "Creating ${scriptpath}"
    echo "#!/usr/bin/env bash" > ${scriptpath}
    echo >> ${scriptpath}
    echo "${SNORKTEST}/src6/snork.py splitpops \\" >> ${scriptpath}
    echo "-bin ${SNORKTEST}/src6 \\" >> ${scriptpath}
    echo "-profile None \\" >> ${scriptpath}
    echo "-config ${SNORKTEST}/src6/snork.config.wtchg \\" >> ${scriptpath}
    echo "-orgid Hepc \\" >> ${scriptpath}
    echo "-dataid ${dataid} \\" >> ${scriptpath}
    echo "-samplename ${dataid} \\" >> ${scriptpath}
    echo "-bampath ${SNORKTEST}/raw/${dataid}.bam \\" >> ${scriptpath}
    echo "-targetrefid wtchgR00000071 \\" >> ${scriptpath}
    echo "-targetrefpath ${SNORKTEST}/ref/wtchgR00000071/wtchgR00000071.fasta \\" >> ${scriptpath}
    echo "-outdir ${SNORKTEST}/dat/SnorkHCVT/${dataid} \\" >> ${scriptpath}
    echo "-logdir ${SNORKTEST}/dat/SnorkHCVT/${dataid} \\" >> ${scriptpath}
    echo "-overwrite False \\" >> ${scriptpath}
    echo "-deleteints True \\" >> ${scriptpath}
    echo "-verbosity INFO" >> ${scriptpath}
    chmod a+x ${scriptpath}
done
