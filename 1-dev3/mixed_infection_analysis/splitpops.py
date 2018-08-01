#!/usr/bin/env python

import os
import sys
import time

import snorkversion

_processname = 'splitpops'
_outfilesuffixL = [
    '.targetpop.readclass.bam',	# BAM file with reads tagged with RG:Z:readclass
    '.targetpop.readclass.txt',	# BAM file with reads tagged with RG:Z:readclass
    '.targetpop.stats.txt',     # TSV file of stats: dataid readclass readcount pct
    '.targetpop.sentinel.txt'	# Empty file indicating sub-tool finished
]
# Also produces _targetpop_POPID.bam(s)

def Prerequisites(args, P, mylogger, myhandler, processname):
    mylogger.debug('Checking pre-requisites')
    mylogger.debug('Checking pre-requisites finished')
    return 0

def Convert_Input_Reads_BAM_To_FASTQ(args, P, mylogger, myhandler, _processname):
    'Name-sort the BAM file, then convert to a pair of FASTQ files. Exit program on error.'
  # Output files
    fq1 = os.path.join(args.outdir, os.path.basename(args.bampath).replace('.bam', '_1.fastq'))
    fq2 = os.path.join(args.outdir, os.path.basename(args.bampath).replace('.bam', '_2.fastq'))
  # Exit if outfiles exist and no overwrite
    if os.path.exists(fq1) and os.path.getsize(fq1) > 0 and \
       os.path.exists(fq2) and os.path.getsize(fq2) > 0 and \
       not args.overwrite:
        mylogger.info('Skipping BAM to FASTQ conversion')
        return fq1, fq2
    mylogger.info('Converting input BAM to FASTQ')
  # Variables
    samtools = P.option['prog']['samtools']
    inbam = args.bampath
    tmpbamstempath = os.path.join(args.outdir,
        os.path.basename(args.bampath).replace('.bam', '_namesorted_tmp'))
    java = P.option['prog']['java']
    samtofastq = P.option['prog']['samtofastq']
    intermediateL = [tmpbamstempath+'.bam']
  # Sort BAM reads by name
    cmd = '{samtools} sort -n {inbam} {tmpbamstempath}'.format(
        samtools=samtools, inbam=inbam, tmpbamstempath=tmpbamstempath)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to name-sort BAM ({0})'.format(inbam))
        sys.exit(P.err_code('ErrorSysCall'))
  # Convert BAM to FASTQ, losing any readgroups associated with each BAM record
    cmd = '{java} -jar {samtofastq} INPUT={inbam} FASTQ={fq1} SECOND_END_FASTQ={fq2}'.format(
        java=java, samtofastq=samtofastq, inbam=inbam, fq1=fq1, fq2=fq2)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to convert BAM to FASTQ ({0}, {1})'.format(
            inbam, fq1.replace('_1.fastq', '[1|2].fastq')))
        mylogger.error('Command: {0}'.format(cmd))
        sys.exit(P.err_code('ErrorSysCall'))
  # Remove intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return FASTQ paths
    return fq1, fq2

def Map_Reads_To_Target_References(args, P, mylogger, myhandler, _processname, fq1, fq2):
    'Map the FASTQ reads to the single FASTA file of references. Exit program on error.'
  # Output files
    mapbam = os.path.join(args.outdir, os.path.basename(fq1).replace('_1.fastq', '.targetpop.mapped.bam'))
    mapbamstem = mapbam.replace('.bam', '')
  # Exit if outfiles exist and no overwrite
    if os.path.exists(mapbam) and os.path.getsize(mapbam) > 0 and \
       not args.overwrite:
        mylogger.info('Skipping read mapping to reference')
        return mapbam
    mylogger.info('Mapping reads to references of known genotype classification')
  # Variables
    samtools = P.option['prog']['samtools']
    mapsam = os.path.join(args.outdir, os.path.basename(fq1).replace('_1.fastq', '.targetpop.mapped.sam'))
    bwa = P.option['prog']['bwa']
    refpath = args.targetrefpath.replace('.fasta', '.fa')
    intermediateL = [mapsam]
  # Run "bwa mem" keeping shorter, non-primary split hits
    cmd = '{bwa} mem -M {refpath} {fq1} {fq2}'.format(
        bwa=bwa, refpath=refpath, fq1=fq1, fq2=fq2)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to map BAM to target refs ({0}, {1})'.format(
            fq1.replace('_1.fastq', '[1|2].fastq'), refpath))
        mylogger.error('{0}'.format(re))
        sys.exit(P.err_code('ErrorSysCall'))
    with open(mapsam, 'w') as sam_fp:
        sam_fp.write('{0}\n'.format(ro.rstrip('\n')))
  # Convert SAM output to BAM, excluding (-F) non-primary (256) and supplementary (2048) alignments
    #cmd = '{samtools} view -F 2304 -Sb -o {outbam} {insam}'.format(
    #    samtools=samtools, outbam=mapbam, insam=mapsam)
    cmd = '{samtools} view -Sb -F 2304 {insam} | {samtools} sort -n - {outbamstem}'.format(
        samtools=samtools, insam=mapsam, outbamstem=mapbamstem)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to convert SAM to BAM ({0}, {1})'.format(mapsam, mapbam))
        mylogger.error('Command: {0}'.format(cmd))
        sys.exit(P.err_code('ErrorSysCall'))
#  # Convert name-sorted BAM to SAM format (overwriting existing unsorted SAM of same name)
#    cmd = '{samtools} view -h {inbam}'.format(inbam=mapbam)
#    mylogger.debug('Command: {0}'.format(cmd))
#    rv, ro, re = P.sys_exec(cmd)
#    if rv != 0:
#        mylogger.error('Failed to convert BAM to SAM ({0}, {1})'.format(mapbam, mapsam))
#        mylogger.error('Command: {0}'.format(cmd))
#        sys.exit(P.err_code('ErrorSysCall'))
#    with open(mapsam, 'w') as out_fp:
#        out_fp.write('{0}\n'.format(ro.rstrip('\n')))
  # Remove intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return path of mapped BAM
    return mapbam

def Classify_Read_Pairs_To_Target_Genotypes(args, P, mylogger, myhandler, _processname, mapbam):
    '''
    Read the pairs of name-sorted reads.
    Genotype classification is from the mapped mate with the highest MAPQ.
    If neither mate is mapped the read pair is "unclassified".
    '''
  # Output files
    readclasspath = os.path.join(args.outdir, os.path.basename(mapbam).replace('.targetpop.mapped.bam', '.targetpop.readclass.txt'))
    readclassbam = os.path.join(args.outdir, os.path.basename(mapbam).replace('.targetpop.mapped.bam', '.targetpop.readclass.bam'))
    statspath = os.path.join(args.outdir, os.path.basename(mapbam).replace('.targetpop.mapped.bam', '.targetpop.stats.txt'))
  # Exit if outfiles exist and no overwrite
    if os.path.exists(readclasspath) and os.path.getsize(readclasspath) > 0 and \
       os.path.exists(readclassbam) and os.path.getsize(readclassbam) > 0 and \
       os.path.exists(statspath) and os.path.getsize(statspath) > 0 and \
       not args.overwrite:
        mylogger.info('Skipping read pair classification')
        return readclasspath
    mylogger.info('Classifying each read pair')
  # Variables
    samtools = P.option['prog']['samtools']
    samsubsetpath = readclasspath.replace('.targetpop.readclass.txt', 'targetpop.mapbam.subset')
    readclasssam = os.path.join(args.outdir, os.path.basename(mapbam).replace('.targetpop.mapped.bam', '.targetpop.readclass.sam'))
    intermediateL = [samsubsetpath, readclasssam]
  # Create intermediate file with just readid, refname, mapq
    cmd = '{samtools} view {inbam} | cut -f1,2,3,5'.format(
        samtools=samtools, inbam=mapbam)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to extract QNAME, FLAG, RNAME and MAPQ from BAM ({0})'.format(mapbam))
        mylogger.error('Command: {0}'.format(cmd))
        sys.exit(P.err_code('ErrorSysCall'))
    with open(samsubsetpath, 'w') as out_fp:
        out_fp.write('{0}\n'.format(ro.rstrip('\n')))
  # Read each pair of records (in name-sorted order) and output allocation.
  # During this pass, also collate the number of read pairs in each readclass.
    classlistcnt = {}
    with open(readclasspath, 'w') as out_fp:
        R = ['readid', 'readclass']
        out_fp.write('{0}\n'.format('\t'.join(R)))
        try:
            with open(samsubsetpath) as in_fp:
                for line1 in in_fp:
                    line2 = in_fp.next()
                    r1_readid, r1_flag, r1_rname, r1_mapq = line1.split('\t')
                    r2_readid, r2_flag, r2_rname, r2_mapq = line2.split('\t')
                    r1_ismapped = not int(r1_flag)&4
                    r2_ismapped = not int(r2_flag)&4
                    if not r1_ismapped and not r2_ismapped:
                        readclass = "unclassified"
                    elif not r1_ismapped and r2_ismapped:
                        readclass = r2_rname.split('_')[0]
                    elif r1_ismapped and not r2_ismapped:
                        readclass = r1_rname.split('_')[0]
                    else: # if both mapped, return gt of highest scoring mapping,
                          # if a tie and reads map to different genotypes return "ambiguous"
                        r1_mapq = int(r1_mapq)
                        r2_mapq = int(r2_mapq)
                        if r1_mapq > r2_mapq:
                            readclass = r1_rname.split('_')[0]
                        elif r1_mapq < r2_mapq:
                            readclass = r2_rname.split('_')[0]
                        else:
                            r1_readclass = r1_rname.split('_')[0]
                            r2_readclass = r2_rname.split('_')[0]
                            if r1_readclass == r2_readclass:
                                readclass = r1_readclass
                            else:
                                readclass = "ambiguous"
                    R = [r1_readid, readclass]
                    out_fp.write('{0}\n'.format('\t'.join(R)))
                    if readclass not in classlistcnt:
                        classlistcnt[readclass] = 0
                    classlistcnt[readclass] += 1
        except StopIteration:
            pass
  # Print allocations with >args.mintargetpoppct% (typically 1.0) of reads in the
  # .targetpop.list.txt in decreasing order of frequency
    rpscnt = sum([classlistcnt[key] for key in classlistcnt.keys()])
    minrpscnt = int(args.mintargetpoppct / 100.0 * rpscnt)
    keyL_decr = sorted(classlistcnt, key=classlistcnt.get, reverse=True)
    #readclassL = []
    with open(statspath, 'w') as out_fp:
        R = ['readclass', 'readcnt', 'readpct']
        out_fp.write('{0}\n'.format('\t'.join(R)))
        otherclasscnt = 0
        for readclass in keyL_decr:
            if readclass not in ['ambiguous', 'unclassified']:
                if classlistcnt[readclass] > minrpscnt:
                    R = [readclass, classlistcnt[readclass]*2,
                        round(classlistcnt[readclass] / float(rpscnt) * 100.0, 1)]
                    out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
                    #readclassL.append(readclass)
                else:
                    otherclasscnt += classlistcnt[readclass]
        if otherclasscnt:
            R = ['othertargetclasses', otherclasscnt, round(otherclasscnt / float(rpscnt) * 100.0, 1)]
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
            #readclassL.append('othertargetclasses')
        if 'ambiguous' in keyL_decr:
            R = ['ambiguous', classlistcnt['ambiguous']*2,
                round(classlistcnt['ambiguous'] / float(rpscnt) * 100.0, 1)]
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
            #readclassL.append('ambiguous')
        if 'unclassified' in keyL_decr:
            R = ['unclassified', classlistcnt['unclassified']*2,
                round(classlistcnt['unclassified'] / float(rpscnt) * 100.0, 1)]
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in R])))
            #readclassL.append('unclassified')
  # Create a new file .targetpop.readclass.bam with all the mapped reads from
  # .targetpop.mapped.bam annotated with the inferred readclass in the readgroup.
  # Retrieve existing bam header
    cmd = '{samtools} view -H {inbam}'.format(samtools=samtools, inbam=mapbam) 
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to retrieve header from BAM ({0})'.format(mapbam))
        mylogger.error('Command: {0}'.format(cmd))
        sys.exit(P.err_code('ErrorSysCall'))
    mapbamheader = ro.rstrip('\n')
  # Retrieve list of readclasses
    cmd = 'tail -n +2 {readclasspath} | cut -f2 | sort -u'.format(readclasspath=readclasspath)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to retrieve unique list of readclasses ({0})'.format(readclasspath))
        mylogger.error('Command: {0}'.format(cmd))
        sys.exit(P.err_code('ErrorSysCall'))
    readclassL = ro.strip('\n').split('\n') if len(ro.strip('\n')) else []
    readclassL.sort()
  # Retrieve entire file of readclass allocations for each readid
    with open(readclasspath, 'r') as in_fp:
        readclassD = dict(x.rstrip().split(None, 1) for x in in_fp)
  # Print all this info to a SAM file
    with open(readclasssam, 'w') as out_fp:
      # Append existing header from mapped bam containing @HD and @SQ records
        out_fp.write('{0}\n'.format(mapbamheader))
      # Add one @RG line for each readclass in the .targetpop.readclass.txt file
        for readclass in readclassL:
            out_fp.write('@RG\tID:{rgid}\n'.format(rgid=readclass))
      # Iterate through records of the mapped BAM file and output record appended with RG:Z:rgid
        cmd = '{samtools} view {inbam}'.format(samtools=samtools, inbam=mapbam)
        mylogger.debug('Command: {0}'.format(cmd))
        rv, ro, re = P.sys_exec(cmd)
        if rv != 0:
            mylogger.error('Failed to retrieve contents of BAM ({0})'.format(mapbam))
            mylogger.error('Command: {0}'.format(cmd))
            sys.exit(P.err_code('ErrorSysCall'))
        recordL = ro.rstrip('\n').split('\n')[1:]
        for record in recordL:
            R = record.split('\t')
            readid = R[0]
            R.append('RG:Z:{rgid}'.format(rgid=readclassD[readid]))
            out_fp.write('{0}\n'.format('\t'.join(R)))
  # Convert the SAM file to a (name-sorted) BAM file
    cmd = '{samtools} view -Sb -o {outbam} {insam} '.format(
        samtools=samtools, outbam=readclassbam, insam=readclasssam)
    mylogger.debug('Command: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to convert SAM to BAM ({0}, {1})'.format(readclasssam, readclassbam))
        mylogger.error('Command: {0}'.format(cmd))
        sys.exit(P.err_code('ErrorSysCall'))
  # Remove intermediate files
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Return path of file containing classification of each read pair id
    return readclasspath

def Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, fq1, fq2, mapbam):
    mylogger.info('Tidying up')
  # Delete any remaining intermediate files
    intermediateL = [fq1, fq2, mapbam]
    if args.deleteints:
        for path in intermediateL:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except:
                mylogger.warning('Failed to remove intermediate file ({0})'.format(path))
  # Create the sentinel path if it doesn't already exist
    sentinel_path = os.path.join(args.outdir, os.path.basename(args.bampath).replace('.bam','.targetpop.sentinel.txt'))
    if not os.path.exists(sentinel_path) or args.overwrite:
        now = time.localtime()
        with open(sentinel_path, 'w') as out_fp:
            out_fp.write('{0}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', now)))

    return 0

def Process(args, P, mylogger, myhandler, processname, suffixL):
    mylogger.info('Processing')
    alreadydone = True
    for suffix in suffixL:
        outpath = os.path.join(args.outdir, os.path.basename(args.bampath).replace('.bam', suffix))
        if not os.path.exists(outpath) or os.path.getsize(outpath) == 0:
            alreadydone = False
    if alreadydone and not args.overwrite:
        mylogger.info('Skipping processing as all output files already exist and are non-empty')
        return 0
    fq1, fq2 = Convert_Input_Reads_BAM_To_FASTQ(args, P, mylogger, myhandler, _processname)
    mapbam = Map_Reads_To_Target_References(args, P, mylogger, myhandler, _processname, fq1, fq2)
    readclasspath = Classify_Read_Pairs_To_Target_Genotypes(args, P, mylogger, myhandler, _processname, mapbam)
    Tidy_Up_Before_Finishing(args, P, mylogger, myhandler, _processname, fq1, fq2, mapbam)
    mylogger.debug('Processing finished')
    return 0

def main(parser, args, P, mylogger, myhandler, argv):
    mylogger.info('Started')
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname, _outfilesuffixL)
    mylogger.info('Finished')
    return 0
