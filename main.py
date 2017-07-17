#!/usr/bin/env python
# Contact: Tri Le - lequangminhtri@gmail.com

import os
import urllib2
import json
import sys
import shutil

###########################
### Support functions   ###
###########################

def process_metadata(metafile):
    ### Format metadata into a json of samples
    # Set up
    name_tag = 'sample_alias'
    status_tag = 'status'
    fastq_tag = 'fastq_ftp'
    # Process
    data = [x.strip().split('\t') for x in open(metafile).readlines()]
    headers, data = data[0], zip(*data[1:])
    ddata = dict(zip(headers, data))
    meta = dict()
    for i in xrange(len(ddata[name_tag])):
        name = ddata[name_tag][i]
        info = dict()
        info['fastq'] = ['ftp://%s' % (x) \
                        for x in ddata[fastq_tag][i].split(';')]
        info['status'] =  ddata[status_tag][i]
        meta[name] = info
    return meta



def try_mkdir(dirname):
    # Create a directory
    if os.path.exists(dirname):
        print '[WARNING] %s exists' % dirname
    else:
        os.mkdir(dirname)



###################
### Download    ###
###################

def download(url, filename):
    ### Download a file
    tag = '[downnload]'
    size = 1024 * 8
    request = urllib2.urlopen(url)
    text = request.read(size)
    print tag, filename
    with open(filename, 'w') as f:
        while text:
            f.write(text)
            text = request.read(size)



def download_wrapper(meta, dirname):
    ### Download all data in the meta
    try_mkdir(dirname)
    tag = '[download_wrapper]'
    new_meta = dict()
    for name, info in meta.iteritems():
        print tag, name
        filenames = list()
        for i in range(2):
            filename = os.path.join(dirname, '%s_%s.fastq.gz' % (name, i + 1))
            filenames.append(filename)
            if not os.path.exists(filename):
                download(info['fastq'][i], filename)
            else:
                print '[WARNING] %s exists' % filename

        info['data'] = filenames
        meta[name] = info
    return meta



###############
### Trim    ###
###############

def trim(fastqfiles, outfiles, path_tool, path_adapter):
    ### Trim a paired-end data with trimmomatic
    unpaired_files = ['.tmp.fastq.gz', '.tmp_2.fastq.gz']
    cmd = 'java -jar %s PE %s %s %s %s %s %s' % (path_tool, fastqfiles[0],
                            fastqfiles[1], outfiles[0], unpaired_files[0],
                            outfiles[1], unpaired_files[1])
    cmd += ' ILLUMINACLIP:%s:2:30:10 ' % path_adapter
    cmd += 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
    os.system('%s >/dev/null 2>&1' % cmd)
    for f in unpaired_files:
    	os.remove(f)



def trim_wrapper(meta, dir_output, path_tool, path_adapter):
    ### Trim all data from a given metadata of multiple samples
    tag = '[trim_wrapper]'
    try_mkdir(dir_output)
    new_meta = dict()
    for name, info in meta.iteritems():
        print tag, name
        outfiles = [os.path.join(dir_output,
                    '%s_%s.fastq.gz' % (name, x + 1)) for x in range(2)]
        if os.path.exists(outfiles[0]) & os.path.exists(outfiles[1]):
            print '[WARNING] %s is trimmed' % name
        else:
            trim(info['data'], outfiles, path_tool, path_adapter)
        info['trim'] = outfiles
        new_meta[name] = info
    return new_meta



###################
### Assemble    ###
###################

def assembly(path_fastq, path_contig, path_spades):
    ### Assemble a sample
    dir_tmp = '.assem_tmp'
    cmd = '%s --careful ' % path_spades
    cmd += '-1 %s -2 %s -o %s' % (path_fastq[0], path_fastq[1], dir_tmp)
    os.system('%s >/dev/null 2>&1' % cmd)
    shutil.move(os.path.join(dir_tmp, 'scaffolds.fasta'), path_contig)
    shutil.rmtree(dir_tmp)



def assembly_wrapper(meta, dir_assem, path_spades):
    ### Assemble all data
    tag = '[assembly_wrapper]'
    try_mkdir(dir_assem)
    new_meta = dict()
    for name, info in meta.iteritems():
        path_contig = os.path.join(dir_assem, '%s.fasta' % name)
        if os.path.exists(path_contig):
            print '[WARNING] %s is assembled' % name
        else:
            print tag, name
            assembly(info['trim'], path_contig, path_spades)
        info['contig'] = path_contig
        new_meta[name] = info
    return info



###########################
### Call non-syntenic   ###
###########################

def call_nsr(path_contig, path_nsr, path_ref, path_sibelia):
    ### Call NSR of a sample
    dir_tmp, path_tmp = '.nsr_tmp', 'tmp.fasta'
    cmd = '%s -o %s -v %s ' % (path_sibelia, dir_tmp, path_tmp)
    cmd += '%s %s' % (path_ref, path_contig)
    os.system('%s >/dev/null 2>&1' % cmd)
    shutil.move(path_tmp, path_nsr)
    shutil.rmtree(dir_tmp)



def call_nsr_wrapper(meta, dir_nsr, path_ref, path_sibelia):
    ### Call NSR of all samples
    tag = '[call_nsr_wrapper]'
    try_mkdir(dir_nsr)
    new_meta = dict()
    for name, info in meta.iteritems():
        path_nsr = os.path.join(dir_nsr, '%s.fasta' % name)
        if os.path.exists(path_nsr):
            print '[WARNING] non-syntenic regions of %s were called' % name
        else:
            print tag, name
            call_nsr(info['contig'], path_nsr, path_ref, path_sibelia)
        info['nsr'] = path_nsr
        new_meta[name] = info
    return new_meta



###############
### Main    ###
###############

def main(path_code):
    ### Main
    # Input
    metafile = 'PRJEB2912.csv'
    dir_output = 'mrsa_process'
    path_ref = 'HE681097.fasta'
    # Set up
    path_spades = 'spades.py'
    path_sibelia = 'C-Sibelia.py'
    dir_download = os.path.join(dir_output, 'download')
    dir_trim = os.path.join(dir_output, 'trim')
    dir_assem = os.path.join(dir_output, 'assembly')
    dir_nsr = os.path.join(dir_output, 'nsr')
    path_trimmo = os.path.join(os.path.dirname(path_code),
                                'Trimmomatic-0.36/trimmomatic-0.36.jar')
    path_adapter = os.path.join(os.path.dirname(path_code),
                                'Trimmomatic-0.36/adapters/TruSeq3-PE.fa')
    # Process
    try_mkdir(dir_output)
    meta = process_metadata(metafile)
    meta = download_wrapper(meta, dir_download)
    meta = trim_wrapper(meta, dir_trim, path_trimmo, path_adapter)
    meta = assembly_wrapper(meta, dir_assem, path_spades)
    meta = call_nsr_wrapper(meta, dir_nsr, path_ref, path_sibelia)



if __name__ == '__main__':
    main(*sys.argv)