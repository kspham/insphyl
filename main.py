#!/usr/bin/env python

import os
import urllib2
import json
import sys
import shutil
from collections import OrderedDict

###########################
### Support functions   ###
###########################

def read_fasta(path_file, name):
    ### Read a fasta file
    header = None
    count = 1
    for line in open(path_file):
        if line.startswith('>'):
            if header:
                header = '%s_%s' % (name, count)
                count += 1
                yield header, seq
            header = line[1:].strip()
            seq = ''
        else:
            seq += line.strip()
    if header:
        header = '%s_%s' % (name, count)
        count += 1
        yield header, seq



def chop(seq, width):
    ### Break a long string into pieces
    new_seq = []
    for i in xrange(0, len(seq), width):
        new_seq.append(seq[i:i+width])
    return '\n'.join(new_seq)



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
                print '[WARNING] %s was already downloaded' % name

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
            print '[WARNING] %s is already trimmed' % name
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
        print tag, name
        path_contig = os.path.join(dir_assem, '%s.fasta' % name)
        if os.path.exists(path_contig):
            print '[WARNING] %s is already assembled' % name
        else:
            assembly(info['trim'], path_contig, path_spades)
        info['contig'] = path_contig
        new_meta[name] = info
    return new_meta



###########################
### Call non-syntenic   ###
###########################

def call_nsr(path_contig, path_nsr, path_ref, path_sibelia):
    ### Call NSR of a sample
    dir_tmp, path_tmp = '.nsr_tmp', 'tmp.fasta'
    cmd = '%s -o %s -u %s ' % (path_sibelia, dir_tmp, path_tmp)
    cmd += '%s %s' % (path_ref, path_contig)
    os.system('%s >/dev/null 2>&1' % cmd)
    path_tmp = os.path.join(dir_tmp, path_tmp)
    shutil.move(path_tmp, path_nsr)
    shutil.rmtree(dir_tmp)



def call_nsr_wrapper(meta, dir_nsr, path_ref, path_sibelia):
    ### Call NSR of all samples
    tag = '[call_nsr_wrapper]'
    try_mkdir(dir_nsr)
    new_meta = dict()
    for name, info in meta.iteritems():
        print tag, name
        path_nsr = os.path.join(dir_nsr, '%s.fasta' % name)
        if os.path.exists(path_nsr):
            print '[WARNING] NSRs of %s were already called' % name
        else:
            call_nsr(info['contig'], path_nsr, path_ref, path_sibelia)
        info['nsr'] = path_nsr
        new_meta[name] = info
    return new_meta



###########################
### Remove duplicate    ###
###########################

def format_fasta(path_in, path_out, prefix='ins'):
    ### Format the header of a fasta file
    count = 1
    out = open(path_out, 'w')
    for line in open(path_in):
        if line.startswith('>'):
            header = '%s.%s' % (prefix, count)
            count += 1
            line = '>%s\n' % header
        out.write(line)
    out.close()



def rm_dup(path_nsr, path_out, path_uclust):
    ### Remove duplicate of a sample
    path_uc = '.nsr.uc'
    cmd = '%s --input %s --uc %s ' % (path_uclust, path_nsr, path_uc)
    cmd += '--id 0.9 --usersort'
    os.system('%s >/dev/null 2>&1' % cmd)
    path_fa = '.nu.fasta'
    cmd = '%s --uc2fasta %s --input %s' % (path_uclust, path_uc, path_nsr)
    cmd += ' --output %s --types S' % (path_fa)
    os.system('%s >/dev/null 2>&1' % cmd)
    os.remove(path_uc)
    format_fasta(path_fa, path_out)
    os.remove(path_fa)


def rm_dup_wrapper(meta, dir_out, path_uclust):
    ### Remove duplicate of all samples
    try_mkdir(dir_out)
    tag = '[rm_dup_wrapper]'
    new_meta = dict()
    for name, info in meta.iteritems():
        print tag, name
        path_out = os.path.join(dir_out, '%s.fasta' % name)
        if os.path.exists(path_out):
            print '[WARNING] %s was already filtered for duplicates' % name
        else:
            rm_dup(info['nsr'], path_out, path_uclust)
        info['ins'] = path_out
        new_meta[name] = info
    return new_meta

################
### Merge    ###
################

def cluster_data(path_all, path_cluster, path_uclust):
    ### Cluster all sequences into groups
    path_uc = '.all.uc'
    cmd = '%s --input %s --uc %s ' % (path_uclust, path_all, path_uc)
    cmd += '--id 0.9 --usersort'
    os.system('%s >/dev/null 2>&1' % cmd)
    cmd = '%s --uc2fasta %s --input %s ' % (path_uclust, path_uc, path_all)
    cmd += '--output %s' % path_cluster
    os.system('%s >/dev/null 2>&1' % cmd)
    os.remove(path_uc)



def merge_fasta(meta, path_out):
    ### Merge all data into a single file
    out = open(path_out, 'w')
    for name, info in meta.iteritems():
        for header, seq in read_fasta(info['ins'], name):
            out.write('>%s\n%s\n' % (header, chop(seq, 80)))
    out.close()



def parse_headers(path_cluster, path_ins):
    ### Parse headers of clustered sequences
    data = open(path_cluster)
    cluster = OrderedDict()
    for line in data:
        line = line[1:].strip()
        seq = next(data).strip()
        index, score, header = line.split('|', 2)
        index = int(index)
        sample = '_'.join(header.split('_')[:-1])
        if score == '*':
            cluster[index + 1] = {'samples': [sample], 'sequence': seq}
        else:
            cluster[index + 1]['samples'].append(sample)
    return cluster



def merge_data(meta, dir_result, path_uclust):
    ### Merge data and obtain correlation between samples
    try_mkdir(dir_result)
    path_ins = os.path.join(dir_result, 'cluster.json')
    tag = '[merge_data]'
    if os.path.exists(path_ins):
        print '[WARNING] Data was already merged'
    else:
        # Merge
        print tag, 'Merge'
        path_all = '.all.fasta'
        merge_fasta(meta, path_all)
        # Cluster
        print tag, 'Cluster'
        path_cluster = '.clusters.fasta'
        cluster_data(path_all, path_cluster, path_uclust)
        # Parse headers
        print tag, 'Parse'
        cluster = parse_headers(path_cluster, path_ins)
        json.dump(cluster, open(path_ins, 'w'), indent = 4, sort_keys = True)
        # Clean up
        for f in [path_all, path_cluster]:
            os.remove(f)
    return path_ins



#######################
### Parse boolean   ###
#######################

def plot_boolean(meta_json, path_bcor, path_rmain, dir_result):
    ### Visualize boolean heatmap
    path_plot = os.path.join(dir_result, 'boolean.pdf')
    cmd = 'Rscript %s plot_boolean %s %s %s' % (path_rmain,
                                                meta_json, path_bcor, path_plot)
    os.system('%s >/dev/null 2>&1' % cmd)



def parse_boolean(meta, path_ins, dir_result, path_rmain):
    ### Generate correlation from clusters
    tag = '[parse_boolean]'
    ins = json.load(open(path_ins), object_pairs_hook=OrderedDict)
    all_sam = meta.keys()
    path_bcor = os.path.join(dir_result, 'cor.csv')
    if os.path.exists(path_bcor):
        print '[WARNING] Boolean correlation is already genereated'
    else:
        # Create a boolean matrix
        print tag, 'boolean matrix'
        out = open(path_bcor, 'w')
        out.write('%s\n' % ','.join([''] + all_sam))
        for _id, info in ins.iteritems():
            row = [x in info['samples'] for x in all_sam]
            row = map(str, map(int, row))
            row = [_id] + row
            out.write('%s\n' % ','.join(row))
        out.close()
        # Visualize boolean matrix
        print tag, 'visualize'
        meta_json = '.meta.json'
        json.dump(meta, open(meta_json, 'w'))
        plot_boolean(meta_json, path_bcor, path_rmain, dir_result)
        os.remove(meta_json)
        return path_bcor



###############
### Main    ###
###############

def main(path_code, metafile, path_ref, dir_output):
    ### Main
    # Input
    # metafile = 'PRJEB2912.csv'
    # dir_output = 'mrsa_process'
    # path_ref = 'HE681097.fasta'
    # Set up
    path_spades = 'spades.py'
    path_sibelia = 'C-Sibelia.py'
    path_uclust = 'uclust'
    dir_download = os.path.join(dir_output, 'download')
    dir_trim = os.path.join(dir_output, 'trim')
    dir_assem = os.path.join(dir_output, 'assembly')
    dir_nsr = os.path.join(dir_output, 'nsr')
    dir_fnsr = os.path.join(dir_output, 'insertions')
    dir_result = os.path.join(dir_output, 'result')
    path_rmain = os.path.join(os.path.dirname(path_code),
                                'r_main.R')
    path_trimmo = os.path.join(os.path.dirname(path_code),
                                'Trimmomatic-0.36/trimmomatic-0.36.jar')
    path_adapter = os.path.join(os.path.dirname(path_code),
                                'Trimmomatic-0.36/adapters/TruSeq3-PE.fa')
    # From fastq to NSR
    try_mkdir(dir_output)
    meta = process_metadata(metafile)
    meta = download_wrapper(meta, dir_download)
    meta = trim_wrapper(meta, dir_trim, path_trimmo, path_adapter)
    meta = assembly_wrapper(meta, dir_assem, path_spades)
    meta = call_nsr_wrapper(meta, dir_nsr, path_ref, path_sibelia)
    meta = rm_dup_wrapper(meta, dir_fnsr, path_uclust)
    # Post-NSR
    path_ins = merge_data(meta, dir_result, path_uclust)
    path_bcor = parse_boolean(meta, path_ins, dir_result, path_rmain)
    print "Success"



if __name__ == '__main__':
    main(*sys.argv)
