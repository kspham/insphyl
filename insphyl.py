#!/usr/bin/env python

# Requirements:
# - Java (7.0)
# - SPAdes assembler (3.9.0)
# - QUAST (4.5)
# - Sibelia (3.0.6)

import os
import json
import urllib2 as u
import shutil
import multiprocessing as mp

def download_data(URL,pathOut,bSize = 8192):
	Resp = u.urlopen(URL,'rb'); Out = open(pathOut,'wb')
	TotalSize = int(Resp.info().getheaders('Content-length')[0]); Size = 0
	while True:
		Buffer = Resp.read(bSize); Size += len(Buffer)
		if not Buffer: break
		Out.write(Buffer)
		Status = '....%s %s/%s'%(pathOut,Size,TotalSize)
		Status += chr(8) * (len(Status)+1)
		print Status,
	print '\n',

def download_mrsa(pathIn,pathOut,pathOut2):
	try: os.mkdir(pathOut)
	except: pass
	In = open(pathIn,'rb'); next(In); Out = {}
	for Line in In:
		Cols = Line.strip().split('\t')
		Links = ['ftp://'+x for x in Cols[9].split(';')]
		FileNames = [os.path.join(pathOut,'%s_%s.fastq.gz'%(Cols[17].replace('.','_'),i)) for i in [1,2]]
		Out[Cols[17].replace('.','_')] = FileNames
		for i in range(2): download_data(Links[i],FileNames[i])
	json.dump(Out,open(pathOut2,'w'))
	return Out

def trim_data(RawData,pathOut,pathOut2):
	try: os.mkdir(pathOut)
	except: pass
	OldWD = os.getcwd(); os.chdir('Trimmomatic-0.36'); Out = {}
	for Name,Files in RawData.items():
		print '....%s'%Name
		Files = [os.path.join(OldWD,x) for x in Files]
		OutFiles = [os.path.basename(x) for x in Files]
		CMD = 'java -jar trimmomatic-0.36.jar PE -quiet %s %s .%s .unpair_1 .%s .unpair_2 ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(Files[0],Files[1],OutFiles[0],OutFiles[1])
		os.system(CMD)
		os.remove('.unpair_1'); os.remove('.unpair_2')
		Out[Name] = OutFiles
	os.chdir(OldWD)
	for Name,Files in Out.items():
		nFiles = [os.path.join(pathOut,File) for File in Files]
		for i in range(2): shutil.move('Trimmomatic-0.36/.'+Files[i],nFiles[i])
		Out[Name] = nFiles
	json.dump(Out,open(pathOut2,'w'))
	return Out

def assemble_data(TrimData,pathOut,pathOut2):
	try: os.mkdir(pathOut)
	except: pass
	nCPU = mp.cpu_count()/2; Out = {}
	for Name,Files in TrimData.items():
		print '....%s'%Name
		aData = os.path.join(pathOut,'%s.fasta'%Name)
		CMD = 'spades.py --careful -t %s -1 %s -2 %s -o .assembled_data > .log'%(nCPU,Files[0],Files[1])
		os.system(CMD)
		os.remove('.log')
		shutil.move('.assembled_data/scaffolds.fasta',aData)
		Out[Name] = aData
	json.dump(Out,open(pathOut2,'w'))
	shutil.rmtree('.assembled_data')
	return Out

def quast_data(aData,pathOut):
	CMD = 'quast.py -o %s '%(pathOut)
	CMD+= ' '.join(aData.values()) + '> .log'
	os.system(CMD)
	os.remove('.log')

def call_nonsyntenics(aData,pathRef,pathOut,pathOut2,MinLen=500):
	try: os.mkdir(pathOut)
	except: pass
	Out = {}
	for Name,File in aData.items():
		fileOut = Name+'.fasta'; nFileOut = os.path.join(pathOut,fileOut)
		CMD = 'C-Sibelia.py -m %s -o .call_nonsyntenics -u %s %s %s'%(MinLen,fileOut,pathRef,File)
		os.system(CMD)
		shutil.move('.call_nonsyntenics/%s'%fileOut,nFileOut)
		Out[Name] = nFileOut
	json.dump(Out,open(pathOut2,'w'))
	shutil.rmtree('.call_nonsyntenics')
	return Out

def format_nsr(NSR,pathOut,pathOut2):
	try: os.mkdir(pathOut)
	except: pass
	Out = {}
	for Name,File in NSR.items():
		nFile = os.path.join(pathOut,os.path.basename(File))
		Save = open(nFile,'w'); Count = 1
		for Line in open(File):
			if Line[0]=='>':
				Line = '>%s.%s'%(Name,Count)+'\n'
				Count += 1
			Save.write(Line)
		Save.close()
		Out[Name] = nFile
	json.dump(Out,open(pathOut2,'w'))
	return Out

def remove_duplicate(fNSR,pathOut,pathOut2):
	try: os.mkdir(pathOut)
	except: pass
	Out = {}
	for Name,File in fNSR.items():
		print '....%s'%Name
		nFile = os.path.join(pathOut,os.path.basename(File))
		CMD = 'uclust --quiet --sort %s --output .sorted_data.fasta'%File
		os.system(CMD)
		CMD = 'uclust --quiet --input .sorted_data.fasta --uc .index.uc'
		os.system(CMD)
		CMD = 'uclust --quiet --input .sorted_data.fasta --uc2fasta .index.uc --types S --output %s'%nFile
		os.system(CMD)
		Out[Name] = nFile
	for x in ['.sorted_data.fasta','.index.uc']: os.remove(x)
	json.dump(Out,open(pathOut2,'w'))
	return Out

def combine_data(Files,pathOut):
	Out = open(pathOut,'w')
	for File in Files:
		for Line in open(File): Out.write(Line)
	Out.close()

def combine_header(File):
	Out = {}
	for Line in open(File):
		if Line[0]=='>':
			Id,Status,Sample = Line.strip()[1:].split('|')
			Sample = Sample.split('.')[0]
			if Id in Out and Sample not in Out[Id]: Out[Id].append(Sample)
			elif Id not in Out: Out[Id] = [Sample]
	return Out

def compare_nsr(fuNSR,pathOut,pathOut2):
	combine_data(fuNSR.values(),'.merge.fasta')
	CMD = 'uclust --quiet --sort .merge.fasta --output .sorted_data.fasta'
	os.system(CMD)
	CMD = 'uclust --quiet --input .sorted_data.fasta --uc .index.uc'
	os.system(CMD)
	CMD = 'uclust --quiet --input .sorted_data.fasta --uc2fasta .index.uc --output .final.fasta'
	os.system(CMD)
	dHeader = combine_header('.final.fasta')
	Out = open(pathOut,'w'); In = open('.final.fasta'); dIns = {}; NameList = fuNSR.keys()+['reference']
	for Line in In:
		if Line[0]=='>':
			Id,Status,Sample = Line.strip()[1:].split('|')
			if Status=='*':
				nHead = '>Insertion_%s_%s'%(int(Id)+1,'|'.join(dHeader[Id]))
				Out.write(nHead+'\n'+next(In)+'\n')
				dIns['_'.join(nHead[1:].split('_')[:2])] = dict(zip(NameList,map(int,[(x in dHeader[Id]) for x in NameList])))
	Out.close()
	for x in ['.sorted_data.fasta','.merge.fasta','.index.uc','.final.fasta']: os.remove(x)
	json.dump(dIns,open(pathOut2,'w'))
	return dIns

def boolean_heatmap(pathIn,pathOut):
	CMD = 'Rscript plotting.R boolean_heatmap %s %s'%(pathIn,pathOut)
	os.system(CMD)

if __name__ == '__main__':
	print '..download data of MRSA'
	#RawData = download_mrsa('PRJEB2912.txt','mrsa_data','mrsa_data.json')
	RawData = json.load(open('mrsa_data.json'))

	print '..trim data of MRSA'
	#TrimData = trim_data(RawData,'mrsa_trimmed_data','mrsa_trimmed_data.json')
	TrimData = json.load(open('mrsa_trimmed_data.json'))

	print '..assemble data of MRSA'
	#aData = assemble_data(TrimData,'mrsa_assembled_data','mrsa_assembled_data.json')
	aData = json.load(open('mrsa_assembled_data.json'))

	print '..quast data of MRSA'
	#quast_data(aData,'mrsa_quast')

	print '..call nonsyntenic regions from data of MRSA'
	#NSR = call_nonsyntenics(aData,'HE681097.fasta','mrsa_nsr','mrsa_nsr.json')
	NSR = json.load(open('mrsa_nsr.json'))

	print '..format NSR of data of MRSA'
	#fNSR = format_nsr(NSR,'mrsa_formatted_nsr','mrsa_formatted_nsr.json')
	fNSR = json.load(open('mrsa_formatted_nsr.json'))

	print '..remove duplicate from NSR of MRSA'
	#uNSR = remove_duplicate(fNSR,'mrsa_unique_nsr','mrsa_unique_nsr.json')
	uNSR = json.load(open('mrsa_unique_nsr.json'))

	print '..reformat NSR of MRSA'
	#fuNSR = format_nsr(uNSR,'mrsa_reformatted_nsr','mrsa_reformatted_nsr.json')
	fuNSR = json.load(open('mrsa_reformatted_nsr.json'))

	print '..compare NSR of MRSA'
	#cNSR = compare_nsr(fuNSR,'mrsa_comparison.fasta','mrsa_comparison.json')
	cNSR = json.load(open('mrsa_comparison.json'))

	print '..plot booleanheatmap'
	boolean_heatmap('mrsa_comparison.json','mrsa_comparison_1.pdf')
