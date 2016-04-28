trimmomatic = 'java -jar /home/valencianaplop/program/trimmomatic/trimmomatic-0.35.jar'
wdir = '/home/valencianaplop/mrsa/analysis'
import os
import glob

def main():
	try:
		os.mkdir(os.path.join(wdir,'trim_data'))
	except:
		pass
	os.chdir(os.path.join(wdir,'trim_data'))
	id_list = os.path.join(wdir,'id_list.txt')
	id_list = open(id_list).read().splitlines()
	for line in id_list:
		col = line.split('\t')
		file_list = col[1].replace('\n','').split(';')
		for i in range(len(file_list)):
			file_list[i]=os.path.basename(file_list[i])
			file_list[i]=os.path.join(wdir,'raw_data',file_list[i])
		p1 = col[-1].replace('\n','')+'_1.fastq.gz'
		u1 = col[-1].replace('\n','')+'_u1.fasta.gz'
		p2 = col[-1].replace('\n','')+'_2.fastq.gz'
		u2 = col[-1].replace('\n','')+'_u2.fastq.gz'
		cmd = ' '.join(file_list+[p1,u1,p2,u2])
		cmd = trimmomatic + ' PE -phred33 ' + cmd
		cmd+= ' ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
		os.system(cmd)
		os.system('rm '+u1+' '+u2)

main()
