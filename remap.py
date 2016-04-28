bwa = 'bwa'
samtools = 'samtools'
wdir = '/home/valencianaplop/mrsa/analysis'
import os
import glob

def remap(raw,ref,Id):
	name = os.path.basename(raw[0])
	cmd = bwa+' index '+ref
	os.system(cmd)

	for i in range(2):
		cmd = bwa+' aln '+ref+' '+raw[i]
		cmd+= ' > '+os.path.basename(raw[i]).split('.')[0]+'.sai'
		os.system(cmd)
	
	cmd = bwa+' sampe '+ref
	cmd+= ' '+os.path.basename(raw[0]).split('.')[0]+'.sai'
	cmd+= ' '+os.path.basename(raw[1]).split('.')[0]+'.sai'
	cmd+= ' '+raw[0]
	cmd+= ' '+raw[1]
	cmd+= ' > '+Id+'.sam'
	os.system(cmd)

	file_list = glob.glob(ref+'.*')+glob.glob('*.sai')
	os.system('rm '+' '.join(file_list))

	cmd = samtools+' faidx '+ref
	os.system(cmd)

	cmd = samtools+' import '+ref+'.fai '
	cmd+= Id+'.sam '+Id+'.bam'
	os.system(cmd)
	
	cmd = samtools+' sort '
	cmd+= Id+'.bam '+Id+'.sorted'
	os.system(cmd)
	
	cmd = samtools+' index '+Id+'.sorted.bam'
	os.system(cmd)

	os.system('rm '+Id+'.sam '+Id+'.bam')

def main():
	try:
		os.mkdir(os.path.join(wdir,'remap'))
	except:
		pass
	os.chdir(os.path.join(wdir,'remap'))
	ref = os.path.join(wdir,'all.fasta')
	id_list = open(os.path.join(wdir,'id_list.txt'))
	file_list = glob.glob(os.path.join(wdir,'trim_data','*'))
	for line in id_list:
		Id = line.split('\t')[-1].replace('\n','')
		raw = list()
		for file1 in file_list:
			if Id in file1:
				raw.append(file1)
		raw.sort()
		remap(raw,ref,Id)

main()
		
