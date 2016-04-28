genomeCoverageBed = 'genomeCoverageBed'
mergeBed = 'mergeBed'
wdir = '/home/valencianaplop/mrsa/analysis'

import os
import glob

def create_len_dict(file1):
	read = open(file1)
	len_dict = dict()
	for line in read:
		if line!='' and line[0]=='>':
			len_dict[line[1:-1]]=len(next(read))-1
	return len_dict

def main():
	os.chdir(wdir)
	try:
		os.mkdir('coverage')
	except:
		pass

	len_dict = create_len_dict('all.fasta')
	sort_head = len_dict.keys()
	sort_head.sort()
	os.chdir(os.path.join(wdir,'coverage'))	
	file_list=glob.glob(os.path.join(wdir,'remap','*.sorted.bam'))

	save = open(os.path.join(wdir,'heatmap.tsv'),'w')
	for head in sort_head:
		save.write('\t'+head)

	for file1 in file_list:
		name = os.path.basename(file1).split('.')[0]
		cmd = genomeCoverageBed
		cmd+= ' -ibam '+file1
		cmd+= ' -bg > '+name+'.bed'
		os.system(cmd)

		cmd = mergeBed
		cmd+= ' -i '+name+'.bed'
		cmd+= ' > '+name+'2.bed'
		os.system(cmd)

		os.system('mv '+name+'2.bed '+name+'.bed')
		
		save.write('\n'+name)
		read = open(name+'.bed')
		base,cov,cov_dict = str(),int(),dict()
		for line in read:
			col = line.split('\t')
			if col[0]!=base:
				if base in sort_head:
					cov = float(cov)/len_dict.get(base)*100
					cov_dict[base]=cov

				cov = int(col[2])-int(col[1])
				base = col[0]
			else:
				cov+= int(col[2])-int(col[1])
		for head in sort_head:
			if head in cov_dict.keys():
				save.write('\t'+str(cov_dict.get(head)))
			else:
				save.write('\t0')
		

main()
