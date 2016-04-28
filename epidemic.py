wdir = '/home/valencianaplop/mrsa/analysis'

import os
import trile

def main():
	os.chdir(wdir)
	save = open('epidemic.fasta','w')
	read = open('heatmap.tsv')
	head = next(read).split()
	seq_dict = trile.read_fasta('all.fasta')
	bool_dict = dict()
	for line in read:
		col = line.split()
		name = col[0]
		col = col[1:]
		for i in range(len(col)):
			if float(col[i])>20.0:
				if head[i] not in bool_dict.keys():
					bool_dict[head[i]]=name
				else:
					bool_dict[head[i]]+='.'+name
	count = 1
	for key in bool_dict.keys():
		val = bool_dict.get(key)
		if val.count('O')==7 and 'N' not in val:
			save.write('>E'+str(count)+'\n')
			save.write(seq_dict.get(key)+'\n')
			count+=1

main()
			
