wdir = '/home/valencianaplop/mrsa/analysis'
import os

def main():
	os.chdir(wdir)
	read = open('all.fasta')
	save = open('all.tsv','w')
	for line in read:
		if line!='' and line[0]=='>':
			Id = line.split('|')[0][1:]
			save.write('\t'+Id)
	
	id_list = open('id_list.txt').read().splitlines()
	for line in id_list+['reference']:
		Id = line.split('\t')[-1]
		save.write('\n'+Id)
		read = open('all.fasta')
		score = 0
		for line2 in read:
			if line2!='' and line2[0]=='>':
				if '|'+Id in line2:
					score = 1
				else:
					score = 0
				save.write('\t'+str(score))
	save.close()

main()
