wdir = '/home/valencianaplop/mrsa/analysis'
import os
import glob
import trile

def main():
	os.chdir(wdir)
	read = open('all.fasta')
	save = open('all.fasta.tmp','w')
	count = 0
	heads,seqs,head,seq = list(),list(),list(),str()
	for line in read:
		if line!='' and line[0]=='>':
			if '|*|' in line:
				Id = '>I'+str(count)+'|'
				final_head = Id+'|'.join(head)
				heads.append(final_head)
				seqs.append(seq)
				count+=1
				head=[(line.split('|')[-1].split('.')[0])]
				seq = next(read).replace('\n','')
			else:
				head2 = line.split('|')[-1].split('.')[0]
				if head2 not in head:
					head.append(head2)
	heads = heads[1:]
	seqs = seqs[1:]
	for i in range(len(heads)):
		save.write(heads[i]+'\n'+seqs[i]+'\n')
	os.system('mv all.fasta.tmp all.fasta')

main()					
