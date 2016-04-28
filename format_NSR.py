wdir = '/home/valencianaplop/mrsa/analysis'
import os
import glob
import trile

def main():
	try:
		os.mkdir(os.path.join(wdir,'format_NSR'))
	except:
		pass
	file_list = glob.glob(os.path.join(wdir,'NSR','*.fasta'))
	for file1 in file_list:
		name = os.path.basename(file1).split('.fasta')[0]
		fasta = trile.read_fasta(file1)
		save = open(os.path.join(wdir,'format_NSR',name+'.fasta'),'w')
		count = 0
		for head in fasta.keys():
			count+=1
			save.write('>'+name+'.'+str(count)+'\n')
			save.write(fasta.get(head)+'\n')
		save.close()

main()
