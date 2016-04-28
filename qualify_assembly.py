quast = '/home/valencianaplop/program/quast/quast.py'
wdir = '/home/valencianaplop/mrsa/analysis'
import os

def main():
	id_list = open(os.path.join(wdir,'id_list.txt')).read().splitlines()
	assembly = list()
	for line in id_list:
		name = line.split('\t')[-1].replace('\n','')
		assembly.append(os.path.join(wdir,'assembly',name+'.fasta'))

	cmd = ' -o '+os.path.join(wdir,'quast_evaluation')
	cmd = quast + cmd + ' '+ ' '.join(assembly)
	os.system(cmd)

main()
