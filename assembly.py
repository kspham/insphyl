spades = '/home/valencianaplop/program/spades/bin/spades.py'
wdir = '/home/valencianaplop/mrsa/analysis'
import os

def main():
	try:
		os.mkdir(os.path.join(wdir,'assembly'))
	except:
		pass
	id_list = open(os.path.join(wdir,'id_list.txt')).read().splitlines()
	file_list = os.listdir(os.path.join(wdir,'trim_data'))
	for line in id_list:
		name = line.split('\t')[-1].replace('\n','')
		file_list = list()
		for file1 in os.listdir(os.path.join(wdir,'trim_data')):
			if name in file1:
				file_list.append(os.path.join(wdir,'trim_data',file1))
		cmd = spades + ' -o ' + os.path.join(wdir,'assembly',name)
		cmd+= ' -1 '+file_list[0]+' -2 '+file_list[1]
		cmd+= ' --careful'
		os.system(cmd)
		os.system('mv '+os.path.join(wdir,'assembly',name,'scaffolds.fasta')+' '+os.path.join(wdir,'assembly',name+'.fasta'))
		os.system('rm -r '+os.path.join(wdir,'assembly',name))

main()
