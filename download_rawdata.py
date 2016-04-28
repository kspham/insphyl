wdir = '/home/valencianaplop/mrsa/analysis'
import os

def main():
	id_list = os.path.join(wdir,'id_list.txt')
	id_list = open(id_list).read().splitlines()
	try:
		os.mkdir(os.path.join(wdir,'raw_data'))
	except:
		pass
	os.chdir(os.path.join(wdir,'raw_data'))
	for line in id_list:
		col = line.split('\t')
		file_list = col[1].replace('\n','').split(';')
		for file1 in file_list:
			os.system('wget -c '+file1)

main()
