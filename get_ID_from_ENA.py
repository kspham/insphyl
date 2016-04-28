wdir = '/home/valencianaplop/mrsa/analysis'
id_dict = {'MRSA_10C':'O5','MRSA_11C':'O6','MRSA_12C':'O7','MRSA_14C':'N1','MRSA_15C':'N2','MRSA_16B':'N3','MRSA_17B':'N4','MRSA_19B':'N6','MRSA_1B':'O1','MRSA_20B':'N7','MRSA_6C':'O2','MRSA_7C':'O3','MRSA_8C':'O4','MRSA_18B':'N5'}
import os

def main():
	info = os.path.join(wdir,'ERP001256.txt')
	info = open(info).read().splitlines()[1:]
	output = open(os.path.join(wdir,'id_list.txt'),'w')
	for line in info:
		col = line.split('\t')
		print col[7]
		output.write(col[7]+'\t'+col[8]+'\t'+id_dict.get(col[7])+'\n')

main()
