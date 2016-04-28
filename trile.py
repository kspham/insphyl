import sys

def get_argv(key,type1,default):
	output = default
	for item in sys.argv:
		i = sys.argv.index(item)
		if '-'+key==item:
			if type1==int or type1==float:
				output = type1(sys.argv[i+1])
			elif type1==str:
				output = sys.argv[i+1]
			elif type1==list:
				output = list()
				for ii in range(i+1,len(sys.argv)):
					if sys.argv[ii][0]=='-':
						break
					else:
						output.append(sys.argv[ii])
			break
	return output

def read_fasta(file1):
	fasta_dict = dict()
	read = open(file1)
	seq_name,seq_content,seq = list(),list(),str()
	for line in read:
		if line!='' and line[0]=='>':
			seq_name.append(line[1:].replace('\n',''))
			seq_content.append(seq)
			seq = str()
		else:
			seq+= line.replace('\n','')
	seq_content.append(seq)
	seq_content = seq_content[1:]
	return dict(zip(seq_name,seq_content))
