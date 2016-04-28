csibelia = '/home/valencianaplop/program/sibelia/bin/C-Sibelia.py'
wdir = '/home/valencianaplop/mrsa/analysis'
ref_file = 'HE681097.fasta'
mlen = 500
import os
import glob

def main():
	ref = os.path.join(wdir,ref_file)
	file_list = glob.glob(os.path.join(wdir,'assembly','*.fasta'))
	try:
		os.mkdir(os.path.join(wdir,'NSR'))
	except:
		pass
	os.chdir(os.path.join(wdir,'NSR'))
	for file1 in file_list:
		name = os.path.basename(file1).split('.fasta')[0]
		cmd = csibelia+' -o '+name
		cmd+= ' -u '+name+'.fasta'
		cmd+= ' -m '+str(mlen)
		cmd+= ' -v '+name+'.vcf'
		cmd+= ' '+ref +' '+file1
		os.system(cmd)
		os.system('mv '+os.path.join(name,name+'.fasta')+' ./')
		os.system('mv '+os.path.join(name,name+'.vcf')+' ./')
		os.system('rm -r '+name)

main()	
