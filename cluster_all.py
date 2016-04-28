uclust = 'uclustq1.2.22_i86linux64'
wdir = '/home/valencianaplop/mrsa/analysis'
threshold = 90
import os
import glob
import trile

def main():
	os.chdir(wdir)
	file_list = glob.glob(os.path.join(wdir,'cluster_NSR','*.fasta'))
	save=open('all.fasta','w')
	for file1 in file_list:
		fasta = trile.read_fasta(file1)
		for head in fasta.keys():
			save.write('>'+head+'\n')
			save.write(fasta.get(head)+'\n')
	save.close()

	file1 = 'all.fasta'
	name = 'all'
        cmd = uclust+' --sort '+file1
        cmd+= ' --output '+file1+'.tmp'
        os.system(cmd)

        cmd = uclust+' --input '+file1+'.tmp'
        cmd+= ' --uc '+name+'.uc'
        os.system(cmd)

        cmd = uclust+' --input '+name+'.fasta --uc2fasta '
        cmd+= name+'.uc --output '
        cmd+= name+'.fasta.tmp2'
        os.system(cmd)
        os.system('mv '+name+'.fasta.tmp2 '+name+'.fasta')
        os.system('rm '+name+'.uc '+name+'.fasta.tmp')

main()
