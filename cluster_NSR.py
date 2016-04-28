uclust = 'uclustq1.2.22_i86linux64'
wdir = '/home/valencianaplop/mrsa/analysis'
threshold = 90
import os
import glob
import trile

def main():
	try:
		os.mkdir(os.path.join(wdir,'cluster_NSR'))
	except:
		pass
	os.chdir(os.path.join(wdir,'cluster_NSR'))
	file_list = glob.glob(os.path.join(wdir,'format_NSR','*.fasta'))
	for file1 in file_list:
		name = os.path.basename(file1).split('.fasta')[0]
		cmd = uclust+' --sort '+file1
		cmd+= ' --output '+os.path.basename(file1)
		os.system(cmd)

		cmd = uclust+' --input '+os.path.basename(file1)
		cmd+= ' --uc '+name+'.uc'
		os.system(cmd)

		cmd = uclust+' --input '+name+'.fasta --uc2fasta '
		cmd+= name+'.uc --types S --output '
		cmd+= name+'.fasta.tmp'
		os.system(cmd)
		os.system('mv '+name+'.fasta.tmp '+name+'.fasta')
		os.system('rm '+name+'.uc')
		
                fasta = trile.read_fasta(name+'.fasta')
                save = open(name+'.fasta','w')
                count = 0
                for head in fasta.keys():
                        count+=1
                        save.write('>'+name+'.'+str(count)+'\n')
                        save.write(fasta.get(head)+'\n')
                save.close()


main()
