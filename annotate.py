wdir = '/home/valencianaplop/mrsa/analysis'
prokka = '/home/valencianaplop/program/prokka/bin/prokka'

import os

def main():
	os.chdir(wdir)
	cmd = prokka
	cmd+= ' --force --quiet --outdir annotate '
	cmd+= 'epidemic.fasta'
	os.system(cmd)

main()
