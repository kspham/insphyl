import os
import ftplib

wdir = '/home/trile/Documents/insphyl'

def main():
	os.chdir(wdir)
	os.makedirs('RawData')
	IdList = open('ERP001256.txt')
	for Line in IdList:
		Name,Links = Line.strip().split()[1:]
		Links = Links.split(';')
		for Link in Links:
			

main()
 
