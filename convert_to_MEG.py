wdir = '/home/valencianaplop/mrsa/analysis'
import os

def Matrix2Mega(path, outpath): #Thanks Tan
        def comp(row1, row2):
                d = []
                result = 0
                for trace in range(1, len(row1)):
                        r1 = float(row1[trace])
                        r2 = float(row2[trace])
                        d.append(abs(r1-r2))
                for diff in d:
                        result+=diff #xem lai float va int
                return result
        def read_matrix(path):
                matrix = open(path)
                data = matrix.readlines()[1:]
                records = []
                for row in data:
                        record = row.split('\t')
                        records.append(record[:])
                return records
        def write_mega(hash,outpath):
                output = open(outpath, 'w')
                output.write("#mega\n!Title: matrix.matrix;\n\n")
                s = ''
                for i in range(len(hash)):
                        s += '['+str(i+1)+'] ' +'#' + hash[i] + '\n'
                output.write(s+'\n')
                s = ''
                for i in range(len(hash)):
                        s += '\t'+str(i+1)
                output.write('[' + s+']\n')
                for i in range(len(dmatrix)):
                        s = []
                        for j in range(i):
                                s.append(str(dmatrix[i][j]))
                        string = '\t'.join(s)
                        output.write('['+str(i+1)+']\t'+string+'\n')
	def generate_matrix(records):
                result = []
                for i in range(len(records)):
                        result.append([])
                        for j in range(len(records)):
                                result[i].append(0)
                for i in range(len(records)):
                        for j in range (i+1, len(records)):
                                result[i][j] = comp(records[i], records[j])
                                result[j][i] = comp(records[i], records[j])
                return result
        matrix = read_matrix(path)
        dmatrix=generate_matrix(matrix)
        hash = []
        for record in matrix:
                hash.append(record[0])
        write_mega(hash, outpath)

def main():
	os.chdir(wdir)
	Matrix2Mega('all.tsv','all.meg')

main()
