# Executable pipeline
import sys
import os

if len(sys.argv) != 8 or sys.argv[1] == '-h':
    if sys.argv[1] != '-h':
        print 'wrong number of files'

    print '185.py reference sample1forward sample1reverse sample2forward sample2reverse -o outputdir'

    sys.exit()

outputdir = sys.argv[7]
os.system('mkdir ' + outputdir)

with open('./'+outputdir+'/info.txt', 'w') as file:
    file.write('reference file: ' + sys.argv[1])
    file.write('\nsample 1 forward file: ' + sys.argv[2])
    file.write('\nsample 1 reverse file: ' + sys.argv[3])
    file.write('\nsample 2 forward file: ' + sys.argv[4])
    file.write('\nsample 2 reverse file: ' + sys.argv[5])

for i in [1, 2, 3, 4, 5]:
    os.system('cp ' + sys.argv[i] + ' ' + outputdir)

os.system('mv ' + ' ./' + outputdir + '/' + sys.argv[1] + ' ./' + outputdir + '/reference.fna')
os.system('mv ' + ' ./' + outputdir + '/' + sys.argv[2] + ' ./' + outputdir + '/sample1_forward.fastq')
os.system('mv ' + ' ./' + outputdir + '/' + sys.argv[3] + ' ./' + outputdir + '/sample1_reverse.fastq')
os.system('mv ' + ' ./' + outputdir + '/' + sys.argv[4] + ' ./' + outputdir + '/sample2_forward.fastq')
os.system('mv ' + ' ./' + outputdir + '/' + sys.argv[5] + ' ./' + outputdir + '/sample2_reverse.fastq')


