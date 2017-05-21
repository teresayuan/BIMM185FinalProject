f = raw_input('name of file or path: \n')
o = raw_input('name of output file name: \n')

lines = [line.rstrip('\r\n') for line in open(f)]

with open(o, 'w') as file:
    for i in range(0, len(lines)):
        if i % 4 == 0:
            file.write(lines[i][:-2] + '/' + lines[i][-1] + '\n')
        else:
            file.write(lines[i] + '\n')