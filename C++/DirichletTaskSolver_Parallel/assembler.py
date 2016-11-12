import os
import sys

if len(sys.argv) != 2:
    raise Exception("set the only parameter path to info file");

num_points = -1
num_procs = -1
num_row_procs = -1
with open(sys.argv[1]) as fin:
    for i, line in enumerate(fin):
        if (i == 1):
            num_procs = int(line)

        if (i == 9):
            num_row_procs = int(line)

        if (i == 11):
            num_points = int(line)

for t in ['VALUE', 'TRUE']:
    values = []
    for _ in xrange(num_points):
        values.append([])
        for _ in xrange(num_points):
            values[-1].append(0.0)

    cur_row_index = 0
    cur_col_index = 0
    for proc in xrange(1, num_procs + 1):
        with open('{0}_PART_POINTS_{1}_PROC_{2}'.format(t, num_points, proc), 'r') as fin:
            temp = 0
            for line in fin:
                vals = line.split(',')[: -1]
                temp = len(vals)
                for i, v in enumerate(vals):
                    values[cur_row_index][cur_col_index + i] = v
                cur_row_index += 1

            if proc % num_row_procs == 0:
                cur_col_index += temp
                cur_row_index = 0

    with open('FINAL_{0}_POINTS_{1}_PROC_{2}'.format(t, num_points, num_procs), 'w') as fout:
        for row in values:
            for c in row:
                fout.write('{},'.format(c))
            fout.write('\n')
