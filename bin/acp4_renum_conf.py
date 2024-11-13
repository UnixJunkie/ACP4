#!/usr/bin/env python3

# append a conformer number (like _1) to each molecule in a .sdf file
# 1st conformer is assigned _0

import sys

sdf_input_fn = sys.argv[1]
seen_names = {}
# 2nd line in a .sdf file is supposed to be a molecule name
mol_name_line_num = 1

with open(sdf_input_fn) as input:
    for i, line in enumerate(input.readlines()):
        if i == mol_name_line_num:
            mol_name = line.strip()
            prev_count = seen_names.get(mol_name, 0)
            print('%s_%d' % (mol_name, prev_count))
            seen_names[mol_name] = prev_count + 1
            mol_name_flag = False
        else:
            print(line, end='')
            if line.strip() == "$$$$":
                # two lines down will be a molecule name
                mol_name_line_num = i + 2
