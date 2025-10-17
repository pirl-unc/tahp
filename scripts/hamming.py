import argparse
from pprint import pprint
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--counts-input', '-ci', required = True)
    parser.add_argument('--vertices-input', '-vi', required = True)
    parser.add_argument('--size-input', '-si', required = True)
    parser.add_argument('--output', '-el', required = True)
    parser.add_argument('--size-output', '-so', required = True)
    parser.add_argument('--vertex-cta-output', '-vco', required = True)
    parser.add_argument('--cta-counts-output', '-cco', required = True)
    return parser.parse_args()

def hamming_distance(string1, string2):
    if len(string1) != len(string2):
        raise ValueError("Strings must be of equal length")
    distance = 0
    for char1, char2 in zip(string1, string2):
        if char1 != char2:
            distance += 1
    return distance

def get_node_size(args):
    node_sizes = {}
    with open(args.size_input, 'r') as ifo:
        for line in ifo.readlines():
            line = line.strip()
            line = line.split('\t')
            print(line)
            size = line[0]
            node_id = line[1]
            node_sizes[node_id] = int(size)
    total_cells = sum(node_sizes.values())
    return node_sizes, total_cells

def parse_counts(args):
    cols_to_idx = {}
    cta_counts = {}
    total_cells = 0
    no_ctas_detected = 0
    with open(args.counts_input, 'r') as cfo:
        for line_idx, line in enumerate(cfo.readlines()):
            line = line.strip()
            line = line.replace('"', '')
            line = line.split('\t')
            if line_idx == 0:
                for col_idx, col in enumerate(line):
                    if col:
                        cols_to_idx[col] = col_idx
                        cta_counts[col] = 0
                    
            else:
                for k,v in cols_to_idx.items():
                    if line[v] == 'TRUE':
                        cta_counts[k] += 1
                if 'TRUE' not in line:
                    no_ctas_detected += 1  
                total_cells += 1
    for k,v in cta_counts.items():
        cta_counts[k] = int(v)/total_cells * 100

    cta_counts["None"] = no_ctas_detected/total_cells * 100

    pprint(cta_counts)

    with open(args.cta_counts_output, 'w') as ccfo:
        for k,v in cta_counts.items():
           ccfo.write("{}	{}\n".format(k, v)) 

    vertex_to_cta = {}

    with open(args.counts_input, 'r') as cfo:
        for line_idx, line in enumerate(cfo.readlines()):
            line = line.strip()
            line = line.replace('"', '')
            line = line.split('\t')
            if line_idx > 0:
                cta_list = []
                for k,v in cols_to_idx.items():
                        if line[v] == 'TRUE':
                            cta_list.append(k)
                cta_freq = 0
                print(cta_list)
                cta_source_set = ''
                for cta in cta_list:
                    if not(cta.startswith('ERV')) and not(cta.startswith('chr')) and 'C' not in cta_source_set:
                        cta_source_set += ('C')
                    elif cta.startswith('ERV') and 'E' not in cta_source_set:
                        cta_source_set += ('E')
                    elif cta.startswith('chr') and 'S' not in cta_source_set:
                        cta_source_set += ('S')
                if cta_source_set == '':
                    cta_source_set = 'None'
                vertex_code = ''.join(line[1:]).replace('TRUE', '1').replace('FALSE', '0')
                print("{}	{}".format(vertex_code, cta_source_set))
#                if vertex_code not in vertex_to_cta and cta_source_set != 'None':
                if cta_source_set != 'None':
                    print("unsorted: {}".format(cta_source_set))
                    print("sorted: {}".format(''.join(sorted(cta_source_set))))
                    vertex_to_cta[vertex_code] = ''.join(sorted(cta_source_set))
                else:
                    vertex_to_cta[vertex_code] = cta_source_set
    print("vertex_to_cta:")
    pprint(vertex_to_cta)

    return cta_counts, vertex_to_cta
              

def get_hammings(args, node_sizes, total_cells, vertex_to_cta, cta_counts):

    cta_colors = {
"None": "#FFFFFF",
"C": "#31688e",
"E": "#440154",
"S": "#35b779",
"CE": "#31688E",
"CS": "#339084",
"ES": "#3D5c67",
"CES": "#000000"}

    tups = []

    hammings = {}

    revised_node_sizes = {}

    revised_node_ctas = {}

    revised_node_colors = {}

    with open(args.vertices_input, 'r') as ifo:
        for line in ifo.readlines():
            line = line.strip()
            if line not in tups:
                tups.append(line)

    print("tups count: {}".format(len(tups)))

    tups2 = tups[:]

    for tup_idx, tup in enumerate(tups):
        print("tup: {}".format(tup))
        hammings[tup_idx] = []
        revised_node_sizes[tup_idx+1] = node_sizes[tup]
        revised_node_colors[tup_idx+1] = cta_colors[vertex_to_cta[tup]]
#        revised_node_colors[tup_idx+1] = "#000000"
        revised_node_ctas[tup_idx+1] = vertex_to_cta[tup]
        for tup2_idx, tup2 in enumerate(tups2):
           if hamming_distance(tup, tup2) == 1:
               hammings[tup_idx].append(tup2_idx)

    to_print = []

    for k,v in hammings.items():
       for recip in v:
               to_print.append("{}\t{}\n".format(k, recip))

    with open(args.output, 'w') as ofo:
        for entry in to_print:
            ofo.write("{}".format(entry))

    with open(args.size_output, 'w') as sofo:
#        sofo.write("1	0	0	NA	NA	NA\n")
        sum_node_size = sum(revised_node_sizes.values())
        source_cat_count = {}
        for k,v in revised_node_sizes.items():
            human_revised_node_ctas = revised_node_ctas[k].replace('C', 'CTA, ').replace('E', 'ERV, ').replace('S', 'SNV').strip(', ')    
            if human_revised_node_ctas not in source_cat_count.keys():
                source_cat_count[human_revised_node_ctas] = v
            else :
                source_cat_count[human_revised_node_ctas] += v

        for k,v in revised_node_sizes.items():
            print(k)
            print(revised_node_colors[k])
            print(revised_node_ctas[k])
            human_revised_node_ctas = revised_node_ctas[k].replace('C', 'CTA, ').replace('E', 'ERV, ').replace('S', 'SNV').strip(', ')
            sofo.write("{}	{}	{}	{}	{} ({:.2f}%)\n".format(k, v, v/total_cells*100, revised_node_colors[k], human_revised_node_ctas, float(source_cat_count[human_revised_node_ctas]/float(sum_node_size) * 100)))
#            sofo.write("{}	{}	{}	{}	{} ({:.2f}%)\n".format(k, v, v/total_cells*100, revised_node_colors[k], revised_node_ctas[k], float(cta_counts[revised_node_ctas[k]])))
#            sofo.write("{}	{}	{}	{} {}\n".format(k, v, v/total_cells*100, revised_node_colors[k], revised_node_ctas[k]))

    with open(args.vertex_cta_output, 'w') as vcofo:
        for k,v in revised_node_ctas.items():
            if v:
                vcofo.write("{}	{}\n".format(k, v))
            else:
                vcofo.write("{}	None\n".format(k))

def main():
  args = get_args()
  cta_counts, vertex_to_cta = parse_counts(args)
  node_sizes, total_cells = get_node_size(args)
  get_hammings(args, node_sizes, total_cells, vertex_to_cta, cta_counts)

if __name__=='__main__':
    main()
