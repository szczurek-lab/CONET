#Set this variable to directory containing CONET executable
bin_dir = './'
import sys
sys.path.append('../..')
import pandas as pd
import networkx as nx
import numpy as np
import conet 
import conet.src.data_converter.data_converter as dc
import conet.src.conet as c
import conet.src.conet_parameters as cp
import conet.src.inference_result as ir
from numpy import genfromtxt

def read_conet_to_digraph(path):
    with open(path) as f:
        content = f.readlines()
        
    content = [x.strip().replace('1_', '').replace(".000000", '') for x in content] 
    tree = nx.DiGraph()
    for x in content:
        parent = x.split('-')[0]
        child = x.split('-')[1]
        tree.add_edge(parent, child)
    return tree
    
def get_conet_attachment(path, nodes):
    with open(path) as f:
        content = f.readlines()
        
    content = [x.strip().replace('1_', '').replace(".000000", '') for x in content] 
    content = [(int(x.split(';')[1]), int(x.split(';')[2])) for x in content]
    content = [str(x) for x in content]
    content= [x.replace(' ', '') for x in content]
    
    node_to_cells = {}
    for node in nodes:
        node_to_cells[node] = []
    for cell in range(0, len(content)):
        node_to_cells[content[cell]].append(cell)
        
    return node_to_cells

def get_conet_inferred_counts(corrected_counts, conet_tree, attachment):
    post_order_node = list(reversed(list(nx.topological_sort(conet_tree))))
    conet_counts = np.zeros(corrected_counts.shape, dtype=np.float64)
    conet_counts.fill(-1)
    for i in range(0, len(post_order_node) - 1):
        node = post_order_node[i].replace('(', '').replace(')', '')
        start = int(node.split(',')[0])
        end = int(node.split(',')[1])
        cells = attachment[post_order_node[i]]
        sum_ = 0
        count = 0
        for c in cells:
            for j in range(start, end):
                if conet_counts[c,j] < 0:
                    count = count + 1
                    sum_ = sum_ + corrected_counts[c,j]
        
        if count > 0:
            inf_cn = round(sum_ /count)
            for c in cells:
                for j in range(start, end):
                    if conet_counts[c,j] < 0:
                        conet_counts[c,j] = inf_cn
                
        parent_node = list(conet_tree.predecessors(post_order_node[i]))[0]
        attachment[parent_node].extend(attachment[post_order_node[i]])
    
    for i in range(0, conet_counts.shape[0]):
        for j in range(0, conet_counts.shape[1]):
            if conet_counts[i,j] < 0:
                conet_counts[i,j] = 2
    return conet_counts       
    
    
def convert_counts_to_breakpoint_matrix(counts):
    counts = np.transpose(counts)
    brkps = np.copy(counts)
    brkps.fill(0)
    for i in range(0, counts.shape[0]):
        for j in range(0, counts.shape[1]):
            if (j == 0 and counts[i,j] != 2) or (counts[i,j] != counts[i, j -1]):
                brkps[i,j] = 1
    return brkps

def extract_events(cell, counts, result, result_index):
    i = 1
    while i < counts.shape[0]:
        count = counts[i, cell]
        if count != 2:
            start = i
            i = i +1
            while i < counts.shape[0] and counts[i, cell] == count:
                i = i + 1
            end = i - 1
            result[result_index, ] = [cell, count, start, end]
            result_index = result_index + 1
        else:
            i = i +1
    return result, result_index
  
  
from collections import Counter
def CONET_inference_pipeline( id_, file):
    corr_reads = genfromtxt('models/corrected_counts_' + id_, delimiter=';')
    cells = corr_reads.shape[0]
    loci = corr_reads.shape[1]
    print(loci)
    save_counts_in_CONET_format(bin_dir+"counts_synthetic", corr_reads, loci, cells)
    data_converter = dc.DataConverter(bin_dir+"counts_synthetic", 
                                  delimiter= ';', 
                                  default_bin_length = 1, 
                                  event_length_normalizer = loci,
                                  add_chromosome_ends = False,
                                  neutral_cn = 2.0)
    #breakpoint_candidates_indices = range(0, loci)
    breakpoint_candidates_indices = list(map(lambda x: int(x), list(genfromtxt('models/indices_' + id_, delimiter=';'))))
    data_converter.create_CoNET_input_files(breakpoint_candidates_indices, bin_dir, add_chr_ends_to_indices=False)

    conet = c.CONET(bin_dir + "CONET")
    params = cp.CONETParameters(data_size_prior_c = 0.05,data_dir = bin_dir, counts_penalty_c=100000, 
                                param_inf_iters=200000, seed = 2167, mixture_size=2, pt_inf_iters=800000,
                               use_event_lengths_in_attachment=False,
                               event_length_penalty_c = 1)
    conet.infer_tree(params)

    tree = read_conet_to_digraph(bin_dir+"inferred_tree")
    real_tree = nx.read_edgelist("models/tree_" + id_)


    conet_counts = get_conet_inferred_counts(corr_reads, tree, get_conet_attachment(bin_dir + "inferred_attachment", list(map(lambda x : str(x), list(tree.nodes)))))

    inferred_counts = np.transpose(conet_counts)
    real_counts = np.transpose(genfromtxt('models/counts_' + id_, delimiter=';'))
    MSE = ((inferred_counts - real_counts)**2).mean(axis=None)

    real_brkps = convert_counts_to_breakpoint_matrix(real_counts)

    inferred_brkps = convert_counts_to_breakpoint_matrix(inferred_counts)

    SD = (np.abs(real_brkps - inferred_brkps)).sum(axis=None) / inferred_counts.shape[1]

    FP = np.sum(np.logical_and(inferred_brkps == 1, real_brkps == 0)) / np.sum(inferred_brkps == 1)

    FN = np.sum(np.logical_and(inferred_brkps == 0, real_brkps == 1)) / np.sum(real_brkps == 1)

    result = np.zeros([1000000, 4])

    ind = 0
    for i in range(0, inferred_counts.shape[1]): 
         result, ind = extract_events(i, inferred_counts, result, ind)


    ind = ind-1
    EN = ind/inferred_counts.shape[1]

    result = result[range(0, ind), :]
    result_hash = []
    for i in range(0, ind):
        result_hash.append(str(result[i,1]) + "_" + str(result[i,2]) + "_" + str(result[i,3]))
    x = list(Counter(result_hash).values())
    SINGLETONS = x.count(1)


    real_edges = set(real_tree.edges)

    inferred_edges = set(tree.edges)


    intersection = real_edges.intersection(inferred_edges)
    EdgePrec =  len(intersection) / len(inferred_edges)
    EdgeSen = len(intersection) / len(real_edges)

    events = set(map(lambda x:x[1], inferred_edges))

    r_events = set(map(lambda x:x[1], real_edges))

    EventPrecision = len(r_events.intersection(events)) / len(events)

    EventFN = len(r_events.intersection(events)) / len(r_events)

    results = id_ + ";" + str(MSE) + ";"+ str(SD) + ";"+ str(FP) + ";"+ str(FN) + ";"+ str(EN) + ";"+ str(SINGLETONS)
    results = results + ";"+ str(len(inferred_edges)) + ";" + str(EdgePrec) + ";" + str(EdgeSen) + ";"
    results = results + str(EventPrecision) + ";" + str(EventFN) +"\n"
    print(results)
    file.write(results)
    

def save_counts_in_CONET_format(path, counts, no_loci, no_cells):
    counts = np.transpose(counts)
    add = np.zeros([no_loci, 4], dtype=np.float64)
    add.fill(1)
    add[:,1] = range(0, no_loci)
    add[:,2] = range(1, no_loci+1)
    full_counts = np.hstack([add, counts])
    names = np.zeros([1, 4 + no_cells], dtype=np.float64)
    names[0, 4:] = range(0, no_cells)
    full_counts = np.vstack([names, full_counts])
    np.savetxt(path, full_counts, delimiter=";")
    
thread = int(sys.argv[1])
THREADS = 8
print("Thread " + str(thread))
    
file = open("conet_results_" + str(thread), "a")
cells = [200, 1000]
tree_size = [20, 40]
for ts in tree_size:
    for cell in cells:
        for i in range(0,8):
            if i % THREADS == thread:
                loci = 1500
                if ts == 40:
                    loci = 10000 
                CONET_inference_pipeline("NEW_NOISE" + str(cell) +"_" + str(ts) + "_" + str(i), file)
                CONET_inference_pipeline("NEW_NOISE" + str(cell) +"_" + str(ts) + "_" + str(i), file)
                CONET_inference_pipeline("NO_NOISE" + str(cell) +"_" + str(ts) + "_" + str(i), file)
                CONET_inference_pipeline("NO_NOISE" + str(cell) +"_" + str(ts) + "_" + str(i), file)												
file.close()
