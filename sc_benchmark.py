from GENIE3 import *
from numpy import *
from multiprocessing import cpu_count


def run(exp_path, tf_path, out, max_regulations=100000):
    print "running GENIE3  on", exp_path

    start_time = time.clock()

    # parse the gene names from the file
    exp_file = open(exp_path)
    exp_first_row = exp_file.readline().rstrip('\n').split('\t')
    nr_cols = len(exp_first_row)
    gene_names = [line.rstrip('\n').split('\t')[0] for line in exp_file.readlines()]
    exp_file.close()

    # load the matrix of expression data
    exp_data = loadtxt(exp_path, skiprows=1, usecols=range(1, nr_cols)).transpose()

    # load the TF
    tf_file = open(tf_path)
    tr_factors = [line.rstrip('\n') for line in tf_file.readlines()]
    tf_file.close()

    # compute the result
    (VIM, prediction_score, treeEstimators) = GENIE3(expr_data=exp_data,
                                                     gene_names=gene_names,
                                                     regulators=tr_factors,
                                                     nthreads=cpu_count())
    get_link_list(VIM,
                  gene_names=gene_names,
                  regulators=tr_factors,
                  file_name=out,
                  maxcount=max_regulations)

    print "elapsed time for", exp_path, ":", start_time - time.clock(), "\t threads:", cpu_count()
