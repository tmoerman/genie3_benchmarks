from GENIE3 import *
from numpy import *
from multiprocessing import cpu_count


def run(exp_path, tf_path, out, k='sqrt', max_regulations=100000):
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
                                                     K=k,
                                                     nthreads=cpu_count())
    get_link_list(VIM,
                  gene_names=gene_names,
                  regulators=tr_factors,
                  file_name=out,
                  maxcount=max_regulations)

    print "elapsed time for", exp_path, ":", start_time - time.clock(), "\tthreads:", cpu_count(), "\tK:", k


TF_path = "/media/tmo/data/work/datasets/TF/mm9_TFs.txt"


def run_zeisel():
    in_path = "/media/tmo/data/work/datasets/zeisel/expression_sara_filtered.txt"
    out_dir = "/media/tmo/data/work/datasets/benchmarks/genie3/zeisel/"
    run(in_path, TF_path, out_dir + "zeisel.filtered.genie3.k.sqrt.txt", k='sqrt')
    run(in_path, TF_path, out_dir + "zeisel.filtered.genie3.k.all.txt",  k='all')


def run_macosko_sampled():
    in_path = "/media/tmo/data/work/datasets/macosko/in/sampledEsetMR.tsv"
    out_dir = "/media/tmo/data/work/datasets/benchmarks/genie3/macosko/"
    run(in_path, TF_path, out_dir + "macosko.sampled.genie3.k.sqrt.txt", k='sqrt')
    run(in_path, TF_path, out_dir + "macosko.sampled.genie3.k.all.txt", k='all')


def run_macosko_full():
    in_path = "/media/tmo/data/work/datasets/macosko/in/allEsetMR.tsv"
    out_dir = "/media/tmo/data/work/datasets/benchmarks/genie3/macosko/"
    run(in_path, TF_path, out_dir + "macosko.full.genie3.k.sqrt.txt", k='sqrt')
    run(in_path, TF_path, out_dir + "macosko.full.genie3.k.all.txt", k='all')