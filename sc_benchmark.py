from GENIE3 import *
from numpy import *
from multiprocessing import cpu_count


def run(exp_path, tf_path, out, max_regulations=100000, nr_threads=cpu_count()):
    """
    :param exp_path: The path to the expression matrix.
    :param tf_path: The path to the file with the transcription factor list.
    :param out: The output path.
    :param max_regulations: The maximum number of regulatory connections to output.
    :param nr_threads: The nr of threads to use to run GENIE3, defaults to the nr of available processors.
    :return: void
    """
    print "running GENIE3 on", exp_path

    T_wall = time.time()

    # parse the gene names from the file
    exp_file = open(exp_path)
    exp_first_row = exp_file.readline().rstrip('\n').split('\t')
    nr_cols = len(exp_first_row)
    gene_names = [line.rstrip('\n').split('\t')[0] for line in exp_file.readlines()]
    exp_file.close()

    # load the matrix of expression data
    print "transposing the expression matrix"
    T_start = time.time()
    exp_data = loadtxt(exp_path, skiprows=1, usecols=range(1, nr_cols)).transpose()
    print "matrix.transpose() took:", time.time() - T_start

    # load the TF
    tf_file = open(tf_path)
    tr_factors = [line.rstrip('\n') for line in tf_file.readlines()]
    tf_file.close()

    # compute the result
    print "calling GENIE3 on the transposed matrix"
    (VIM, prediction_score, treeEstimators) = GENIE3(expr_data=exp_data,
                                                     gene_names=gene_names,
                                                     regulators=tr_factors,
                                                     nthreads=nr_threads)

    T_list = time.time()
    get_link_list(VIM,
                  gene_names=gene_names,
                  regulators=tr_factors,
                  file_name=out,
                  maxcount=max_regulations)
    print "calculating the link list:", time.time() - T_list

    print "GENIE3 wall time", exp_path, ":", time.time() - T_wall, "\t threads:", cpu_count()
