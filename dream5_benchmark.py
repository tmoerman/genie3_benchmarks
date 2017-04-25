from GENIE3 import *
from numpy import *
from dreamtools import D5C4

# wd = "/Users/tmo/Work/ghb2016/data/Dream_challenge/GNI/training data/"
wd = "/media/tmo/data/work/datasets/dream5/training data/"

nw1 = "Network 1 - in silico"
nw2 = "Network 2 - S. aureus"
nw3 = "Network 3 - E. coli"
nw4 = "Network 4 - S. cerevisiae"

nw1_exp = wd + nw1 + "/net1_expression_data.tsv"
nw1_TFs = wd + nw1 + "/net1_transcription_factors.tsv"
nw2_exp = wd + nw2 + "/net2_expression_data.tsv"
nw2_TFs = wd + nw2 + "/net2_transcription_factors.tsv"
nw3_exp = wd + nw3 + "/net3_expression_data.tsv"
nw3_TFs = wd + nw3 + "/net3_transcription_factors.tsv"
nw4_exp = wd + nw4 + "/net4_expression_data.tsv"
nw4_TFs = wd + nw4 + "/net4_transcription_factors.tsv"

out_path = "/out.genie3"
suffix   = ".txt"


def run(exp_path, tf_path, out, k='sqrt', max_regulations=100000):
    #
    exp_file = open(exp_path)
    gene_names = exp_file.readline().rstrip('\n').split('\t')
    exp_file.close()
    #
    tf_file = open(tf_path)
    tr_factors = [line.rstrip('\n') for line in tf_file.readlines()]
    tf_file.close()
    #
    exp_data = loadtxt(exp_path, skiprows=1)
    (VIM, prediction_score, treeEstimators) = GENIE3(expr_data=exp_data,
                                                     gene_names=gene_names,
                                                     regulators=tr_factors,
                                                     K=k,
                                                     nthreads=88)
    get_link_list(VIM,
                  gene_names=gene_names,
                  regulators=tr_factors,
                  file_name=out,
                  maxcount=max_regulations)

def run_nw1():
    run(nw1_exp, nw1_TFs, wd + nw1 + out_path + suffix)

def run_nw2():
    run(nw2_exp, nw2_TFs, wd + nw2 + out_path + suffix)

def run_nw3():
    run(nw3_exp, nw3_TFs, wd + nw3 + out_path + suffix)

def run_nw4():
    run(nw4_exp, nw4_TFs, wd + nw4 + out_path + suffix)

def eval_template():
    s = D5C4()
    filenames = s.download_template()
    return s.score(filenames)

def eval_raw():
    s = D5C4()
    out1 = wd + nw1 + out_path + suffix
    out3 = wd + nw3 + out_path + suffix
    out4 = wd + nw4 + out_path + suffix
    filenames = [out1, out3, out4]
    return s.score(filenames)

def eval_100k():
    s = D5C4()
    out1_100k = wd + nw1 + out_path + ".100k" + suffix
    out3_100k = wd + nw3 + out_path + ".100k" + suffix
    out4_100k = wd + nw4 + out_path + ".100k" + suffix
    filenames = [out1_100k, out3_100k, out4_100k]
    return s.score(filenames)