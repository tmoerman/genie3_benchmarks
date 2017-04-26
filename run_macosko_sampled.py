from sc_benchmark import run

TF_path = "/media/tmo/data/work/datasets/TF/mm9_TFs.txt"
in_path = "/media/tmo/data/work/datasets/macosko/in/sampledEsetMR.tsv"
out_path = "/media/tmo/data/work/datasets/benchmarks/genie3/macosko/macosko.sampled.genie3.txt"

# nr of threads for Nostromo
nr_threads = 44

run(in_path, TF_path, out_path, nr_threads=nr_threads)
