from sc_benchmark import run

TF_path = "/media/tmo/data/work/datasets/TF/mm9_TFs.txt"
in_path = "/media/tmo/data/work/datasets/macosko/in/allEsetMR.tsv"
out_path = "/media/tmo/data/work/datasets/benchmarks/genie3/macosko/macosko.full.genie3.txt"

run(in_path, TF_path, out_path)
