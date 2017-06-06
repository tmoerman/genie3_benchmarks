from sc_benchmark import run

TF_path  = "/media/tmo/data/work/datasets/katina/threeLines_TFs_10x.lst"
in_path  = "/media/tmo/data/work/datasets/katina/threeLines_10x.tsv"
out_path = "/media/tmo/data/work/datasets/katina/out/threeLines_10x.out.tsv"

nr_threads = 16

run(in_path, TF_path, out_path, nr_threads=nr_threads)
