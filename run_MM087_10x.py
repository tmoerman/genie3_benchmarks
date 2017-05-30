from sc_benchmark import run

TF_path  = "/media/tmo/data1/work/datasets/katina/MM087_TFs_10x.lst"
in_path  = "/media/tmo/data1/work/datasets/katina/MM087_10x.tsv"
out_path = "/media/tmo/data1/work/datasets/katina/out/MM087_10x.out.tsv"

nr_threads = 44

run(in_path, TF_path, out_path, nr_threads=nr_threads)