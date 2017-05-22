from sc_benchmark import run

TF_path  = "/media/tmo/data/work/datasets/TF/mm9_TFs.txt"
in_path  = "/media/tmo/data/work/datasets/megacell/out/matrix.subsets/3k/nr.rounds.250/from.exclude/subset.3k.from.exclude.txt/part-00000"
out_path = "/media/tmo/data/work/datasets/benchmarks/genie3/megacell/megacell.3k.genie3.txt"

nr_threads = 88

run(in_path, TF_path, out_path, nr_threads=nr_threads)