from sc_benchmark import run

TF_path = "/media/tmo/data/work/datasets/TF/mm9_TFs.txt"
in_path = "/media/tmo/data/work/datasets/zeisel/expression_sara_filtered.txt"
out_path = "/media/tmo/data/work/datasets/benchmarks/genie3/zeisel/zeisel.filtered.genie3.txt"

run(in_path, TF_path, out_path)
