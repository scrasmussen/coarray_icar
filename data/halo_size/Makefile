lc=python3.8
type=halo_times

input=cray_timing_results.txt
# input=cray_timing_results_optimized.txt

all: cray

cray:
	$(lc) plot_$(type).py $(input) -t "Halo size communication time, 500x500x30 domain, Cray"
