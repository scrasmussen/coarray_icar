lc=python3.8
type=halo_times

input=cray_timing_results.txt
input=cray_timing_results_optimized.txt

all: cray

cray:
	$(lc) plot_$(type).py $(input) -t "Halo depth communication time, Cray"

no:
	$(lc) plot_$(type).py $(input)
