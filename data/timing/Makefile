lc=python3.8
type=scaling
# type=speedup

# type=timing_particles
# file=_particles

all: plot

plot:
	$(lc) plot_$(type).py timing_results$(file).txt

cray:
	$(lc) plot_$(type).py cray_timing_results_optimized.txt -t "Strong Scaling Cray"
	# $(lc) plot_$(type).py cray_timing_results.txt -t "Strong Scaling, Cray"

weak:
	$(lc) plot_weak_scaling.py cray_weak_scaling.txt -t "Weak Scaling Cray"

new: gcc
gcc:
	$(lc) plot_$(type).py gcc9_timing_results.txt -t "Strong Scaling, New GCC 9"
