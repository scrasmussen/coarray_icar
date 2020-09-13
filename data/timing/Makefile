lc=python3.8
type=scaling
# type=speedup

# type=timing_particles
# file=_particles

all: plot

plot:
	$(lc) plot_$(type).py timing_results$(file).txt

cray:
	$(lc) plot_$(type).py cray_timing_results.txt -t "Strong Scaling, Cray"

new: gcc
gcc:
	$(lc) plot_$(type).py gcc9_timing_results.txt -t "Strong Scaling, New GCC 9"
