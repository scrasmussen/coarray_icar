lc=python3.8
type=strong_scaling
# type=speedup

# type=timing_particles
# file=_particles

all: plot

plot:
	$(lc) plot_$(type).py timing_results$(file).txt

c: cray
cray:
	$(lc) plot_$(type).py cray_strong_scaling.txt -t "Strong Scaling Cray"
	# $(lc) plot_$(type).py cray_timing_results.txt -t "Strong Scaling, Cray"


new: gcc
gcc:
	$(lc) plot_$(type).py gcc9_timing_results.txt -t "Strong Scaling, New GCC 9"


s: strong
strong:
	$(lc) plot_strong_scaling.py cray_strong_scaling.txt -t "Strong Scaling Cray"

w: weak
weak:
	$(lc) plot_weak_scaling.py cray_weak_scaling.txt -t "Weak Scaling Cray"
