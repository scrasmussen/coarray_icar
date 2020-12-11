lc=python3.8
type=strong_scaling
# type=speedup

type=timing_particles
file=_particles
machine=cray
# machine=cheyenne

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

h: hist
hist:
	$(lc) hist_wind_comm.py ../../../icar_data/particles

ha: histall
histall:
	$(lc) hist_all_wind_comm.py ../../../icar_data/particles/wind

m: multinode
multinode:
	$(lc) plot_multinode.py

p: particles
particles:
	$(lc) plot_timing_particles.py cray_particle_timing.txt -t "Particle Scaling"

se: strongeff
strongeff:
	$(lc) plot_strong_efficiency.py

ss: strongspeed
strongspeed:
	$(lc) plot_strong_speedup.py

ws: weakspeed
weakspeed:
	$(lc) plot_weak_speedup.py

s: strong
strong:
	$(lc) plot_strong_scaling.py -t "Strong Scaling"

w: weak
weak:
	$(lc) plot_weak_scaling.py  -t "Weak Scaling"

wind:
	$(lc) plot_wind_scaling.py

wc:
	$(lc) plot_wind_comm.py

we:
	$(lc) plot_weak_efficiency.py

wn:
	$(lc) plot_wind_node_scaling.py
