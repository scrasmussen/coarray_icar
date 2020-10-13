lc=python3

# ----- one particle runs for testing -----
input=particles/single-200ts-no-sounding-dry.txt
input=particles/single-200ts-sounding-dry.txt
# input=particles/single-200ts-no-sounding-wet.txt
# input=particles/single-200ts-no-sounding-wet-16p.txt
input=particles/single-sounding-wet.txt
input=particles/single-sounding-dry.txt
input=particles/eight-wet-dry.txt
input=particles/eight-dry.txt


all: build

build:
	$(lc) plot.py $(input)

b: brunt
brunt:
	$(lc) brunt_vaisala.py $(input)

f: fft
fft:
	$(lc) fft.py $(input)

g: gif
gif:
	$(lc) plot.py $(input) gif

p: pub
pub:
	$(lc) pub.py $(input)

p2: plot_2d
plot_2d:
	$(lc) plot_2d.py $(input)



.PHONY: s sounding
s: sounding
sounding:
	$(lc) get_sounding.py

t: temp
temp:
	$(lc) temperature_plot.py $(input)

a: altitude
altitude:
	$(lc) altitude_plot.py $(input)

r: random

random:
	$(lc) random_plot.py $(input)
