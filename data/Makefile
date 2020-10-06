lc=python3
input=particles/ideal-no-rh-40p-1000ts-rinit-5velocity.txt
input=particles/ideal-rh-500ts.txt
# input=particles/ideal-rh-test-pre.txt
input=particles/ideal-rh-test.txt
# input=particles/ideal-rh-test.txt.1
input=particles/ideal-no-rh-test.txt
# input=particles/ideal-no-rh-test2.txt
input=particles/ideal-no-rh-test-check.txt


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
