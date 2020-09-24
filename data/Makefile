lc=python3

input=particles/ideal-no-rh-40p-1000ts-rinit-5velocity.txt
input=particles/ideal-rh-500ts.txt
input=particles/ideal-rh-test.txt
# input=particles/ideal-rh-test-pre.txt
input=particles/ideal-no-rh-test.txt

all: build

build:
	$(lc) plot.py $(input)

g: gif
gif:
	$(lc) plot.py $(input) gif

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
