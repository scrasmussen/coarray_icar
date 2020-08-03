lc=python3

input=ideal-rh-100p-500ts.txt
# input=ideal-no-rh-100p-500ts.txt
input=ideal-no-rh-100p-500ts-replacement.txt
input=ideal-no-rh-100p-500ts-both-replacement.txt


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
