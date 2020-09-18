lc=python3

input=ideal-rh-100p-500ts.txt
# input=ideal-no-rh-100p-500ts.txt
input=particles/old/ideal-no-rh-100p-500ts-replacement.txt
# input=ideal-no-rh-100p-500ts-both-replacement.txt
# input=cray-ideal-no-rh-100p-500ts.txt
# input=ideal-no-rh-100p-200ts.txt
# input=ideal-no-rh-10p-300ts.txt
input=particles/ideal-no-rh-40p-400ts.txt
input=particles/ideal-no-rh-40p-500ts.txt
input=particles/ideal-no-rh-40p-500ts-rinit.txt
input=particles/ideal-no-rh-40p-1000ts-rinit-5velocity.txt

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
