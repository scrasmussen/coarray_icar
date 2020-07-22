lc=python3

input=ideal-rh-100p-500ts.txt
input=ideal-no-rh-100p-500ts.txt


all: build

build:
	$(lc) plot.py $(input)

t: temp
temp:
	$(lc) temp_plot.py cafimage/$(input)
