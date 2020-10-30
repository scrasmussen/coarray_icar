lc=python3
data=../../../icar_data/particles
files=data_no_rh.txt data_no_rh_all.txt


n=single
n=eight
files=$(n)-dry.txt $(data)/$(n)-dry.txt



all: run

run:
	$(lc) brunt_vaisala.py $(files)
