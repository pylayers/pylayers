SRC= $(wildcard *.ipynb)
OBJ= $(SRC:.ipynb=.ipyex)

all: $(OBJ)

%.ipyex: %.ipynb
	runipy --pylab $< $@ 
