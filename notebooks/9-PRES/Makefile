SRC= $(wildcard *.ipynb)
OBJ= $(SRC:.ipynb=.html)

all: $(OBJ)

%.html: %.ipynb
	ipython nbconvert $< $@ 
