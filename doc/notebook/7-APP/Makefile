SRC= $(wildcard *.ipynb)
OBJ= $(SRC:.ipynb=.rst)

all: $(OBJ)

%.rst: %.ipynb
	echo "Running" $<
	runipy --pylab $<  2> /dev/null
	ipython nbconvert $< --to rst 
