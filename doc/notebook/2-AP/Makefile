SRC= $(wildcard *.md)
IPYNB= $(SRC:.md=.ipynb)
RST= $(SRC:.md=.rst)

clean:
	rm *.ipynb
	rm *.rst 

ipynb: $(IPYNB) 

rst: $(RST)

%.ipynb: %.md
	ipymd --from markdown --to notebook $<

%.rst: %.ipynb
	runipy -o --skip-exception $<  >/dev/null 2> /dev/null
	ipython nbconvert $< --to rst 
	rm $<
