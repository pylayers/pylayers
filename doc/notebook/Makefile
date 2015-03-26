all: 
	cd 1-GIS && make rst  
	cd 2-AP && make rst  
	cd 3-PHY && make rst 
	cd 4-MOB && make rst  
	cd 5-SIM && make rst 
	cd 7-APP && make rst
	cd 8-MISC && make rst  
	cd 9-PRES && make rst

conv:
	for i in `find . -name "*.ipynb"`;\
	do \
		a=`basename $$i`;\
		b=`dirname $$i`;\
		ipymd --from notebook --to markdown --overwrite $$b/$$a;\
	done	

