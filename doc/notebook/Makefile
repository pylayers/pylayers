all: 
	cd 1-GIS && make all
	cd 2-AP && make all
	cd 3-PHY && make all
	cd 4-MOB && make all
	cd 5-SIM && make all
	cd 6-NET && make all
	cd 7-APP && make all
	cd 8-MISC && make all
	cd 9-PRES && make all
web:	
	for i in `find . -name "*.ipyex"`;do  a=`basename $$i`;b=`dirname $$i`; cp $$i $$NOTEWEB/$$b/$${a%%.ipyex}.ipynb ; done
	cp TOC.ipynb $$NOTEWEB

