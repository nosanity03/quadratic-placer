build: qplacer.cpp solver.cpp solver.h
	# sample usage: `make build N=10`
ifdef N
	sed -i 's/if (n >= .*/if (n >= '$$N') return;/' qplacer.cpp
endif
	g++ -std=c++11 -o qplacer qplacer.cpp solver.cpp

show: out.txt outpad.txt ./visualization/visualize.py
	python ./visualization/visualize.py out.txt outpad.txt

clean:
	rm -f qplacer a.out out.txt outpad.txt
	cd doc; make clean

leakcheck: qplacer
	# sample usage: `make leakcheck input=<input filename>`
ifdef input
	valgrind --leak-check=yes ./qplacer $$input
endif

cachecheck: qplacer
	# sample usage: `make cachecheck input=<input filename>`
ifdef input
	valgrind --tool=cachegrind ./qplacer $$input
endif

doc: ./doc/docconfig ./doc/Makefile
	cd doc; make generate

.PHONY: doc clean
