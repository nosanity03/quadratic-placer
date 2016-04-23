build: test.cpp solver.cpp solver.h
	g++ -std=c++11 -o test test.cpp solver.cpp

show: out.txt ../visualization/visualize.py
	python ../visualization/visualize.py out.txt

clean:
	rm -f test a.out out.txt
