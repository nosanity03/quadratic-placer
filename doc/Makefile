generate: ../qplacer.cpp ../solver.cpp ../solver.h
	cp ../qplacer.cpp .
	cp ../solver.cpp .
	cp ../solver.h .
	doxygen docconfig
	google-chrome ./html/index.html
	cd latex/; make
	cp ./latex/refman.pdf reference_manual.pdf
	rm -f *.log

clean:
	rm -rf qplacer.cpp solver.h solver.cpp *.log
	cd latex/; make clean

new: 
	doxygen -g docconfig
