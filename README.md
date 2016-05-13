# Quadratic-Placer

A quadratic placer for placements of chips with less than 100000 gates. Can be 
extended to a higher number of gates as well. The ASIC has been assumed to have
x-coordinate from 0 to 100, and y-coordinate from 0-100.

# Usage

## Build

Build the project using the Makefile provided. An optional input is the value
of N. It decides the number of recursions that the placer need to take. If N =
10, it means that the gates need to be distributed into 2^10 basic cells of
equal area.

`make build` or `make build N=10`

## Run

The executable needs an input file in a specific format. Some input files are
presented in the folder named "benchmarks". Use one of these input files to run.

`./qplacer ./benchmarks/industry1`

## Visualize

The final placement is written back to an output file. The output can be
visualized using a python script. Use the Makefile for the same.

`make show`

# Dependencies

## Build

* **C++11 std library**: Use `-std=c++11` flag while using g++.

## Visualize

* **Python 2**: Python 2 is mostly preinstalled in linux systems.
* **pyplot**: pyplot is present in the python library matplotlib. Install 
`python-matplotlib` in your machine.
