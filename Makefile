all : enigma

enigma : enigma.cc
	g++ -Wall -O3 -mssse3 -o enigma enigma.cc
