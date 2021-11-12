.DEFAULT_GOAL := compile

compile:
	g++ -o ./bin/main -I ./include ./src/*.cpp
run:
	./bin/main