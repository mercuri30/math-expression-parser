all:
	g++ -Wall -Wextra -Wconversion unit.cpp expression.cpp -lcppunit -ldl

clean:
	rm -f a.out
