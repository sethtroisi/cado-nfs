CXX=g++
HEADERS=/usr/local/include
build_tree=/tmp
PATH_LIBS=/usr/local/lib
LIBS=-lgmp
FLAGS=-Wall -Wextra -Werror -pedantic -Wpedantic -Wconversion -O3 -DNDEBUG

all:
	$(CXX) -I$(HEADERS) check_relations.cpp -o $(build_tree)/check_relations -L$(PATH_LIBS) -lntl $(LIBS) $(FLAGS)
