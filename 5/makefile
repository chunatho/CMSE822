
CXX := g++

CPPFLAGS := -std=c++11 -O3 -fopenmp

CC := mpicc

CFLAGS := -std=gnu11 -O3 -fopenmp

SUBSET_SRCDIR := src/subset
SUBSET_CPPFILES := $(filter-out ./$(SUBSET_SRCDIR)/main.cpp,$(wildcard ./$(SUBSET_SRCDIR)/*.cpp))
SUBSET_OFILES := $(SUBSET_CPPFILES:.cpp=.o)

CODEBRK_SRCDIR := src/codebreaker


all: subset encrypter codebreaker
	@echo "Compilation complete."


subset: $(SUBSET_OFILES)
	$(CXX) -o subset.x ./$(SUBSET_SRCDIR)/main.cpp $(SUBSET_OFILES) $(CPPFLAGS)

$(subset_SRCDIR)/%.o: $(SUBSET_SRCDIR)/%.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)


encrypter:
	$(CC) -o encrypter.x ./$(CODEBRK_SRCDIR)/encrypter.c $(CFLAGS)

codebreaker:
	$(CC) -o codebreaker.x ./$(CODEBRK_SRCDIR)/codebreaker.c $(CFLAGS)


clean:
	-rm *.x $(SUBSET_OFILES) $(CS_OFILES)
