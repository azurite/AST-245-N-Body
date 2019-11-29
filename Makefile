# https://www.gnu.org/software/make/manual/make.html

# $@ means everything to the left of the ":"
# $^ means everything to the right of the ":"
# $< means first item to the right of the ":"
# $(patsubst pattern, replacement, text)

IDIR = include
SDIR = src
ODIR = obj
LDIR = lib

CXX = g++
CXXFLAGS = -I /usr/include/eigen3 -I $(IDIR)

_DEPS = data.hpp particle.hpp
DEPS = $(patsubst %, $(IDIR)/%, $(_DEPS))

_OBJ = main.o data.o particle.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

main: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@.out $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *.out
