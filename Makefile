CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

BIN = mcmc_integral
SRC = src/mcmc_integral.cpp

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(SRC)

clean:
	rm -f $(BIN)
