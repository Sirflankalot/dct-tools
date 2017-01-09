CXX := g++

ARGS := -Ofast -Wall -Wextra -Wpedantic -std=c++1z

SOURCES := $(wildcard *.cpp)
ASM     := $(patsubst %.cpp,assembly/%.s,$(SOURCES))
ASMSSE  := $(patsubst %.cpp,assembly/%-sse.s,$(SOURCES))
ASMAVX  := $(patsubst %.cpp,assembly/%-avx.s,$(SOURCES))
OBJ     := $(patsubst %.cpp,%.o,$(SOURCES))
OBJSSE  := $(patsubst %.cpp,%-sse.o,$(SOURCES))
OBJAVX  := $(patsubst %.cpp,%-avx.o,$(SOURCES))

.PHONY: all clean asm execs

all: execs

execs-simd: ARGS += -DSIMD
execs-simd: execs

asm-simd: ARGS += -DSIMD
asm-simd: asm

execs: dct dct-sse dct-avx

asm: $(ASM) $(ASMSSE) $(ASMAVX)

assembly:
	mkdir -p $@

assembly/%.s: %.cpp | assembly
	$(CXX) $(ARGS) $^ -S -masm=intel -o $@

assembly/%-sse.s: %.cpp | assembly
	$(CXX) $(ARGS) $^ -S -masm=intel -o $@

assembly/%-avx.s: %.cpp | assembly
	$(CXX) $(ARGS) $^ -S -masm=intel -o $@

%.o: %.cpp
	$(CXX) $(ARGS) $^ -c -o $@

%-sse.o: %.cpp
	$(CXX) $(ARGS) $^ -c -o $@

%-avx.o: %.cpp
	$(CXX) $(ARGS) $^ -c -o $@

dct: $(OBJ)
dct-sse: ARGS += -msse4.2
dct-sse: $(OBJSSE)
dct-avx: ARGS += -mavx -mavx2 -mfma
dct-avx: $(OBJAVX)

dct dct-sse dct-avx:
	$(CXX) $(ARGS) $^ -o $@

clean:
	rm -rf $(OBJ) $(OBJSSE) $(OBJAVX) $(ASM) $(ASMSSE) $(ASMAVX) assembly dct
