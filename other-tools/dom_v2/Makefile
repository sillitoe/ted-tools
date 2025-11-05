# Makefile for dom_v2
# Builds the standalone domain binary from domain.c

# Toolchain
CC ?= cc

# Compile-time options
# -std=c99: domain.c is written in C99
# -D_XOPEN_SOURCE=600: enable POSIX/XSI extensions used by some libc functions
# -Wall -Wextra: useful warnings (kept reasonable for portability)
CFLAGS ?= -O2 -g -std=c99 -D_XOPEN_SOURCE=600 -Wall -Wextra -Wno-unused-parameter
CPPFLAGS ?=

# Linker options
LDFLAGS ?=
LDLIBS ?= -lm

# Targets and sources
BIN := domain
SRC := domain.c
HDR := domain.h
OBJ := $(SRC:.c=.o)

.PHONY: all clean debug sanitize

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.c $(HDR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(BIN)

# Build with debug info and without optimizations
# Usage: make clean && make debug
debug: CFLAGS := -O0 -g3 -std=c99 -D_XOPEN_SOURCE=600 -Wall -Wextra -DDEBUG=1
debug: LDLIBS := -lm
debug: $(BIN)

# Address/UBSan build (clang/gcc)
# Usage: make clean && make sanitize
sanitize: CFLAGS := -O1 -g -std=c99 -D_XOPEN_SOURCE=600 -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra
sanitize: LDFLAGS := -fsanitize=address,undefined
sanitize: LDLIBS := -lm
sanitize: $(BIN)
