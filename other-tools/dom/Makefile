CC=gcc
# Use a legacy C dialect for K&R-style sources and relax strict checks on modern compilers (clang on macOS)
# -std=gnu89: allow implicit-int and old-style (K&R) function definitions
# -fcommon: avoid multiple-definition errors with tentative definitions on newer toolchains
# -Wno-implicit-int/-Wno-implicit-function-declaration: suppress errors promoted by newer standards
# -Wno-deprecated-non-prototype: silence warnings about K&R prototypes
# -Wno-return-type: some functions in legacy code omit explicit return types
CFLAGS=-I. -O -std=gnu89 -fcommon -Wno-implicit-int -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -Wno-return-type

dom: dom.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o
	$(CC) -o dom dom.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o -lm

clean:
	/bin/rm -f dom dom.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o
