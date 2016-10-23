CC = mpic++
CFLAGS = -W -Wall -std=c++11 -lboost_system -lboost_filesystem
BUILDDIR = build
DEPS = $(wildcard ./src/*.cpp)
OBJ = srcmain

srcmain: $(DEPS)
	mkdir -p $(BUILDDIR) && $(CC) -o $(BUILDDIR)/$(OBJ) $(DEPS) $(CFLAGS)

clean:
	rm -f $(BUILDDIR)/$(OBJ) && rm -r $(BUILDDIR)

