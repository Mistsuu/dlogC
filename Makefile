# You can change this to release|memcheck|verbose|static
BUILD = release

# Settings for folders
OBJECT_DIR  = obj
LIBRARY_DIR = src
SOURCE_DIR  = src

# Compilers options
CC              = gcc
CFLAGS.release  = -I$(LIBRARY_DIR) -O2
CFLAGS.memcheck = -I$(LIBRARY_DIR) -O2 -ggdb -fsanitize=address
CFLAGS.verbose  = -I$(LIBRARY_DIR) -O2 -ggdb -D DLOG_VERBOSE=1
CFLAGS.static   = -I$(LIBRARY_DIR) -O2 -static
CFLAGS          = $(CFLAGS.$(BUILD))
LIBS            = -lgmp -lpthread

# Requirements and stuffs
FULLDEPS := $(shell find $(LIBRARY_DIR) -name '*.h')
FULLOBJS := $(shell find $(SOURCE_DIR) -name '*.c' | sed -e "s/^$(SOURCE_DIR)/$(OBJECT_DIR)/" | sed -e "s/\\.c$$/.o/")

dlog: $(FULLOBJS) $(OBJECT_DIR)/main.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(OBJECT_DIR):
	mkdir -p $(OBJECT_DIR)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.c $(FULLDEPS) | $(OBJECT_DIR)
	$(CC) -c -o $@ $< $(CFLAGS)

# For commands
.PHONY: clean run dbg
clean:
	rm -f $(SOURCE_DIR)/*.o $(OBJECT_DIR)/*.o ./dlog
