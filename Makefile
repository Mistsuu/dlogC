# Settings for folders
OBJECT_DIR  = obj
LIBRARY_DIR = src
SOURCE_DIR  = src

# Compilers options
CC          = gcc
CFLAGS      = -I$(LIBRARY_DIR) -O3 -ggdb -fsanitize=address
LIBS        = -lgmp

# Requirements and stuffs
FULLDEPS := $(shell find $(LIBRARY_DIR) -name '*.h')
FULLOBJS := $(shell find $(SOURCE_DIR) -name '*.c' | sed -e "s/^$(SOURCE_DIR)/$(OBJECT_DIR)/" | sed -e "s/\\.c$$/.o/")

$(OBJECT_DIR):
	mkdir -p $(OBJECT_DIR)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.c $(FULLDEPS) | $(OBJECT_DIR)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(FULLOBJS) $(OBJECT_DIR)/main.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# For commands
.PHONY: clean run dbg
clean:
	rm -f $(SOURCE_DIR)/*.o $(OBJECT_DIR)/*.o ./main

run: main
	./main

dbg: main
	gdb ./main