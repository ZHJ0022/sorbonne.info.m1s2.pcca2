CC = gcc
CFLAGS = -Wall -Wextra -g
INCLUDES = -Iinclude

TEST ?= test1

TARGET = build/$(TEST)

SRC = src/sturmNaif.c src/util.c tests/$(TEST).c

$(TARGET): $(SRC)
	mkdir -p build
	$(CC) $(CFLAGS) $(INCLUDES) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f build/*
