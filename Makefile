CC = gcc
CFLAGS = -Wall -Wextra -g
INCLUDES = -Iinclude
LDLIBS = -lgmp -lm

BUILD_DIR = build
SRC=src/sturmNaif.c \
    src/util.c \
    src/poly.c \
    src/calculate.c

TEST1 = $(BUILD_DIR)/test1
TEST2 = $(BUILD_DIR)/test2
TEST3 = $(BUILD_DIR)/test3

all: $(TEST1) $(TEST2) $(TEST3)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TEST1): $(SRC) tests/test1.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LDLIBS)

$(TEST2): $(SRC) tests/test2.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LDLIBS)

$(TEST3): $(SRC) tests/test3.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LDLIBS)

run_t1: $(TEST1)
	@echo "Running test1..."
	@./$(TEST1)

run_t2: $(TEST2)
	@echo "Running test2..."
	@./$(TEST2)

run_t3: $(TEST3)
	@echo "Running test3..."
	@./$(TEST3)


run: run_t1 run_t2 run_t3
	@echo "All tests completed!"

clean:
	rm -f $(BUILD_DIR)/*