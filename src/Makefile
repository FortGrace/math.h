CC=gcc
CFLAGS= -Wall -Werror -Wextra
PLATFORM := $(shell uname)

all: test gcov_report

ifeq ($(PLATFORM), Darwin)
test: s21_math.a
	$(CC) --coverage s21_test.c s21_math.c -o test s21_math.a -lcheck
	./test
else
test: s21_math.a
	$(CC) --coverage s21_test.c s21_math.c -o test s21_math.a -lcheck_pic -lpthread -lrt -lm -lsubunit
	./test
endif
gcov_report: test
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info

s21_math.o: s21_math.c s21_math.h
	$(CC) $(CFLAGS) -c s21_math.c

s21_math.a: s21_math.o
	ar rc s21_math.a s21_math.o
	ranlib s21_math.a

clean:
	rm -rf ./*.o ./*.gcda ./*.gcno ./*.info ./*.a test
	rm -rf report

rebuild: clean all
