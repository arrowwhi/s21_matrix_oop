CC=g++
SRC=s21_matrix.cc
OBJ=s21_matrix.o
CFLAGS= -g -Wall -Werror -Wextra -std=c++17
TESTFLAGS=-lgtest

all: gcov_report

my_test:
	$(CC) $(SRC) main.cc -o test.exe
	./test.exe
	rm test.exe

s21_matrix_oop.a:
	$(CC) $(CFLAGS) $(SRC) -c
	ar -rcs s21_matrix_oop.a $(OBJ)

test: s21_matrix_oop.a
	$(CC) $(CFLAGS) s21_matrix_test.cc s21_matrix_oop.a -o test.out $(TESTFLAGS)
	./test.out

gcov_report:
	$(CC) s21_matrix_test.cc -c
	$(CC) --coverage  s21_matrix.cc  s21_matrix_test.o -o test.out $(TESTFLAGS)
	./test.out
	lcov -t "test" -o test.info -c -d ./
	genhtml -o report test.info
	open report/index.html

clean:
	rm -rf *.out *.o s21_matrix_oop.a *.gcda *.gcno *.info 
	-rm -rf report

leaks: test
	leaks --atExit -- ./test.out

clang:
	cp ../materials/linters/.clang-format .
	clang-format -i *.cc *.h
	clang-format -n *.cc *.h
	rm .clang-format
