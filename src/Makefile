SRCS = s21_matrix.c
	
OBJS = ${SRCS:.c=.o}
CC = gcc
CFLAGS = -Wall -Werror -Wextra

all: clean s21_matrix.a test

s21_matrix.a: ${SRCS}
	${CC} ${CFLAGS} -c ${SRCS}
	@ar rcs s21_matrix.a ${OBJS}
	@ranlib s21_matrix.a

test: s21_matrix.a unit_test

unit_test:
	${CC} ${CFLAGS} -c s21_matrix-test.c
	${CC} s21_matrix-test.o s21_matrix.a -lcheck -o test
	./test

gcov_report: s21_matrix.a
	@gcc --coverage s21_matrix-test.c -lcheck ${SRCS} -o unit-test
	@./unit-test
	@lcov -t "test" -o test.info -c -d .
	@genhtml -o report test.info
	
clean:
	@/bin/rm -rf *.o *.a test unit-test *.gcno *gcda report *.info main

checks: cppcheck leaks style

style:
	python3 ../materials/linters/cpplint.py --extensions=c *.c *.h

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem . *.c

leaks: s21_matrix-test.c
	CK_FORK=no leaks --atExit -- ./test
	
main: s21_matrix.a main.c
	${CC} ${CFLAGS} main.c s21_matrix.a -o main

run: main
	./main
