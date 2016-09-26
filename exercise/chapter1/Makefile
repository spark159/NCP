all: prog1_3 prog1_8

CC=gcc

prog1_3: prog1_3.c
	$(CC) -o $@ $<

prog1_8: prog1_8.c
	$(CC) -o $@ $<

clean:
	rm prog1_3 prog1_8 *~

install: prog1_3
	cp $< $<_copy
