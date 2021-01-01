#!/usr/bin/env make
all:
	gcc single.c -o Single -lm
	gcc RPA.c -o RPA -lm

clean:
	rm -f Single
	rm -f RPA

