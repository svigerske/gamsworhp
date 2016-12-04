all : gamsworhp

gamsworhp : gamsworhp.o gmomcc.o gevmcc.o

clean:
	rm -f *.o gamsworhp

install: gamsworhp
	./gamsinst.sh

%.c : gams/apifiles/C/api/%.c
	cp $< $@

LDFLAGS = -ldl -Wl,-rpath,$(realpath gams) -Lworhp/lib -lworhp -Wl,-rpath,$(realpath worhp/lib)
CFLAGS = -Igams/apifiles/C/api -g
