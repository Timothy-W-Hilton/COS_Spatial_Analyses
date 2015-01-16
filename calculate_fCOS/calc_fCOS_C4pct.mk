SHELL=/bin/sh

PROGRAM=calc_fCOS_C4pct

FC=pgf90
LD=pgf90

FFLAGS= -O2 -c -I$(HOME)/local/include 
APINCL=-I$(HOME)/local/include
APILIB= -L$(HOME)/local/lib -lioapi -lnetcdf -ldatetime -lioapi_regrid_tools

CMD=$(PROGRAM).x
SRCDRV=$(PROGRAM).F90
OBJDRV=$(PROGRAM).o

all:  	$(CMD)

$(CMD):  $(OBJDRV)  
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJDRV) \
	$(APILIB) 

%.o: %.F90
	$(FC) $(FFLAGS) $(APINCL) $<

clean:
	rm -f  $(OBJDRV) *.o *.mod

clobber: clean
	rm -f *.x
