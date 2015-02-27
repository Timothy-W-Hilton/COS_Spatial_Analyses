SHELL=/bin/sh

PROGRAM=meld_whelan_kettle_fsoils

FC=pgf90
LD=pgf90

FFLAGS= -O2 -c -C -I$(HOME)/local/include
APINCL=-I$(HOME)/local/include
APILIB= -L$(HOME)/local/lib -lioapi -lnetcdf -lioapi_regrid_tools

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
