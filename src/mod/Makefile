src	= \
	nrtype.f90 \
	nrutil.f90 \
	icatools.f90 \
	ran1.f90 \
	utils.f90 

obj	= $(src:.f90=.o)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

default: libfstica.a

libfstica.a: $(obj)
	$(AR) libfstica.a $(obj)
	@cp -p libfstica.a $(LIBDIR)
	@for p in *.$(MOD); do \
	cp $$p $(INCDIR); \
	done
