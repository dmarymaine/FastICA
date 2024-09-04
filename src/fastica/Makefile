src	= fastICA.f90
obj	= fastICA.o
exec	= $(BINDIR)/fastICA

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@	

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@	

default: $(exec)

$(exec): $(obj)
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS)
