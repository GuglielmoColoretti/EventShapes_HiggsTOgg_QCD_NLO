NAME1    = eerad3
NAME2    = eerad3_combine
NAME3    = eerad3_dist


SOURCEDIR = ./src
OBJDIR = ./obj

VPATH = $(SOURCEDIR)


FFILES1   = eerad3.f histo.f ecuts.f eerad3lib.f phaseee.f sigHG.f aversub0H.f
FFILES2   = eerad3_combine.f
FFILES3   = eerad3_dist.f


#for gfortran compiler 
FC        = gfortran  -ffixed-form  -ffixed-line-length-none
FFLAGS    = -fno-automatic -O
#for ifort compiler
#FC        = ifort
#FFLAGS    = -save -O4


OBJFILES1 = $(addprefix $(OBJDIR)/,$(patsubst %.f,%.o,$(FFILES1)))
OBJFILES2 = $(addprefix $(OBJDIR)/,$(patsubst %.f,%.o,$(FFILES2)))
OBJFILES3 = $(addprefix $(OBJDIR)/,$(patsubst %.f,%.o,$(FFILES3)))

$(OBJDIR)/%.o:	%.f
	$(FC) $(FFLAGS) -c $< -o $@

$(NAME1): $(OBJFILES1)
	$(FC) $(FFLAGS) -o $@ $(OBJFILES1)
$(NAME2): $(OBJFILES2)
	$(FC) $(FFLAGS) -o $@ $(OBJFILES2)
$(NAME3): $(OBJFILES3)
	$(FC) $(FFLAGS) -o $@ $(OBJFILES3)

clean:
	rm -f obj/* eerad3 eerad3_dist eerad3_combine


