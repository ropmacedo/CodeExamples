BASE = $(shell /bin/pwd)
EXED = $(BASE)/$(EXECTBL)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src
INCD = $(BASE)/header

#Compiler of C
CC=gcc

#Flags -I tells which directory to look for the .h files 
CFLAGS =-Wall -I$(INCD) -O3

#Libraries (Math)
LIBS=

#set of headers .h files
DEPS = $(INCD)/$(wildcard *.h)


#set of source .c files
SRC = SDIRK_preconditioner_fin_diff_preconditioner.c initial_guess_fin_diff_preconditioner.c butcher_tableau.c SDIRK_preconditioner.c newton_bicgstab.c output.c coordinates.c problem_specific_functions.c field_and_bound_eqs.c func_and_jacobian.c newton_direct.c solve_equations.c initial_guess.c utilities.c spectral_utilities.c 2+1_Free_Boundaries.c parameters.c derivatives.c debug.c allocate_and_free_derivs.c initial_data.c index_routines.c
  
#set of object .o files
OBJS = $(SRC:.c=.o)

EXECTBL  = Wave

#rule to say that .o files depends on .c and files on DEPS. All those
#files are going to be checked on the execution for the up-to-date version


all:$(EXECTBL)


$(EXECTBL): $(OBJS)
	@echo '' 
	@echo '============ ...building the executable ============'
	@echo ''
	cd $(OBJD); $(CC) $(CFLAGS) -o ../$(EXECTBL) $(OBJS) $(LIBS)
	@echo ''
	@echo '============ ...Done! ============'


vpath %.o $(OBJD)
vpath %.c $(SRCD)
%.o:%.c $(DEPS)
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	cd $(OBJD); $(CC) -o $@ -c $(SRCD)/$(*).c $(LIBS) $(CFLAGS)



flush:
	rm -f $(EXECTBL) *.dat fort.* $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc*

clean:
	rm -f *.dat *.yg fort.* $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc* $(SRCD)/*.o $(EXECTBL)
	rm -f *~ $(OBJD)/*~ $(SRCD)/*~ header/*~
