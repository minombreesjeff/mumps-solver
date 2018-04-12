MEX       = mex
MUMPS_DIR = $(HOME)/MUMPS

LMETISDIR = $(HOME)/METIS/metis
LMETIS = -L$(LMETISDIR) -lmetisSUN
LPORDDIR = $(MUMPS_DIR)/PORD/lib
LPORD = -L$(LPORDDIR) -lpordSUN
LORDERINGS = $(LPORD) $(LMETIS)

LIBSEQ  = -L$(MUMPS_DIR)/libseq -lmpiseqSUN

LIBFORT = -L/opt2/SUNWspro7/lib -lsunperf -lfminvai -lfai2 -lfsu -lfmaxvai -lfmaxlai -lfai -lfsumai -lsunmath -lfui

OPTC    = -v -g 
