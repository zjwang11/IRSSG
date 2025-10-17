# Unified Makefile for IRSSG project

.SUFFIXES: .o .f .f90

# Compiler settings
F90 = ifort
FLAGS = -fpp #-g -traceback -check all
LDFLAGS = $(FLAGS) -qmkl
AR = ar
ARFLAGS = rcs

# Directories
LIBDIR = lib
SRCDIR = src_pw
SRCDIR_WANN = src_wann

# Where to place generated .mod files
MODDIR = build/mod
MODFLAGS = -module $(MODDIR) -I$(MODDIR)
MKDIR_MOD = @mkdir -p $(MODDIR)

# Data path settings (contains kLittleGroups)
# If USE_ABS_DATA_PATH=1, lib_bilbao is compiled with a fixed absolute path,
# otherwise it will read IRVSPDATA env var at runtime (recommended for packaging).
USE_ABS_DATA_PATH ?= 1
IRSSGDATA_ABS := $(abspath $(LIBDIR))
IRSSGDATA_CFLAGS :=
ifeq ($(USE_ABS_DATA_PATH),1)
IRSSGDATA_CFLAGS += -DIRSSGDATA_PATH='"$(IRSSGDATA_ABS)"'
endif

# Library objects (from lib_ssg)
LIB_OBJS = $(LIBDIR)/lib_params.o $(LIBDIR)/lib_comms.o $(LIBDIR)/mathlib.o $(LIBDIR)/lib_bilbao.o $(LIBDIR)/lib_chrct.o \
           $(LIBDIR)/invreal33.o $(LIBDIR)/invmati.o $(LIBDIR)/kgroup.o $(LIBDIR)/irrep_ssg.o \
           $(LIBDIR)/get_ssg.o $(LIBDIR)/comprel.o $(LIBDIR)/linear_rep.o

# Source objects (PW + Wann drivers + unified main)
PW_OBJS = $(SRCDIR)/comms.o $(SRCDIR)/init.o $(SRCDIR)/wave_data_pw.o $(SRCDIR)/main.o
WANN_OBJS = $(SRCDIR_WANN)/file_util.o $(SRCDIR_WANN)/comms.o $(SRCDIR_WANN)/init.o $(SRCDIR_WANN)/wave_data.o $(SRCDIR_WANN)/main.o
SRC_MAIN = $(SRCDIR)/unified_main.o

SRC_OBJS = $(PW_OBJS) $(WANN_OBJS) $(SRC_MAIN)

# Library name
LIBRARY = $(LIBDIR)/libIRSSG.a

# Default target
default: irssg

# Build the main executable
irssg: $(LIBRARY) $(SRC_OBJS)
	$(F90) $(LDFLAGS) -o irssg $(SRC_OBJS) -L$(LIBDIR) -lIRSSG

# Build the library
$(LIBRARY): $(LIB_OBJS)
	$(AR) $(ARFLAGS) $@ $^

# Compile library source files
$(LIBDIR)/lib_params.o: $(LIBDIR)/lib_params.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/lib_comms.o: $(LIBDIR)/lib_comms.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/mathlib.o: $(LIBDIR)/mathlib.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/lib_bilbao.o: $(LIBDIR)/lib_bilbao.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -DIRVSPDATA $(IRSSGDATA_CFLAGS) -o $@ $<

$(LIBDIR)/lib_chrct.o: $(LIBDIR)/lib_chrct.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/invreal33.o: $(LIBDIR)/invreal33.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/invmati.o: $(LIBDIR)/invmati.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/kgroup.o: $(LIBDIR)/kgroup.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/irrep_ssg.o: $(LIBDIR)/irrep_ssg.f90 $(LIBDIR)/lib_chrct.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/get_ssg.o: $(LIBDIR)/get_ssg.f90 $(LIBDIR)/lib_params.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/comprel.o: $(LIBDIR)/comprel.f90 $(LIBDIR)/lib_params.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(LIBDIR)/linear_rep.o: $(LIBDIR)/linear_rep.f90 $(LIBDIR)/lib_params.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

# Compile source files
$(SRCDIR)/comms.o: $(SRCDIR)/comms.f90 $(LIBDIR)/lib_params.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR)/init.o: $(SRCDIR)/init.f90 $(SRCDIR)/comms.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<


$(SRCDIR)/wave_data_pw.o: $(SRCDIR)/wave_data_pw.f90 $(SRCDIR)/comms.o
	$(MKDIR_MOD)
	ifort -g -traceback -assume byterecl $(MODFLAGS) -c -o $@ $<

$(SRCDIR)/main.o: $(SRCDIR)/main.f90 $(SRCDIR)/comms.o $(SRCDIR)/init.o $(SRCDIR)/wave_data_pw.o \
                   $(LIBDIR)/get_ssg.o $(LIBDIR)/linear_rep.o $(LIBDIR)/comprel.o $(LIBDIR)/kgroup.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR)/%.o: $(SRCDIR)/%.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

# Compile Wann sources
$(SRCDIR_WANN)/file_util.o: $(SRCDIR_WANN)/file_util.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR_WANN)/comms.o: $(SRCDIR_WANN)/comms.f90 $(LIBDIR)/lib_params.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR_WANN)/init.o: $(SRCDIR_WANN)/init.f90 $(SRCDIR_WANN)/comms.o $(SRCDIR_WANN)/file_util.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR_WANN)/wave_data.o: $(SRCDIR_WANN)/wave_data.f90 $(SRCDIR_WANN)/comms.o
	$(MKDIR_MOD)
	ifort -g -traceback -assume byterecl $(MODFLAGS) -c -o $@ $<

$(SRCDIR_WANN)/main.o: $(SRCDIR_WANN)/main.f90 $(SRCDIR_WANN)/comms.o $(SRCDIR_WANN)/init.o $(SRCDIR_WANN)/wave_data.o \
                        $(LIBDIR)/get_ssg.o $(LIBDIR)/linear_rep.o $(LIBDIR)/comprel.o $(LIBDIR)/kgroup.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR_WANN)/%.o: $(SRCDIR_WANN)/%.f90
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

$(SRCDIR)/unified_main.o: $(SRCDIR)/unified_main.f90 $(SRCDIR)/main.o $(SRCDIR_WANN)/main.o
	$(MKDIR_MOD)
	$(F90) -c $(FLAGS) $(MODFLAGS) -o $@ $<

# Clean targets
clean:
	rm -f $(LIB_OBJS) $(SRC_OBJS) $(LIBRARY) irssg libIRSSG.a
	rm -f $(LIBDIR)/*.mod $(SRCDIR)/*.mod $(SRCDIR_WANN)/*.mod
	rm -f $(LIBDIR)/*.o $(SRCDIR)/*.o $(SRCDIR_WANN)/*.o
	rm -f ./*.mod ./*.o
	rm -rf $(MODDIR)


.PHONY: default irssg clean
