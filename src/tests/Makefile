# Usage:
#
#    make USE_ASSERTIONS=.false. # Defaults to .true. if not provided

# Example of printing a variable named "var" delimited by square braces
# $(info $$var is [${var}])

ifeq ($(USE_ASSERTIONS),)
  USE_ASSERTIONS:=.false.
endif

ifndef NETCDF
	NETCDF=/usr/local
endif
NCDF_PATH = ${NETCDF}
LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
INCNETCDF = -I$(NCDF_PATH)/include

BUILD=../build
objects_dir=../objects
utilities=../utilities
physics=../physics
tests_dir=.

coarray_icar_root := ../..
root := $(HOME)

lib_dirs := #-L $(root)/lib
inc_dirs := #-I $(root)/inc
install_dir := $(root)/bin


# Options: cray, gnu, intel
ifndef COMPILER
	COMPILER := cray
endif

# Specify the compiler invocation name and flags
ifeq ($(COMPILER),gnu)
  compile := caf
  link := $(compile)
  parallel := -fopenmp -lgomp
  opt := $(parallel) -fcheck=all -Wall -fbacktrace -g -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/
  link_flags := $(parallel) -fcheck=all -Wall -fbacktrace -g
  ifeq ($(MODE),fast)
	  opt := $(parallel) -Ofast -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/
	  link_flags := $(parallel) -Ofast
  endif
  ifeq ($(MODE),profile)
	  opt := $(parallel) -g -Ofast -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/
	  link_flags := $(parallel) -g -Ofast
  endif
endif

ifeq ($(COMPILER),oshfort)
  compile := oshfort
  link := $(compile)
  OSH_FLAGS := -fcoarray=lib -L/glade/u/home/elfanfa/Cheyenne/opencoarrays_openshmem/src/openshmem -lcaf_openshmem
  opt := -fcheck=all -Wall -fbacktrace -g -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/ $(OSH_FLAGS)
  link_flags := -fcheck=all -Wall -fbacktrace -g $(OSH_FLAGS)
  ifeq ($(MODE),fast)
	  opt := -Ofast -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/ $(OSH_FLAGS)
	  link_flags := -Ofast $(OSH_FLAGS)
  endif
endif


ifeq ($(COMPILER),intel)
	compile := ifort
	link := $(compile)
	parallel := -qopenmp -coarray-config-file=./cafconfig.txt -coarray=distributed
	opt := ${parallel} -standard-semantics -check all -traceback -g -fpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -module $(BUILD)/
	link_flags := ${parallel} -standard-semantics -check all -traceback -g

	ifeq ($(MODE),fast)
		opt := ${parallel} -O3 -xHost -fpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -module $(BUILD)/
	    link_flags := ${parallel} -O3 -xHost
	endif
endif

ifeq ($(COMPILER),cray)
      compile := ftn
	  link := $(compile)
      link_flags :=
      opt := -e T -DUSE_ASSERTIONS=$(USE_ASSERTIONS)
	  # note to disable openmp on cray add -h noomp
endif

ifeq ($(COMPILER),nvidia)
  compile := nvfortran
  link := $(compile)
  parallel := -O3 -stdpar=gpu -Minfo
  # parallel := -O3 -stdpar=multicore -Minfo
	convection_parallel :=
	ifeq ($(LOOP),DO_LOOP)
	  convection_parallel :=
	endif
	ifeq ($(LOOP),DC_LOOP)
    # have to build who project to use gpu
	  convection_parallel := -stdpar=gpu -O3 -Minfo
	  # convection_parallel := -stdpar=multicore -O3 -Minfo
	endif
	ifeq ($(LOOP),OMP_LOOP)
	  convection_parallel += -mp
	endif

  # opt := $(parallel) -fcheck=all -Wall -fbacktrace -g -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/
  opt := $(parallel) -I./ -D${LOOP}=true -DDEBUG_VAL=${DEBUGVAL} -traceback -Mchkstk -DNO_COARRAYS -g -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS)  $(BUILD)/
  link_flags := $(parallel) -Wall -g -mp
  ifeq ($(MODE),fast)
	  opt := $(parallel) -O3 -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -J $(BUILD)/
	  link_flags := $(parallel) -O3
  endif
  ifeq ($(MODE),profile)
	  opt := $(parallel) -g -O3 -cpp $(BUILD)/
	  link_flags := $(parallel) -g -O3
  endif
endif


comp_flags := $(opt) $(inc_dirs)
link_flags := $(link_flags) $(lib_dirs)

executable := test-ideal test-initialization test-timer
extra_flags := -g



# Locate required source files
# objects_source := $(shell ls $(objects_dir)/*.f90)
# bare_objects_source := $(subst $(objects_dir)/,,$(objects_source)) # Eliminate paths
#
# # Swap .f90 file extensions for .o extensions
# object_files := $(subst .f90,.o,$(build_source_files))

all: setup test-ideal

tests: test-ideal test-initialization test-timer

# Executable depends on object files
test-ideal: \
			$(BUILD)/convection_type_interface.o		\
			$(BUILD)/convection_exchangeable_implementation.o	\
			$(BUILD)/convection_exchangeable_interface.o \
			$(BUILD)/test-ideal.o					\
			$(BUILD)/domain_implementation.o		\
			$(BUILD)/domain_interface.o				\
			$(BUILD)/exchangeable_implementation.o	\
			$(BUILD)/exchangeable_interface.o		\
			$(BUILD)/timer_implementation.o			\
			$(BUILD)/timer_interface.o				\
			$(BUILD)/assertions_implementation.o	\
			$(BUILD)/assertions_interface.o			\
			$(BUILD)/configuration_implementation.o	\
			$(BUILD)/configuration_interface.o		\
			$(BUILD)/grid_interface.o				\
			$(BUILD)/grid_implementation.o			\
			$(BUILD)/co_utilities.o					\
			$(BUILD)/mp_driver.o					\
			$(BUILD)/mp_thompson.o
	$(link) $^ -o $@ $(link_flags)

test-initialization: \
					 $(BUILD)/test-initialization.o			\
					 $(BUILD)/domain_implementation.o 		\
					 $(BUILD)/domain_interface.o 			\
					 $(BUILD)/exchangeable_implementation.o	\
					 $(BUILD)/exchangeable_interface.o		\
					 $(BUILD)/assertions_implementation.o	\
					 $(BUILD)/assertions_interface.o		\
					 $(BUILD)/configuration_implementation.o\
					 $(BUILD)/configuration_interface.o		\
					 $(BUILD)/timer_implementation.o		\
		 			 $(BUILD)/timer_interface.o				\
					 $(BUILD)/grid_interface.o				\
					 $(BUILD)/grid_implementation.o			\
					 $(BUILD)/co_utilities.o				\
					 $(BUILD)/mp_driver.o					\
					 $(BUILD)/mp_thompson.o
	$(link) $^ -o $@ $(link_flags)

test-timer: \
			$(BUILD)/timer_interface.o	\
			$(BUILD)/test-timer.o		\
			$(BUILD)/timer_implementation.o
	$(link) $^ -o $@ $(link_flags)


install: $(executable)
	mkdir -p $(install_dir)
	cp $(executable) $(install_dir)

# All dependencies are interfaces or libraries (no implementation should depend on another implementation)
$(BUILD)/test-initialization.o: test-initialization.f90 		\
								$(BUILD)/domain_interface.o     \
								$(BUILD)/assertions_interface.o \
								$(BUILD)/mp_driver.o			\
								$(BUILD)/grid_interface.o		\
								$(BUILD)/timer_interface.o
	$(compile) $(comp_flags) -c $< -o $@

# $(BUILD)/test-bcast.o:  test-bcast.f90		\
# 						$(BUILD)/co_utilities.o

$(BUILD)/test-ideal.o: 	test-ideal.f90					\
						$(BUILD)/domain_interface.o 	\
						$(BUILD)/assertions_interface.o \
						$(BUILD)/mp_driver.o			\
						$(BUILD)/grid_interface.o		\
						$(BUILD)/timer_interface.o
	$(compile) $(comp_flags) -c $< -o $@

$(BUILD)/test-timer.o: test-timer.f90 $(BUILD)/timer_interface.o
	$(compile) $(comp_flags) -c $< -o $@


$(BUILD)/mp_thompson.o: 				$(BUILD)/timer_interface.o			\
										$(BUILD)/co_utilities.o

$(BUILD)/mp_driver.o: 					$(BUILD)/mp_thompson.o 				\
										$(BUILD)/convection_exchangeable_interface.o			\
										$(BUILD)/domain_interface.o

$(BUILD)/grid_implementation.o: 		$(BUILD)/grid_interface.o

$(BUILD)/timer_implementation.o: 		$(BUILD)/timer_interface.o

$(BUILD)/assertions_implementation.o: 	$(BUILD)/assertions_interface.o

$(BUILD)/exchangeable_interface.o: 		$(BUILD)/grid_interface.o

$(BUILD)/exchangeable_implementation.o:	convection_flags :=
$(BUILD)/exchangeable_implementation.o: $(BUILD)/exchangeable_interface.o 	\
										$(BUILD)/grid_interface.o

$(BUILD)/convection_type_interface.o: 		$(BUILD)/grid_interface.o

$(BUILD)/convection_exchangeable_interface.o: 		\
										$(BUILD)/convection_type_interface.o 	\
										$(BUILD)/exchangeable_interface.o 	\
										$(BUILD)/grid_interface.o

$(BUILD)/convection_exchangeable_implementation.o: convection_flags :=  $(convection_parallel)
$(BUILD)/convection_exchangeable_implementation.o: \
									  $(BUILD)/convection_exchangeable_interface.o 	\
										$(BUILD)/assertions_interface.o		\
										$(BUILD)/domain_interface.o \
										$(BUILD)/grid_interface.o

$(BUILD)/configuration_implementation.o: $(BUILD)/configuration_interface.o

$(BUILD)/domain_interface.o: 			$(BUILD)/configuration_interface.o 	\
										$(BUILD)/exchangeable_interface.o	\
										$(BUILD)/grid_interface.o

$(BUILD)/domain_implementation.o:	convection_flags :=
$(BUILD)/domain_implementation.o: 		$(BUILD)/domain_interface.o 		\
										$(BUILD)/assertions_interface.o		\
										$(BUILD)/grid_interface.o


$(info $$BUILD        = ${BUILD})
$(info $$tests_dir    = ${tests_dir})
$(info $$objects_dir  = ${objects_dir})
$(info $$utilities    = ${utilities})
$(info $$physics      = ${physics})

# Prefer files in the present working directory (picks the local main over another main)
# %.o: %.f90
# 	$(compile) $(comp_flags) -c $< -o $@

%.o: $(tests_dir)/%.f90
	$(compile) $(comp_flags) -c $< -o $@

# Compile dependencies from utilities directory
%.o: $(utilities)/%.f90
	$(compile) $(comp_flags) -c $< -o $@

# Compile dependencies from physics directory
%.o: $(physics)/%.f90
	$(compile) $(comp_flags) -c $< -o $@

# Compile dependencies from objects
%.o: $(objects_dir)/%.f90
	$(compile) $(convection_flags) $(comp_flags) -c $< -o $@


.PHONY : clean uninstall install tests all

setup:
	mkdir -p ${BUILD}

clean :
	-rm -f  $(BUILD)/*.o $(BUILD)/*.mod $(BUILD)/*.smod $(executable) *.mod
uninstall :
	-rm -f  $(install_dir)/$(executable)
