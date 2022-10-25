#
#------------------------------------------------------------------------------
#     LINUX 64 bit Definitions
#------------------------------------------------------------------------------
#
linux64:
	@echo
	@echo "---------------------------------------------------------"
	@echo "Configuring for Linux 64 bit..."
	@echo "CMP_F90 = /usr/bin/gfortran -c" > Makefile.config
	@echo "LNK_F90 = /usr/bin/gfortran" >> Makefile.config
#	@echo "FFLAGS = -Bstatic -DCONVERT_BIG -fpp -DMAIN_UCASE -O3 -vec-report0 -I../incl/" >> Makefile.config
	@echo "FFLAGS = -Bstatic -DCONVERT_BIG -DMAIN_UCASE -O3 -I../incl/" >> Makefile.config
	@echo "LIBS = ../libs/liblapack.a ../libs/libskit.a ../libs/libblas.a" >> Makefile.config
	@echo "INCL = -lg2c -lstdc++ -lgfortran" >> Makefile.config
#
	@echo "Building UHMESH for Linux 64 bit..."
	@echo "---------------------------------------------------------"
	@echo
	@cd code/ ; $(MAKE)
	@echo
#
#------------------------------------------------------------------------------
#     Cleaning definitions
#------------------------------------------------------------------------------
#
clean:
	@echo
	@echo "---------------------------------------------------------"
	@echo "Cleaning UHMESH ..."
	@echo
	@cd code/ ; $(MAKE) clean
	rm -f Makefile.config
	@echo "---------------------------------------------------------"
	@echo
