CC=g++
#FC=g77  -g -r8
FC=gfortran -O3 -fbounds-check -fdefault-real-8
#FC=gfortran -fopenmp -O3 -fbounds-check -fdefault-real-8
#FC=pathf90 -g -static
OBJS= main.o gentab.o DISTILL.o liqDat.o  CORE_ANC.o  CORE_CPP.o  CORE_ECS.o  CORE_MLT.o  CORE_QUI.o  CORE_PR.o	\
      FLASH2.o  IDEALGAS.o  MIX_HMX.o   REALGAS.o  SETUP2.o  TRNS_ECS.o  TRNS_TCX.o	\
      UTILITY.o CORE_BWR.o  CORE_DE.o   CORE_FEQ.o  CORE_PH0.o  CORE_STN.o  FLSH_SUB.o \
      MIX_AGA8.o PROP_SUB.o  SAT_SUB.o  SETUP.o   TRNSP.o TRNS_VIS.o 

#link to appropriate gfortran libraries depending on the system.

multProp.exe: $(OBJS) 
	$(FC) -o multProp.exe $(OBJS) 

main.o: main.f90
	$(FC) -c main.f90

gentab.o: gentab.f90
	$(FC) -c gentab.f90

DISTILL.o: DISTILL.FOR
	$(FC) -c DISTILL.FOR

liqDat.o: liqDat.f90
	$(FC) -c liqDat.f90



CORE_ANC.o: CORE_ANC.FOR
	$(FC) -c CORE_ANC.FOR

CORE_CPP.o: CORE_CPP.FOR
	$(FC) -c CORE_CPP.FOR

CORE_ECS.o: CORE_ECS.FOR
	$(FC) -c CORE_ECS.FOR

CORE_MLT.o: CORE_MLT.FOR
	$(FC) -c CORE_MLT.FOR


CORE_QUI.o: CORE_QUI.FOR
	$(FC) -c CORE_QUI.FOR

CORE_PR.o: CORE_PR.FOR
	$(FC) -c CORE_PR.FOR

FLASH2.o: FLASH2.FOR
	$(FC) -c FLASH2.FOR

IDEALGAS.o: IDEALGAS.FOR
	$(FC) -c IDEALGAS.FOR

MIX_HMX.o: MIX_HMX.FOR
	$(FC) -c MIX_HMX.FOR

REALGAS.o: REALGAS.FOR
	$(FC) -c REALGAS.FOR

SETUP2.o: SETUP2.FOR
	$(FC) -c SETUP2.FOR 

TRNS_ECS.o: TRNS_ECS.FOR
	$(FC) -c TRNS_ECS.FOR

TRNS_TCX.o: TRNS_TCX.FOR
	$(FC) -c TRNS_TCX.FOR

UTILITY.o: UTILITY.FOR
	$(FC) -c UTILITY.FOR

CORE_BWR.o: CORE_BWR.FOR
	$(FC) -c CORE_BWR.FOR

CORE_DE.o: CORE_DE.FOR
	$(FC) -c CORE_DE.FOR

CORE_FEQ.o: CORE_FEQ.FOR
	$(FC) -c CORE_FEQ.FOR

CORE_PH0.o: CORE_PH0.FOR
	$(FC) -c CORE_PH0.FOR

CORE_STN.o: CORE_STN.FOR
	$(FC) -c CORE_STN.FOR

FLSH_SUB.o: FLSH_SUB.FOR
	$(FC) -c FLSH_SUB.FOR

MIX_AGA8.o: MIX_AGA8.FOR
	$(FC) -c MIX_AGA8.FOR

PROP_SUB.o: PROP_SUB.FOR
	$(FC) -c PROP_SUB.FOR

SAT_SUB.o: SAT_SUB.FOR
	$(FC) -c SAT_SUB.FOR

SETUP.o: SETUP.FOR
	$(FC) -c SETUP.FOR

TRNSP.o: TRNSP.FOR
	$(FC) -c TRNSP.FOR

TRNS_VIS.o: TRNS_VIS.FOR
	$(FC) -c TRNS_VIS.FOR

lib:   $(OBJS)
	ar r librefprop.a $(OBJS)
	ranlib librefprop.a

clean:
	rm -f *.o *.exe
