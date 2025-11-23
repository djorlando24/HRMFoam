! Calls the subroutines for calculating the Distillation curve,
! the fluid.dat, and the liquid.dat - Kshitij Neroorkar

Program main
implicit none
double precision ::p_st,p_end,h_st,h_end 
double precision :: z(20), zkg(20), molwt, step_p, step_h
integer :: nc,ierr,i,flag
character  :: hrf*3,herr*255,hfmix*255,junk*10
character*255 :: hfiles(20)
logical random
logical generateTable, DISTILL, liqDat



!Reading in the inputs for generating table 
Open(unit=10, FILE = "table.info", Form = "Formatted", Action = "Read")

read(10,*)
read(10,*)flag
read(10,*)
read(10,*)nc  !number of components, transport properties of which component
read(10,*) (hfiles(i) , i=1,nc) !Fluid selection
read(10,*) (zkg(i) , i=1,nc) !moleFractions of each component
read(10,*)junk,p_st,junk,p_end,junk,step_p !Pressure range and step
read(10,*)junk,h_st,junk,h_end,junk,step_h !Enthalpy range and step


!reference state for thermodynamic calculations
!'DEF' - default reference state as specified in fluid file
hrf="DEF"

!file including mixture coefficients
hfmix="fluids/HMX.BNC"



!Set path and initialize to obtain fluid information
call SETPATH('.')


call SETUP(nc,hfiles,hfmix,hrf,ierr,herr) 

if (ierr .ne. 0) then
   write(*,*)'Error msg: ',herr
   if((ierr .ne.(-117)) .and. (ierr .ne.(-28))) then
      write(*,*)'Exiting'
      stop
   endif
endif


!convert from mass fraction to model fraction
call XMOLE (zkg, z, molwt)

!this is required according to Eric Lemmon, it will be removed in the future
!Call SATDATA(z)
!Call SATDATA(z)


call SETREF(hrf, 2, z, 0.0, 0.0,0.0, 0.0, ierr, herr)


!Check for error in initialization
if (ierr .ne. 0) then
   write(*,*)'Error msg:',herr
   if((ierr .ne.(-117)) .and. (ierr .ne.(-28))) then
   write(*,*)'Exiting'
     stop
   endif 
endif


if (flag .eq.1)   then

    random = DISTILL(nc, hfiles, zkg)
endif

if (flag .eq. 2) then

    random = generateTable(nc, hfiles, zkg,z,molwt,p_st,p_end, step_p, h_st,h_end,step_h)

!    random = liqDat(nc, hfiles, zkg, z, molwt)

!    random = DISTILL(nc, hfiles, zkg) 
endif

if (flag .eq.3) random = liqDat(nc, hfiles, zkg, z, molwt)

end program main




