
!subroutine to modify liquid.dat for spray simulations
!Written by Kshitij Neroorkar
logical function liqDat(nc, hfiles, zkg, z, molwt)
Implicit none


double precision :: z(20), zkg(20) ,x(20),y(20),xkg(20),ykg(20),ybub(20),xdew(20)
double precision :: tcrit, pcrit, Dcrit, molwt, t,p, rhol, rhov
double precision :: tmin, tmax, Dmax, pmax, eta, tcx, sigma, e, hl, hv, s, Cv, Cpl
double precision :: w, hjt, Cpv,temp,pv, tprev, rholConst
integer :: nc,ierr,nroot,i, nopoints,kph,j
character  :: hrf*3,herr*255,hfmix*255,junk*10, header1*100, header2*100, header3*100, header4*100
character*255 :: hfiles(20)
logical SatErr

!write fluid.dat header
!Open(unit=11, FILE = "liquid.dat" , Form = "Formatted", Action = "Read")

!read(11,*) 

Open(unit=21, FILE = "liquid_new.dat" , Form = "Formatted", Action = "Write")



!get critical properties
call CRITP(z,tcrit, pcrit,Dcrit,ierr,herr)
if (ierr .ne. 0) then
   write(*,*) "Critical pressure value", pcrit
   write(*,*) ierr, herr
!   stop
endif



write(21,1011)"!Temperature  viscosity surf tension laten heat    vapor pres  cond        density    spec heat"
write(21,1011)"!    K        (N.s/m^2)      (N/m)     (j/kg)       (pascal)   (w/m.k)     (g/m^3)    (j/kg.k)"
1011 Format (A)

write(21,*)"  ic8h18"
write(21,1012)tcrit
1012 Format (2X,F10.6)

nopoints = tcrit/10+2

call LIMITS ('EOS',z,tmin,tmax,Dmax,pmax)


t=0
temp =t
kph =1
SatErr = .false.

!converge requires number of points so since starting from 0 we go to nopoints-1
do i=0,nopoints-1

   tprev =temp

!keep temp between tmin and tcrit
   if (t < tmin) then
       temp =tmin
      else if (t > tcrit) then
           temp =tcrit
           else
               temp=t
   endif

!density needs to be at a constant temperature for spray
!taking at 15 degC since thats where the experiments calculate density.
!DD Feb2024 changed the constant to make it explicitly Real(8)

   CALL SATT(288D0,z,kph,pv,rholConst,rhov,x,y,ierr,herr)
 
!vapor pressure 
   CALL SATT(temp,z,kph,pv,rhol,rhov,x,y,ierr,herr)
   if (ierr .ne. 0) then
       write(*,*) "inside SATT"
       write(*,*) ierr, herr

         
!if error =131 means that satt is not converged
!in this case just calculate properties at the last converged point
       if(ierr .eq. -131) then
          temp =tprev
          CALL SATT(temp,z,kph,pv,rhol,rhov,x,y,ierr,herr)
          SatErr = .true.
       endif
   endif

!viscosity and thermal conductivity
   CALL TRNPRP(temp,rhol,z,eta,tcx,ierr,herr)
   if (ierr .ne. 0) then
       write(*,*) "inside TRNPRP"
       write(*,*) ierr, herr
   endif

!surface tension
   CALL SURTEN(temp,rhol,rhov, x,y, sigma,ierr,herr)
    if (ierr .ne. 0) then
       write(*,*) "inside SURTEN"
       write(*,*) ierr, herr
   endif

!specific heat and enthalpy  of  liquid 
   CALL THERM(temp, rhol, x, p, e, hl, s, Cv, Cpl, w, hjt) 

!enthalpy of vapor for enthalpy of vaporization
   CALL THERM(temp, rhov, y, p, e, hv, s, Cv, Cpv, w, hjt) 

   write(21,1013) t, (eta*1E-6), sigma, ((hv-hl)/molwt*1000.0), pv*1000.0, tcx, rholConst*molwt, Cpl/molwt*1000.0
1013 Format (8(2X,E10.4))

   if (SatErr .eqv. .true.) exit
!   write(21,*)t,  p*1000.0
   t=t+10

enddo


liqDat = .TRUE.

close(21)

return
end 

