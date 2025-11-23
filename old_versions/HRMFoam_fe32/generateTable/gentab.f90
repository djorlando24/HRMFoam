
!subroutine to generate property lookup table from REFPROP
!Written by Shiva Gopalakrishnan <sgopalak@ecs.umass.edu>
!modified by Kshitij Neroorkar for hydrocarbon  mixtures
logical function generateTable(nc, hfiles, zkg, z, molwt, p_st, p_end, step_p, h_st, h_end, step_h)
Implicit none


double precision ::p_st,p_end,p,h_st,h_end,h, hres, pres,t1,p1,d1,t2,p2,d2
double precision :: z(20), zkg(20) ,x(20),y(20),xkg(20),ykg(20),ybub(20),xdew(20)
double precision :: pressure, hin, molwt, tbub,tdew, rholbub, rhovdew, rhol, rhov
double precision :: t,qmol, qkg,eta,rho,tcx,e,s,cp,cv,temp,w,step_p,step_h,tcxl,tcxv,etal,etav,alpha
double precision:: ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas, psat, rhol1, rhov1
double precision:: tcrit, pcrit, Dcrit, molwtliq, molwtvap
integer :: nc,ierr,nroot,k1,k2,kph, PHph,i, hcount
integer :: resolution, ntrans
character  :: hrf*3,herr*255,hfmix*255,junk*10
character*255 :: hfiles(20)




!Calculating number of data points
 pres=(p_end - p_st)/step_p + 1
 !resolution=pres
 hres=(h_end -h_st)/step_h + 1
 resolution=int(pres)*int(hres)
 p=p_st
 h=h_st


!write fluid.dat header
Open(unit=11, FILE = "fluid.dat" , Form = "Formatted", Action = "Write")

write(11,'(100(A20,2x,F6.4))')(trim(hfiles(i)), zkg(i), i=1,nc) ,' Fluid Data'
write(11,'(A,f13.1,A, f13.1,A,f13.1)')'P_start',p_st,' P_end',p_end,' Step',step_p
write(11,'(A,f13.1,A, f13.1,A,f13.1)')'H_start',h_st,' H_end',h_end,' Step',step_h
write(11,*)'Points',resolution,' P_res',int(pres),' H_res',int(hres)
!get critical properties
call CRITP(z,tcrit, pcrit,Dcrit,ierr,herr)


if (ierr .ne. 0) then
   write(*,*) "Critical pressure value", pcrit
   write(*,*) ierr, herr
!   stop
endif


write(11,*)'CriticalPressure', pcrit*1000

write(11,'(A)')'Pressure(Pa)     Enthalpy(J/Kg)     Quality        &
Density(kg/m3)  DensityVapor    DensityLiquid   &
Temperature(K)  Viscosity(Pa.s)      Conductivity     Cv(KJ/KgK)     Cp(KJ/KgK)    SoS(m/s)  PSat(Pa)'


!Loop over pressure
Do while (p <= p_end)
h=h_st

!Loop over entahlpy
!Do while (h <= h_end)
!$OMP PARALLEL DO
Do hcount = 1,int(hres)


!Flag for phas in PHFL1 calculations , used only when PHFLSH crashes
   PHph = 0

!initialize values
   cv=0.00
   cp=0.00
   w=0.00


!Input needs to in KPa and KJ/KgMol

h = h_st + step_h*(hcount-1)
   pressure=p/1000
   hin=h*molwt/1000

!Calculating thermodynamic properties*********************************************************************
      

   call PHFLSH(pressure,hin,z,t,rho,rhol,rhov,x,y,qmol,e,s,cv,cp,w,ierr,herr)
   
   

!error 242 means that the PHFLSH exited with an error because the point was too close to critical point
!force liquid state calculations in this case
   if(ierr .ne.0) then
      if(ierr.eq. 242) then
!Phph=1 means liquid
         PHph=1
!PHFL1 is a single phase equivalent of PHFLSH
         call PHFL1(pressure, hin, z, PHph, t, rho, ierr, herr)
         qmol = 0.0
         if (ierr .ne. 0) then
            write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
            write(*,*) ierr, herr
         endif

!error 224 means that the PHFLSH exited with an error because dew point calculation did not converge within 2-phase iteration
!force vapor state calculations in this case

      else if(ierr.eq. 224) then
!PHph=2 means vapor
         PHph=2
!PHFL1 is a single phase equivalent of PHFLSH
         call PHFL1(pressure, hin, z, PHph, t, rho, ierr, herr)
         qmol = 1.0
 
         if(ierr .ne. 0) then
            write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
            write(*,*) ierr, herr
         endif
      else
         write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
         write(*,*)"PHFlash error ", ierr, herr
         stop
      endif
   endif


!here quality (qmol) is in molar basis so converting it to qkg**********************************************

   qmol=max(min(1.0,qmol),0.0)

   call QMass(qmol,x,y, qkg, xkg, ykg, molwtliq, molwtvap, ierr, herr)

   qkg=max(min(1.0,qkg),0.0)


   if (ierr .ne. 0) then
      write(*,*) "Pressure = ", pressure *1000, " Pa" , " Enthalpy = ", hin /molwt*1000, " J/kg"
      write(*,*) ierr, herr
      stop
   endif

   kph=1
!calculating pSat********************************************************************************************
  
   call SATH(hin, z, kph, nroot, k1, t1, p1, d1, k2, t2, p2, d2, ierr, herr)


   if (ierr .ne. 0) then
      write(*,*) "In SATH"
      write(*,*) "Pressure = ", pressure *1000, " Pa" , " Enthalpy = ", hin /molwt*1000, " J/kg"
      write(*,*) ierr, herr
      stop
   endif

   !if the first root is a liquid root, psat is equal to p1
   if (k1 .eq. 1) then
      psat= p1
   else
      write(*,*)"SATH error first root is not liquid ", herr     
      write(*,*)"Exiting"
      stop    
   endif

!Calculating rhoL and rhoV at P_sat***************************************************************************

   kph=1


   if (qkg .gt. 0.0) then
!if in two phase region, calculate rhol and rhov based on pressure
      call SATP(pressure,z,kph,temp,rhol,rhov,x,y,ierr,herr)
    
      if (ierr .ne. 0) then
         write(*,*) "In SATP"
         write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
         write(*,*) ierr, herr
         stop
      endif

   else
!if in the subcooled region, the liquid and vapor densities will be calculated
!at the saturation pressure
      call SATP(psat,z,kph,temp,rhol,rhov,x,y,ierr,herr)
      if (ierr .ne. 0) then
         write(*,*) "In SATP"
         write(*,*) "Pressure = ", pressure *1000 , " Pa "  , " Enthalpy = ", hin /molwt*1000, " J/kg"
         write(*,*) ierr, herr
         stop
      endif
   endif



!since some of the transport properties are not available in REFPROP, we just use the transport properties of the first component
!call PUREFLD(2) 

!if in two-phase region, calculate liquid and vapor tansport properties separately and manually calculate two-phase properties
   if (qkg .gt. 0.0) then
!Calculating transport properties for saturated liquid and vapor and then
!interpolating using the quality or void fraction
      call TRNPRP(t,rhol,x,etal,tcxl,ierr,herr)

      if (ierr .ne. 0) then
         write(*,*) "In TRNPRP for rhol"
         write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
         write(*,*)ierr, herr
         stop
      endif

      call TRNPRP(t,rhov,y,etav,tcxv,ierr,herr)

      if (ierr .ne. 0) then
         write(*,*) "In TRNPRP for rhov"
         write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
         write(*,*) ierr, herr
         stop
      endif

      tcx=(1-qkg)*tcxl + qkg*tcxv
      alpha=max(min((rhol-rho)/(rhol-rhov),1.0),0.0)
      eta=(1-alpha)*etal + alpha*etav
   else

!transport properties will be calculated at the temperature and density 
!obtained from PHFLSH calculation     
      call TRNPRP(t,rho,x,eta,tcx,ierr,herr)

      if (ierr .ne. 0) then
          write(*,*)"TRNPRP error"
          write(*,*) "Pressure = ", pressure *1000, " Pa " , " Enthalpy = ", hin /molwt*1000, " J/kg"
          write(*,*) ierr, herr
          stop
      endif

   endif


!call PUREFLD(0) 


!Order of write - pressure,enthalpy,xbar,density,density_vapor,density_liquid,temperature,viscosity,therrmal conductivity,cv,cp and speed of sound
write(11,'(13(x,e15.6))')p,h,qkg,rho*molwt,rhov*molwt,rhol*molwt,t,eta, tcx, cv/molwt, cp/molwt, w, psat*1000 


h=h+step_h
end do

if (h_end /= h_st) write(6,'(a,f8.3,2a)')'Percent complete %:',100.0*(p-p_st)/(p_end-p_st),char(27),char(77)


p=p+step_p
end do

close(11)
close(10)

write(6,*)


generateTable = .TRUE.
return
end 

