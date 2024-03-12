!Function written by Eric Lemmon
!modified by Kshitij Neroorkar 
      logical function DISTILL(nc, hf, xkg)
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)   !max number of components in mixture
      dimension x(ncmax), xkg(ncmax) ,xliq(ncmax),xvap(ncmax)
      dimension f(ncmax), xvap2(ncmax) 
      character hrf*3, herr*255
      character*255 hf(ncmax),hfmix
      real*8 nmoles,n(ncmax),nlost,nsum,d60(ncmax),tb(ncmax)
      

      do i=1,nc
        t=288.7055d0
        call LIMITK ('EOS',i,t,D,p,tmin,tmax,Dmax,pmax,ierr,herr)
        if (t.lt.tmin) t=tmin
        x(i)=1
C        write(*,*) t,tmin
        call SATT (t,x,1,p,dl,dv,xliq,xvap,ierr,herr)
        if (ierr.ne.0) write (*,*)tmin,  herr
        d60(i)=dl
        p=101.325d0
        call SATP (p,x,1,tt,dl,dv,xliq,xvap,ierr,herr)
        if (ierr.ne.0) write (*,*) herr
        tb(i)=tt
        x(i)=0
      enddo
        

      call XMOLE (xkg, x, xmolwt)

C Calculate the pure component boiling points curve
C components need to be in order of increasing boiling points
      open (3,file='REFPROP_BPT.dat')
      write(3,*)'# variables = "percent vol evap"  "temperature (C)" ' 
      sumv=0
      do i=1,nc
        sumv=sumv+x(i)/d60(i)
      enddo
      v=0
      t=0
      do i=1,nc-1
        tt=10000
        do k=1,nc
          if (tb(k).gt.t .and. tb(k).lt.tt) then
            j=k
            tt=tb(k)
          endif
        enddo
        t=tb(j)
        write (3,1020) v*100.00,t-273.00
        vfrac=x(j)/d60(j)/sumv
        v=v+vfrac
      enddo
      t=tb(nc)
      write (3,1020) v*100.00,t-273.00
      write (3,1020) 100.00,t-273.00

      close(3)


c     normalize xkg, assign mass fraction
C     this is to make sure that the header shows the correct Mass fraction
      xkgsum=0
      do i=1,nc
        xkgsum=xkgsum+xkg(i)
      enddo
      do i=1,nc
        xkg(i)=xkg(i)/xkgsum
      enddo





c     normalize x, assign mole numbers to n
      xsum=0
      do i=1,nc
        xsum=xsum+x(i)
      enddo
      do i=1,nc
        x(i)=x(i)/xsum
        n(i)=x(i)
      enddo

C  Calculate saturated liquid density at 30 degC 


      t=303.00;
      call SATT(t,x,1,p,dl,dv,xliq,xvap,ierr,herr) 


C  Calculate D86 curve    
      open (2,file='REFPROP_D86.dat')
      write(2,'(A20,F10.5,A7,F10.5,A8)')"#liquid density at",t,"K = ",
     >  (dl*xmolwt),"kg/cu.m"

      write(2,'(A1,1x,100(A20,2x,F6.4))') 
     >            '#',(trim(hf(i)), xkg(i), i=1,nc)  
      write(2,*)'# variables = "percent vol evap"  "temperature (C)" ' 
     

c  start off with 1 mole
      nmoles=1
      vv=0
      p=101.325d0
 
      do j=1,5000
c  find temperature of vapor at liquid composition and p
        call SATP (p,x,1,t,dl,dv,xliq,xvap,ierr,herr)
c  get initial volume of liquid at nmoles=1
        if (j.eq.1) v=1.d0/dl
c  find the amount lost (amount of volume in the vapor compared to that in the
c  liquid)
        nlost=nmoles*dv/dl
        nmoles=nmoles-nlost
c  reduce each mole number by the amount lost times its mole faction in the vapor
        nsum=0
        do i=1,nc
          n(i)=n(i)-nlost*xvap(i)
          nsum=nsum+n(i)
        enddo
        do i=1,nc
          x(i)=n(i)/nsum
        enddo

c  find liquid density from the lost vapor that condensed outside the system
c ?use tt or t??
        call SATP (p,xvap,1,tt,dl,dv,xliq,xvap2,ierr,herr)
c  add the volume lost to the total liquid recovered
        vv=vv+nlost/dl
        write (2,1020) (vv/v *100), t-273.15d0 
c   ,t*9.d0/5.d0-459.67d0
       
        if(vv/v .gt. 1)      stop
      
      enddo

 1010 format (1x,20f8.5)
 1020 format (1x,20f12.5)
      DISTILL=.TRUE.

      return
      end
