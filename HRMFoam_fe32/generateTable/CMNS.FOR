c  begin file cmns.f

      subroutine SETINFOdll (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,
     &   acf,dip,Rgas)
c
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_SETINFOdll"::SETINFOdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::SETINFOdll
c     dll_export SETINFOdll
c
      include 'commons.for'
c
      if (abs(icomp).le.nc) then
        wm(icomp)=wmm
        ttpn(icomp)=ttrp
        tnbp(icomp)=tnbpt
        tcrit(icomp)=tc
        pcrit(icomp)=pc
        Dcrit(icomp)=Dc
        Zcrit(icomp)=Zc
        accen(icomp)=acf
        dipole(icomp)=dip
        R=Rgas
      else
        wmm=0.0d0
        ttrp=0.0d0
        tnbpt=0.0d0
        tc=0.0d0
        pc=0.0d0
        Dc=0.0d0
        Zc=0.0d0
        acf=0.0d0
        dip=0.0d0
        Rgas=R
      end if
c

      end
c ======================================================================
      subroutine GETINFOdll (icomp,wmm,ttrp,tnbpt,
     &                       tc,pc,Dc,Zc,acf,dip,Rgas)
c
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_GETINFOdll"::GETINFOdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::GETINFOdll
c     dll_export GETINFOdll
c
      include 'commons.for'
c
      call INFO (icomp,wmm,ttrp,tnbpt,tc,pc,Dc,Zc,acf,dip,Rgas)
      end
c ======================================================================
      subroutine GETCMNdll (tc,rhoc,pc,tred,rhored,
     &                      n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
c
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_GETCMNdll"::GETCMNdll
C  cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_GETCMNdll@48"::GETCMNdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::GETCMNdll
c     dll_export GETCMNdll
c
      include 'commons.for'
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
c
      do i=1,mxcoef
        coefhmx(i)=0.d0
        coefcp0(i)=0.d0
        coefxk(i)=0.d0
      enddo
      tred=tz(1)
      rhored=rhoz(1)
      tc=tcfeq(1)
      rhoc=rhocfeq(1)
      pc=pcfeq(1)
      n1=ntermf(1)
      n2=ncrt(1)
      do i=1,n1+n2
        coefhmx(i)=afeq(1,i)
      enddo
      n3=ntermc(1)
      n4=nterme(1)
      do j=1,n3+n4
        coefcp0(j)=cpc0(1,j)
        coefxk(j)=xkc0(1,j)
      enddo
      end
c ======================================================================
      subroutine SETCMNdll (tc,rhoc,pc,tred,rhored,
     &                      n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
c
cDEC$ ATTRIBUTES DLLEXPORT, Alias: "_SETCMNdll"::SETCMNdll
C  cDEC$ ATTRIBUTES STDCALL, REFERENCE::SETCMNdll
c     dll_export SETCMNdll
c
      include 'commons.for'
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
c
      if (pc.gt.0.d0) pcfeq(1)=pc
      if (tc.gt.0.d0) tcfeq(1)=tc
      if (rhoc.gt.0.d0) rhocfeq(1)=rhoc
      if (tred.gt.0.d0) tz(1)=tred
      if (tred.gt.0.d0) t0feq(1)=tc
      if (rhored.gt.0.d0) rhoz(1)=rhored
      if (rhored.gt.0.d0) rho0feq(1)=rhoc
      if (n1.ne.0) then
        ntermf(1)=n1
        ncrt(1)=n2
        do i=1,n1+n2
          afeq(1,i)=coefhmx(i)
        enddo
        call SETEXP (1)
      endif
      if (n3.ne.0) then
        tsvcp0(1)=0.d0
        ntermc(1)=n3
        nterme(1)=n4
        do j=1,n3+n4
          cpc0(1,j)=coefcp0(j)
          xkc0(1,j)=coefxk(j)
        enddo
      endif
      end
c ======================================================================
      subroutine GETCMN (tc,rhoc,pc,tred,rhored,acf,
     &                   n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
c
cDEC$ ATTRIBUTES DLLEXPORT::GETCMN
c     dll_export GETCMN
c
      include 'commons.for'
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
c
      do i=1,mxcoef
        coefhmx(i)=0.d0
        coefcp0(i)=0.d0
        coefxk(i)=0.d0
      enddo
      tred=tz(1)
      rhored=rhoz(1)
      tc=tcfeq(1)
      rhoc=rhocfeq(1)
      pc=pcfeq(1)
      n1=ntermf(1)
      n2=ncrt(1)
      acf = accen(1)
      do i=1,n1+n2
        coefhmx(i)=afeq(1,i)
      enddo
      n3=ntermc(1)
      n4=nterme(1)
      do j=1,n3+n4
        coefcp0(j)=cpc0(1,j)
        coefxk(j)=xkc0(1,j)
      enddo
      end
c ======================================================================
      subroutine SETCMN (tc,rhoc,pc,tred,rhored,acf,
     &                   n1,n2,coefhmx,n3,n4,coefcp0,coefxk)
c
cDEC$ ATTRIBUTES DLLEXPORT::SETCMN
c     dll_export SETCMN
c
      include 'commons.for'
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
c
c  (re)initialize contents of /FEQSAV/ and /CPPSAV/ when a fluid's parameter are reset
      do 120 i=n0,nx
      delsvf(i)=0.0d0
      tausvf(i)=0.0d0
      cp0sav(i)=0.0d0
      cpisav(i)=0.0d0
      cptsav(i)=0.0d0
      tsvcp0(i)=0.0d0
      do 100 j=1,mxtrm
      phisvf(i,j)=0.0d0
      taup(i,j)=0.0d0
      delp(i,j)=0.0d0
      delli(i,j)=0.0d0
  100 continue
  120 continue
c
      call RESETA
      if (pc.gt.0.d0) pcfeq(1)=pc
      if (tc.gt.0.d0) tcfeq(1)=tc
      if (rhoc.gt.0.d0) rhocfeq(1)=rhoc
      if (tred.gt.0.d0) tz(1)=tred
      if (tred.gt.0.d0) t0feq(1)=tc
      if (rhored.gt.0.d0) rhoz(1)=rhored
      if (rhored.gt.0.d0) rho0feq(1)=rhoc
      if (abs(acf).gt.1.d-20) accen(1)=acf
      if (n1.ne.0) then
        ntermf(1)=n1
        ncrt(1)=n2
        do i=1,n1+n2
          afeq(1,i)=coefhmx(i)
        enddo
        call SETEXP (1)
      endif
      if (n3.ne.0) then
        tsvcp0(1)=0.d0
        ntermc(1)=n3
        nterme(1)=n4
        do j=1,n3+n4
          cpc0(1,j)=coefcp0(j)
          xkc0(1,j)=coefxk(j)
        enddo
      endif
      end
c ======================================================================
      subroutine GETALLCMN (tc,rhoc,pc,tred,rhored,acf,
     &                   n1,n2,coefhmx,tauexp,denexp,denlnexp
     &                      ,n3,n4,coefcp0,coefxk)
c
cDEC$ ATTRIBUTES DLLEXPORT::GETALLCMN
c     dll_export GETALLCMN
c
      include 'commons.for'
      dimension coefhmx(mxcoef),coefcp0(mxcoef),coefxk(mxcoef)
      dimension tauexp(mxcoef),denexp(mxcoef),denlnexp(mxcoef)
c
      do i=1,mxcoef
        coefhmx(i)=0.d0
        tauexp(i)=0.d0
        denexp(i)=0.d0
        denlnexp(i)=0.d0
        coefcp0(i)=0.d0
        coefxk(i)=0.d0
      enddo
      tred=tz(1)
      rhored=rhoz(1)
      tc=tcfeq(1)
      rhoc=rhocfeq(1)
      pc=pcfeq(1)
      acf = accen(1)
      n1=ntermf(1)
      n2=ncrt(1)
      do i=1,n1+n2
        coefhmx(i)=afeq(1,i)
        tauexp(i)=tifeq(1,i)
        denexp(i)=difeq(1,i)
        denlnexp(i)=dlif(1,i)
      enddo
      n3=ntermc(1)
      n4=nterme(1)
      do j=1,n3+n4
        coefcp0(j)=cpc0(1,j)
        coefxk(j)=xkc0(1,j)
      enddo
      end

c ======================================================================
      subroutine GETTYPE(htype)
c
cDEC$ ATTRIBUTES DLLEXPORT::GETTYPE
c     dll_export GETTYPE
c
      include 'commons.for'
      character*3 htype
c
      htype=heos
      end

c ======================================================================
      subroutine LIMITSET (tmin,tmax,Dmax,pmax)
c
cDEC$ ATTRIBUTES DLLEXPORT::LIMITSET
c     dll_export LIMITSET
c
      include 'commons.for'
c
c     reset limit only for the case of a pure component
      if (nc.eq.1) then
c  special case--pure component
        tmeos(1)=tmin
        txeos(1)=tmax
        Deos(1)=Dmax
        peos(1)=pmax
      end if
      end
C END LIMITSet
