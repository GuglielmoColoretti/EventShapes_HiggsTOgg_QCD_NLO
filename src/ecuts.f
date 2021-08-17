c---- routines for event reconstruction, determination of jet cross sections and 
c---- shape variables



      subroutine distrib(wtdis)
c--- optimization of integration for a particular distribution (used only for moments)
      implicit real*8(a-h,o-z)
      common /jetdata/y45J,y34J,y23J,y45D,y34D,y23D,y45G,y34G,y23G
      common /cuts/ycutJ,ycutD,ycutG,Bcut,Ccut,Fcut,Tcut,Scut,em2hcut
      common /evdata/Cpar,Dpar,Spar,Apar,Planar,Tpar,
     .       Tmajor,Tminor,Opar,em2h,em2l,em2d,
     .       bmax,bmin,bsum,bdiff
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      data init/0/
*
* Events are generated uniformly relative to the function
* programmed in distrib. If this function is choosen unity
* events will be generated according to cross section
*
      if(init.eq.0.and.idist.eq.1)then
      if(iaver.eq.0)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' y23J distribution'
      elseif(iaver.eq.1)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' BW '
      elseif(iaver.eq.2)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' C parameter '
      elseif(iaver.eq.3)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' M_H^2'
      elseif(iaver.eq.4)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' 1-T '
      elseif(iaver.eq.5)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' BT parameter '
      elseif(iaver.eq.6)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' y23D '
      elseif(iaver.eq.7)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' y23J '
      elseif(iaver.eq.8)then
        write(*,*)' ***** WARNING ******'
        write(*,*)' Subroutine DISTRIB optimised for',
     . ' y23D (unweighted) '
      endif
        init=1
      endif
      if(idist.eq.1)then
        wtdis=1d0
      else
        wtdis=1d0
      endif

      return
      end
*
************************************************************************
*
      subroutine ecuts(npar,var,ipass)
c--- application of all event selection cuts
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5)
      common /tcuts/ymin,y0 
      common /cuts/ycutJ,ycutD,ycutG,Bcut,Ccut,Fcut,Tcut,Scut,em2hcut
      common /jetdata/y45J,y34J,y23J,y45D,y34D,y23D,y45G,y34G,y23G
      common /evdata/Cpar,Dpar,Spar,Apar,Planar,Tpar,
     .       Tmajor,Tminor,Opar,em2h,em2l,em2d,
     .       bmax,bmin,bsum,bdiff
      common/intech/iaver,imom,idist,iang,idebug
      ipass=0
      call getjet(y45J,y34J,y23J,npar,1,2)
      call getjet(y45D,y34D,y23D,npar,2,1)
      call getjet(y45G,y34G,y23G,npar,3,1)
      call getCD(Cpar,Dpar,npar)
      call getSAP(Spar,Apar,Planar,npar)
      call getT(Tpar,Tmajor,Tminor,Opar,em2h,em2l,em2d,
     .       bmax,bmin,bsum,bdiff,npar)
      call getvar(var)
	
      if(iaver.eq.0)then
        if(Bmax.gt.Bcut)ipass=1
        if(Bsum.gt.Bcut)ipass=1
        if(Cpar.gt.Ccut)ipass=1
        if(Tmajor.gt.Fcut)ipass=1
	if(1d0-Tpar.gt.Tcut)ipass=1
        if(em2h.gt.em2hcut)ipass=1
        if(ycutJ.gt.y34J.and.ycutJ.lt.y23J)ipass=1
      elseif(iaver.eq.1)then
        if(Bmax.gt.Bcut)ipass=1
      elseif(iaver.eq.2)then
        if(Cpar.gt.Ccut)ipass=1
      elseif(iaver.eq.3)then
        if(em2h.gt.em2hcut)ipass=1
      elseif(iaver.eq.4)then
        if(1-Tpar.gt.Tcut)ipass=1
      elseif(iaver.eq.5)then
        if(Bsum.gt.Bcut)ipass=1
      elseif(iaver.eq.6)then
        if(y23D.gt.ycutD)ipass=1
      elseif(iaver.eq.7)then
        if(y23J.gt.ycutJ)ipass=1
      elseif(iaver.eq.8)then
        if(y23D.gt.ycutD)ipass=1
        if(y34D.gt.ycutD)ipass=1
        if(y45D.gt.ycutD)ipass=1
        if(Bmax.gt.Bcut)ipass=1
        if(Bsum.gt.Bcut)ipass=1
        if(Cpar.gt.Ccut)ipass=1
        if(Tmajor.gt.Fcut)ipass=1
	if(1d0-Tpar.gt.Tcut)ipass=1
        if(em2h.gt.em2hcut)ipass=1
      endif
      return
      end
*
************************************************************************
*
      subroutine getvar(var)
c---- provide event shape variable for inclusive cross section
      implicit real*8 (a-h,o-z)
      common /jetdata/y45J,y34J,y23J,y45D,y34D,y23D,y45G,y34G,y23G
      common /cuts/ycutJ,ycutD,ycutG,Bcut,Ccut,Fcut,Tcut,Scut,em2hcut
      common /evdata/Cpar,Dpar,Spar,Apar,Planar,Tpar,
     .       Tmajor,Tminor,Opar,em2h,em2l,em2d,
     .       bmax,bmin,bsum,bdiff
      common/intech/iaver,imom,idist,iang,idebug
      if(iaver.eq.0)then
         var=y23J
      elseif(iaver.eq.1)then
         var=Bmax
      elseif(iaver.eq.2)then
         var=Cpar
      elseif(iaver.eq.4)then
         var=1d0-Tpar
      elseif(iaver.eq.5)then
         var=Bsum
      elseif(iaver.eq.3)then
         var=em2h
      elseif(iaver.eq.7)then
         var=y23J
      elseif(iaver.eq.8)then
         var=y23J
      elseif(iaver.eq.6)then
         var=y23D
      endif	 
      var = var**imom
      return
      end
*
************************************************************************
*
* jetalg = 1; JADE
* jetalg = 2; DURHAM
* jetalg = 3; GENEVA
* jetcom = 1; E scheme
* jetcom = 2; E0 scheme
* jetcom = 3; P scheme
* jetcom = 4; P0 scheme
*
************************************************************************
*
      subroutine getjet(y45,y34,y23,npar,jetalg,jetcom)
c---- determine jet transition parameters
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5) 
      dimension pjet(4,5)
      evis=0d0
      do i=1,npar
      do j=1,4
        pjet(j,i)=ppar(j,i)
      enddo
      evis=evis+pjet(4,i)
      enddo
      y45=0d0
      y34=0d0
      if(npar.ge.5)then
      y45=1d0
      do i=1,npar-1
        do j=i+1,npar
          if(pjet(4,i).gt.0d0.and.pjet(4,j).gt.0d0)then
            vij=v(pjet(1,i),pjet(1,j))
            ei=pjet(4,i)
            ej=pjet(4,j)
            emin=dmin1(ei,ej)
            if(jetalg.eq.1)ytemp=2d0*ei*ej*vij/evis**2
            if(jetalg.eq.2)ytemp=2d0*emin**2*vij/evis**2
            if(jetalg.eq.3)ytemp=8d0*ei*ej*vij/9d0/(ei+ej)**2/evis**2
            if(ytemp.lt.y45)then
              y45=ytemp
              ii=i
              jj=j
            endif
          endif
        enddo
      enddo
      call jetco(pjet(1,ii),pjet(1,jj),evis,jetcom)
      endif
      if(npar.ge.4)then
      y34=1d0
      do i=1,npar-1
        do j=i+1,npar
          if(pjet(4,i).gt.0d0.and.pjet(4,j).gt.0d0)then
c      vij=1-cos(thetaij)	  
            vij=v(pjet(1,i),pjet(1,j))
            ei=pjet(4,i)
            ej=pjet(4,j)
            emin=dmin1(ei,ej)
            if(jetalg.eq.1)ytemp=2d0*ei*ej*vij/evis**2
            if(jetalg.eq.2)ytemp=2d0*emin**2*vij/evis**2
            if(jetalg.eq.3)ytemp=8d0*ei*ej*vij/9d0/(ei+ej)**2 
            if(ytemp.lt.y34)then
              y34=ytemp
              ii=i
              jj=j
            endif
          endif
        enddo
      enddo
      call jetco(pjet(1,ii),pjet(1,jj),evis,jetcom)
      endif
      y23=1d0
      do i=1,npar-1
        do j=i+1,npar
          if(pjet(4,i).gt.0d0.and.pjet(4,j).gt.0d0)then
            vij=v(pjet(1,i),pjet(1,j))
            ei=pjet(4,i)
            ej=pjet(4,j)
            emin=dmin1(ei,ej)
            if(jetalg.eq.1)ytemp=2d0*ei*ej*vij/evis**2
            if(jetalg.eq.2)ytemp=2d0*emin**2*vij/evis**2
            if(jetalg.eq.3)ytemp=8d0*ei*ej*vij/9d0/(ei+ej)**2 
            if(ytemp.lt.y23)then
              y23=ytemp
              ii=i
              jj=j
            endif
          endif
        enddo
      enddo
      return
      end
*
************************************************************************
*
      function v(a,b)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4)
      pa=sqrt(a(1)**2+a(2)**2+a(3)**2)
      pb=sqrt(b(1)**2+b(2)**2+b(3)**2)
      ea=a(4)
      eb=b(4)
      v  = (dot(a(1),b(1))-ea*eb+pa*pb)/pa/pb
      return
      end
*
************************************************************************
*
      subroutine jetco(a,b,evis,jetcom)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4)
*
* different jet clustering algorithms
*
        if(jetcom.eq.1)then
* e-scheme
          do i=1,4
            a(i)=a(i)+b(i)
          enddo
        elseif(jetcom.eq.2)then
* e0-scheme
          pij=dsqrt((a(1)+b(1))**2
     .             +(a(2)+b(2))**2
     .             +(a(3)+b(3))**2)
          fact=(a(4)+b(4))/pij
          do i=1,3
            a(i)=fact*(a(i)+b(i))
          enddo
          a(4)=a(4)+b(4)
        elseif(jetcom.eq.3)then
* p-scheme
          do  i=1,3
            a(i)=a(i)+b(i)
          enddo
          a(4)=dsqrt(a(1)**2+a(2)**2+a(3)**2)
        elseif(jetcom.eq.4)then
* p0-scheme
          evis=evis-a(4)-b(4)
          do  i=1,3
            a(i)=a(i)+b(i)
          enddo
          a(4)=dsqrt(a(1)**2+a(2)**2+a(3)**2)
          evis=evis+a(4)
        endif
        b(4)=-1d0
        return
        end
*
************************************************************************
*
      subroutine getCD(Cpar,Dpar,npar)
c--- compute C and D parameters
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5)
      t11=theta(1,1,npar)
      t12=theta(1,2,npar)
      t13=theta(1,3,npar)
      t22=theta(2,2,npar)
      t23=theta(2,3,npar)
      t33=theta(3,3,npar)
      Dpar=t11*t22*t33-t11*t23**2-t22*t13**2-t33*t12**2+2d0*t12*t13*t23
      Dpar=27d0*Dpar
      Cpar=-t12**2-t13**2-t23**2+t11*t22+t11*t33+t22*t33
      Cpar=3d0*Cpar
      return
      end
*
************************************************************************
*
      function theta(ii,ij,npar)
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5)
      top=0d0
      bot=0d0
      do i=1,npar
        top=top+ppar(ii,i)*ppar(ij,i)/ppar(4,i)
        bot=bot+ppar(4,i)
      enddo
      theta=top/bot
      return
      end
*
************************************************************************
*
      function dot3(a,b)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4)
      dot3=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
*
************************************************************************
*
      subroutine getSAP(Spar,Apar,Planar,npar)
c--- compute S,A and P parameters (historically, not used in precision studies)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793238d0)
      common /pcut/ppar(4,5)
      t11=phi(1,1,npar)
      t12=phi(1,2,npar)
      t13=phi(1,3,npar)
      t22=phi(2,2,npar)
      t23=phi(2,3,npar)
      t33=phi(3,3,npar)
      x=t11+t22+t33
      y=t11*t22+t11*t33+t22*t33-t12**2-t13**2-t23**2
      z=t11*t22*t33-t11*t23**2-t22*t13**2-t33*t12**2+2d0*t12*t13*t23
      C=3d0*y/x**2
      D=27d0*z/x**3
      sqrt1mc=sqrt(1d0-C)
      ctheta=(D-3d0*C+2d0)/2d0/(1d0-C)/sqrt1mc
      if(ctheta.gt.1d0)ctheta=1d0
      if(ctheta.lt.-1d0)ctheta=-1d0
      theta=acos(ctheta)
      ct=cos(theta/3d0)
      st=sin(theta/3d0)
      rt3=dsqrt(3d0)
      s1=1d0/3d0+2d0/3d0*sqrt1mc*ct
      s2=1d0/3d0-1d0/3d0*sqrt1mc*(1d0*ct+rt3*st)
      s3=1d0/3d0-1d0/3d0*sqrt1mc*(1d0*ct-rt3*st)
      s1=2d0/3d0*sqrt1mc*cos(theta/3d0)+1d0/3d0      
      s2=2d0/3d0*sqrt1mc*cos((theta+2d0*pi)/3d0)+1d0/3d0
      s3=2d0/3d0*sqrt1mc*cos((theta+4d0*pi)/3d0)+1d0/3d0
      e3=dmin1(s1,s2,s3)
      if(e3.lt.0d0)e3=0d0
      e1=dmax1(s1,s2,s3)
      e2=1d0-e1-e3
      Spar=3d0/2d0*(e2+e3)
      Apar=3d0/2d0*e3
      Planar=e2-e3
      return
      end
*
************************************************************************
*
      function phi(ii,ij,npar)
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5)
      phi=0d0
      do i=1,npar
        phi=phi+ppar(ii,i)*ppar(ij,i)
      enddo
      return
      end
*
************************************************************************
*
      subroutine getT(Tpar,Tmajor,Tminor,Opar,em2h,em2l,em2d,
     .       bmax,bmin,bsum,bdiff,npar)
c---- compute thrust T and all observables related to it (MH, BW, BT, Tm, TM, etc.) 
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5)
      common /Tdata/ pp(5,15),ta(3),tmja(3),tmna(3)
      do i=1,15
        do j=1,5
          pp(j,i)=0d0
        enddo
      enddo
      n=0
      do j=1,npar
      n=n+1
      do i=1,4
        pp(i,n)=ppar(i,j)
      enddo
      pp(5,n)=pp(4,n)
      enddo
      if(npar.eq.4)then
        do k=2,npar 
        n=n+1
        do i=1,4
          pp(i,n)=ppar(i,1)+ppar(i,k)
        enddo
        pp(5,n)=sqrt(pp(1,n)**2+pp(2,n)**2+pp(3,n)**2)
        enddo
      elseif(npar.eq.5)then
        do j=1,npar-1
        do k=j+1,npar 
        n=n+1
        do i=1,4
          pp(i,n)=ppar(i,j)+ppar(i,k)
        enddo
        pp(5,n)=sqrt(pp(1,n)**2+pp(2,n)**2+pp(3,n)**2)
        enddo
        enddo
      endif
*
* Find thrust axis and calculate thrust
*
      ithrust=0
      Tpar=0d0
      do i=1,15
         t=2d0*pp(5,i)
         if(t.gt.Tpar) then
            Tpar=t
            ithrust=i
         endif
      enddo
* thrust axis (unit vector)
      do j=1,3
         ta(j) = 2d0*pp(j,ithrust)/Tpar
      enddo
*
* Find major axis (maximizes momentum perp to thrust axis)
*
      imajor=0
      Tmajor=0d0
      do i=1,15
         if(i.ne.ithrust)then
            tmp = pp(1,i)*ta(1)+pp(2,i)*ta(2)+pp(3,i)*ta(3)
            t=2d0*sqrt(abs(pp(5,i)**2-tmp**2))
            if(t.gt.Tmajor)then
               Tmajor=t
               imajor=i
               ppta=tmp
            endif
         endif
      enddo
      if(Tmajor.gt.0) then
         do i=1,3
            tmja(i) = 2d0*(pp(i,imajor) - ta(i)*ppta)/Tmajor
         enddo
*
* Find minor axis
*
         tmna(1) = ta(2)*tmja(3)-ta(3)*tmja(2)
         tmna(2) = ta(3)*tmja(1)-ta(1)*tmja(3)
         tmna(3) = ta(1)*tmja(2)-ta(2)*tmja(1)
         Tminor=0d0
         do i=1,npar
           Tminor = Tminor +
     .     abs(pp(1,i)*tmna(1)+pp(2,i)*tmna(2)+pp(3,i)*tmna(3))
         enddo
* oblateness
         Opar=Tmajor-Tminor
      endif
*
* Hemisphere masses
*
      em2h=pp(4,ithrust)**2-pp(5,ithrust)**2
      em2l=(1d0-pp(4,ithrust))**2-pp(5,ithrust)**2
      if(em2l.gt.em2h) then
         tmp=em2l
         em2l=em2h
         em2h=tmp
      endif
      em2d=em2h-em2l
*
* Hemisphere broadening
*
      bplus=0d0
      bminus=0d0
      bot=0d0
      do i=1,npar
         bot=bot+pp(4,i)
         if(i.ne.ithrust)then
            tmp1 = pp(1,i)*ta(1)+pp(2,i)*ta(2)+pp(3,i)*ta(3)
            tmp2 = sqrt(abs(pp(4,i)**2-tmp1**2))
            if(tmp1.gt.0d0)then
               bplus=bplus+tmp2 
            elseif(tmp1.lt.0d0)then
               bminus=bminus+tmp2 
            endif
         endif
      enddo 
      bplus=bplus/bot/2d0
      bminus=bminus/bot/2d0
      bmax=dmax1(bplus,bminus)
      bmin=dmin1(bplus,bminus)
      bsum =bmax+bmin 
      bdiff=bmax-bmin 
      return
      end
