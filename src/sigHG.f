c---- these routines supply the integrands for the vegas integrations in the main program

************************************************************************
*
      function sig3a(x,wgt)
      implicit real*8(a-h,o-z)
      logical plot 
      dimension x(10)
      common /plots/plot
      common /s3/s12,s13,s23
      common /pcut/ppar(4,5)
      common/eventmom3p/pevt3(4,3,1)
      common/eventinv3p/sevt3(3,3,1),ievt3(1)
      parameter(psconv3=3968.80341507837698d0) ! 128*pi^3
c      pscor3 = 0.5d0*s123
      sig3a=0d0
      call phase3ee(x,wtps,ifail)
      if(ifail.eq.1)return
      call fillcommon3pee
      do i=1,3
         do j=1,4
            ppar(j,i) = pevt3(j,i,1)
         enddo 
      enddo
      call ecuts(3,var,ipass)
      if(ipass.eq.0)return
      s12 = sevt3(1,2,1)
      s13 = sevt3(1,3,1)
      s23 = sevt3(2,3,1)
      sig3a=psconv3*sig3(s12,s13,s23)
      sig3a=sig3a*wtps*var

      if(plot)then
         call bino(1,sig3a*wgt,3)
      else
         call distrib(wtdis)
         sig3a=sig3a/wtdis
      endif
      return
      end
*
************************************************************************
*
      function sig4a(x,wgt)
      implicit real*8(a-h,o-z)
      logical plot 
      dimension x(10)
      parameter(i1=1,i2=2,i3=3,i4=4,i5=5)
      common /plots/plot
      common /yij4/y(4,4)
      common /pcut/ppar(4,5)
      common /pmom/p(4,5) 
      common/eventmom4p/pevt4(4,4,12)
      common/eventinv4p/sevt4(4,4,12),ievt4(12)
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      parameter(psconv4=626728.314440256416d0) ! 2048*pi^5
c	pscor4 = 1d0/12d0*s1234^2 or 0.25d0*s1234
      sig4asum = 0d0
      sig4a = 0d0
      call phase4ee(x,wtps,ifail)
      wtps = wtps*psconv4
      if (iang.eq.1) wtps=wtps/2d0
      if(ifail.eq.1)return
      call fillcommon4pee
      if (iang.eq.0) ievtmax = 6
      if (iang.eq.1) ievtmax = 12
      do ievt=1,ievtmax
         sig4w = 0d0
         sig3sw = 0d0
         if (ievt4(ievt).eq.1) then
            do i=1,4
               do j=1,4
                  ppar(j,i) = pevt4(j,i,ievt)
                  p(j,i) = pevt4(j,i,ievt)
               enddo
            enddo
            do i=1,4
               do j=1,4
                  y(i,j) = sevt4(i,j,ievt)
               enddo
            enddo
            call ecuts(4,var,ipass)
            if (ipass.ne.0) then
               wtplot = wtps*wgt
               sig4w=sig4(i1,i2,i3,i4,wtplot,var)
               sig4w=sig4w*wtps*var
               if (njets.eq.3) then
                  sig3sw=sig3s(i1,i2,i3,i4,wtplot) 
                  sig3sw=sig3sw*wtps
               endif
               sig4w=sig4w-sig3sw
            endif
         endif
         sig4asum = sig4asum + sig4w
      enddo
      sig4a = sig4asum
      return
      end
*
************************************************************************
*
      function sig5a(x,wgt)
      implicit real*8(a-h,o-z)
      dimension x(10)
      logical plot 
      common /plots/plot
      common /phase/ips 
      common /yij5/y(5,5)
      common /pcut/ppar(4,5)
      common /pmom/p(4,5) 
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common/invarang5/sijang5(5,5,4),iacc(4)
      common/eventmom5ap/pevt5a(4,5,120)
      common/eventinv5ap/sevt5a(5,5,120),ievt5a(120)
      common/eventmom5bp/pevt5b(4,5,60)
      common/eventinv5bp/sevt5b(5,5,60),ievt5b(60)
      parameter(i1=1,i2=2,i3=3,i4=4,i5=5)
      parameter(psconv5=2.4742242121947481d7) ! 8192*pi^7
      sig5asum = 0d0
      sig5a=0d0
      if (ips.eq.1) then
         call phase5aee(x,wtps,ifail)
         if(ifail.eq.1)return
         wtps = wtps*psconv5
         call fillcommon5apee
         if (iang.eq.0) ievtmax = 30
         if (iang.eq.1) ievtmax = 120
      elseif (ips.eq.2) then
         call phase5bee(x,wtps,ifail)
         if(ifail.eq.1)return
         wtps = wtps*psconv5
         call fillcommon5bpee
         if (iang.eq.0) ievtmax = 15
         if (iang.eq.1) ievtmax = 60
      endif
      if (iang.eq.1) wtps=wtps/4d0
      do ievt=1,ievtmax
         sig5w=0d0
         sig4sw=0d0
         sig3dsw=0d0
         if (ips.eq.1) ievt5 = ievt5a(ievt)
         if (ips.eq.2) ievt5 = ievt5b(ievt)
         if (ievt5.eq.1) then
            do i=1,5
               do j=1,4
                  if (ips.eq.1) ppar(j,i) = pevt5a(j,i,ievt)
                  if (ips.eq.2) ppar(j,i) = pevt5b(j,i,ievt)
                  p(j,i) = ppar(j,i)
               enddo
            enddo
            do i=1,5
               do j=1,5
                  if (ips.eq.1) y(i,j) = sevt5a(i,j,ievt)
                  if (ips.eq.2) y(i,j) = sevt5b(i,j,ievt)
               enddo
            enddo
            call ecuts(5,var,ipass)
            if (ipass.ne.0) then
               wtplot=wtps*wgt
               sig5w=sig5(i1,i2,i3,i4,i5,wtplot,var)
               sig5w=sig5w*wtps*var 
               sig4sw=sig4s(i1,i2,i3,i4,i5,wtplot) 
               sig4sw=sig4sw*wtps  
               if (njets.eq.3) then
                  sig3dsw=sig3ds(i1,i2,i3,i4,i5,wtplot)  
                  sig3dsw=sig3dsw*wtps 
               endif
               sig5w=sig5w-sig4sw-sig3dsw
            endif
         endif
         sig5asum = sig5asum + sig5w
      enddo
      sig5a = sig5asum
      return
      end
*************************************************************************
*
*************************************************************************
* sig3@NLO
*
* nloop = 0: three-parton tree-level matrix elements
* nloop = 1: three-parton 1-loop matrix elements
*************************************************************************

      function sig3(s12,s13,s23)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793238d0)
      parameter(zeta2=1.64493406684822644d0)
      parameter(zeta3=1.20205690315959429d0)
      common /qcd/as,ca,cflo,cf,tr,cn
      common /tcuts/ymin,y0
      common/inphys/nloop,icol,njets
      
      F30_ = 0d0
      G30_ = 0d0
      
      s123 = s12+s13+s23
      yb = 0.004765d0 ! H-bb~ coupling
      heft = -0.000055d0 ! H-gg coupling	
      nf = 5d0
      pscor3 = 0.5d0*s123
c     s123 = 1d0 so actually pscor3 = 0.5d0
	
      GgLO = (ca**2-1d0)*(heft)**2*125.18d0**3/128d0/pi
      BRgLO = GgLO/(GgLO + GbLO)

c     BRgLO = 0.258774d0
      
      if(nloop.eq.0)then
	
      if(icol.eq.0.or.icol.eq.1)then
        F30_ = F30(s12,s13,s23) + F30(s13,s12,s23)
      endif
      if(icol.eq.0.or.icol.eq.2)then
        G30_ = G30n(s12,s13,s23)
      endif
      
      sig3g = (as/2d0/pi)*(1d0/3d0*ca*F30_ + 2d0*nf*G30_)
	
      sig3 = pscor3*sig3g
	
      elseif(abs(nloop).eq.1)then

        s123 = s12+s13+s23
        yb = 0.004765d0 ! H-bb~ coupling
        heft = -0.000055d0 ! H-gg coupling
	nf = 5d0
	pscor3 = 0.5d0*s123
c	s123 = 1d0 so actually pscor3 = 0.5d0

 	
	GgLO = (ca**2-1d0)*heft**2*125.18d0**3/(128d0*pi)
        GgNLO = (1d0 + (alphas/(2d0*pi)*(73d0/6d0*ca - 7d0/3d0*nf)))*GgLO   
        BRgNLO = GgNLO/(GgNLO + GbNLO)
        
c	BRgNLO = 0.296153d0

c	Finite part of one-loop three-parton antennas
	F31 = 0d0
	F31f = 0d0
	G31l = 0d0
	G31sl = 0d0
	G31f = 0d0

c	Finite parts of subtraction terms
	F31s = 0d0 
	F31fs = 0d0 
	G31ls = 0d0
	G31sls = 0d0
	G31fs = 0d0
	
	F31fsE = 0d0	

* N^2 term
      if(icol.eq.0.or.icol.eq.1)then
        F31 = F31_(s12,s13,s23) + F31_(s13,s12,s23)
	F31s = 6d0*4d0*1d0/3d0*73d0/4d0*F30(s12,s13,s23)
      endif    
* N*NF term
      if(icol.eq.0.or.icol.eq.2)then
	F31f = F31hat(s12,s13,s23) + F31hat(s13,s12,s23)
	F31fs = 2d0*(-7d0/6d0)*F30(s12,s13,s23)
	
	G31l = G31(s12,s13,s23)
	G31ls = 2d0*34d0/3d0*G30n(s12,s13,s23)
	
	F31fsE = 2d0*(-1d0)*F30(s12,s13,s23)
      endif     
* -NF/N term
      if(icol.eq.0.or.icol.eq.3)then
	G31sl = G31tilde(s12,s13,s23)
	G31sls = 2d0*19d0/4d0*G30n(s12,s13,s23)
      endif
* NF^2 term
      if(icol.eq.0.or.icol.eq.4)then
	G31f = G31hat(s12,s13,s23)
	G31fs = 2d0*(-1d0)*G30n(s12,s13,s23)
      endif
      
      
        sig3gLO = (as/2d0/pi)*(1d0/3d0*ca*(F30(s12,s13,s23) + F30(s13,s12,s23)) + 2d0*nf*G30n(s12,s13,s23))

        sig3g = (as/2d0/pi)**2*(1d0/3d0*ca**2*F31 + ca*nf*(1d0/3d0*F31f + 2d0*G31l) - 2d0*nf/ca*G31sl + 2d0*nf**2*G31f)
         
        sig3r = (as/2d0/pi)**2*(1d0/12d0*ca**2*F31s + ca*nf*(F31fsE + G31ls) - nf/ca*G31sls + nf**2*G31fs)
        
        sig3 = pscor3*(sig3gLO + sig3g + sig3r)
     
     
      elseif(abs(nloop).eq.2)then
         write(6,*) 'NNLO not implemented'
         stop
         return   
      endif
      return
      end
*************************************************************************
*
*************************************************************************
* sig4@NLO
*
* nloop = 1: four-parton tree-level matrix elements
*************************************************************************
      function sig4(i1,i2,i3,i4,wtplot,var4)
      implicit real*8(a-h,l,o-z)
      parameter(pi=3.141592653589793238d0)
      parameter(zeta2=1.64493406684822644d0)
      parameter(zeta3=1.20205690315959429d0)
      parameter(zeta4=1.08232323371113819d0)
      logical iprint
      logical plot 
      common /plots/plot
      common/inphys/nloop,icol,njets
      common /qcd/as,ca,cflo,cf,tr,cn 
      common /tcuts/ymin,y0
      common /ps/n1,n2,n3,n4
      common /yij4/y(4,4)
      common /print/iprint

      s12=y(i1,i2)
      s13=y(i1,i3)
      s14=y(i1,i4)
      s23=y(i2,i3)
      s24=y(i2,i4)
      s34=y(i3,i4)
      
      s1234 = s12+s13+s14+s23+s24+s34

      yb = 0.004765d0 ! H-bb~ coupling
      heft = -0.000055d0 ! H-gg coupling
      nf = 5d0
      pscor4 = 0.25d0*s1234 
c     s1234 = 1d0 so actually pscor4 = 0.25d0
 	
      GgLO = (ca**2-1d0)*heft**2*125.18d0**3/128d0/pi
      GgNLO = (1d0 + (alphas/(2d0*pi)*(73d0/6d0*ca - 7d0/3d0*nf)))*GgLO
      BRgNLO = GgNLO/(GgNLO + GbNLO)
        
c      BRgNLO = 0.296153d0

      
      if(nloop.lt.0)one=0d0
      if(abs(nloop).eq.1)then

	F40l = 0d0
	F40lA = 0d0
	F40lB = 0d0
	F40lC = 0d0
	F40lD = 0d0
	F40lE = 0d0
	F40lF = 0d0

	G40fl = 0d0
	G40flA = 0d0
	G40flB = 0d0

	G40fsl = 0d0
	
	H40ff = 0d0		 

* N^2 term
	if(icol.eq.0.or.icol.eq.1)then
	   F40lA = F40(s12,s13,s14,s23,s24,s34) ! 1234
	   F40lB = F40(s12,s14,s13,s24,s23,s34) ! 1243
	   F40lC = F40(s13,s12,s14,s23,s34,s24) ! 1324
	   F40lD = F40(s13,s14,s12,s34,s23,s24) ! 1342
	   F40lE = F40(s14,s12,s13,s24,s34,s23) ! 1423
	   F40lF = F40(s14,s13,s12,s34,s24,s23) ! 1432
	
	   F40l = F40lA + F40lB + F40lC + F40lD + F40lE + F40lF
	endif
* NF*N term
	if(icol.eq.0.or.icol.eq.2)then
	   G40flA = G40(s12,s13,s14,s23,s24,s34) ! 1342
	   G40flB = G40(s12,s23,s24,s13,s14,s34) ! 2341
	  
	   G40fl = G40flA + G40flB
	endif
* - NF/N term
	if(icol.eq.0.or.icol.eq.3)then
	   G40fsl = G40tilde(s12,s13,s14,s23,s24,s34)
	endif
* NF^2 term
	if(icol.eq.0.or.icol.eq.4)then
	   H40ff = H40(s12,s13,s14,s23,s24,s34)
	endif

     	sig4g = (as/2d0/pi)**2*( 
     &		   		1d0/12d0*ca**2*F40l   
     &		   		+ nf*ca*G40fl    
     &		   		-nf/ca*G40fsl         
     &        	   		+ nf**2*H40ff  ) 

	sig4 = pscor4*sig4g
      
      elseif(abs(nloop).eq.2)then
         write(6,*) 'NNLO not implemented'
         stop
         return
      endif

      if(plot)then
         call bino(1,sig4*wtplot*var4,4)
      else
         call distrib(wtdis)
         sig4=sig4/wtdis
      endif

      return
      end
*************************************************************************
*
*************************************************************************
* sig3s@NLO
*
* nloop = 1: subtraction terms matrix elements at NLO
*************************************************************************

      function sig3s(i1,i2,i3,i4,wtplot)
c	pscor4 = 0.25d0*s1234
      implicit real*8(a-h,o-z)
      logical plot 
      parameter(pi=3.141592653589793238d0)
      common /qcd/as,ca,cflo,cf,tr,cn 
      common /tcuts/ymin,y0
      common /plots/plot
      common/inphys/nloop,icol,njets
      common /yij4/y(4,4)

      s12=y(i1,i2)
      s13=y(i1,i3)
      s14=y(i1,i4)
      s23=y(i2,i3)
      s24=y(i2,i4)
      s34=y(i3,i4)

      s1234 = s12+s13+s14+s23+s24+s34

      yb = 0.004765d0 ! H-bb~ coupling
      heft = -0.000055d0 ! H-gg coupling
      nf = 5d0
      pscor4 = 0.25d0*s1234 
c     s1234 = 1d0 so actually pscor4 = 0.25d0
 	
      GgLO = (ca**2-1d0)*heft**2*125.18d0**3/(128d0*pi)
      GgNLO = (1d0 + (alphas/(2d0*pi)*(73d0/6d0*ca - 7d0/3d0*nf)))*GgLO
      BRgNLO = GgNLO/(GgNLO + GbNLO)
        
c     BRgNLO = 0.296153d0      

      if(abs(nloop).eq.1)then 

	fac = -pscor4*(as/2d0/pi)**2
	wt1 = 1d0/12d0*ca**2*wtplot*fac
	wt2 = nf*ca*wtplot*fac
	wt3 = -nf/ca*wtplot*fac
	wt4 = nf**2*wtplot*fac

	sF40l = 0d0
	sG40fl = 0d0
	sG40fsl = 0d0
	sH40ff = 0d0	

* N^2 term
      if(icol.eq.0.or.icol.eq.1)then
	 sF40l = F40ls(i1,i2,i3,i4,wt1)
      endif
* NF*N term
      if(icol.eq.0.or.icol.eq.2)then
c	 sG40fl = G40fls(i1,i2,i3,i4,wt2)
	 sG40fl = G40flsE(i1,i2,i3,i4,wt2)
      endif
* -NF/n term
      if(icol.eq.0.or.icol.eq.3)then
        sG40fsl = G40fsls(i1,i2,i3,i4,wt3)
      endif
* NF^2 term
      if(icol.eq.0.or.icol.eq.4)then
        sH40ff = H40ffs(i1,i2,i3,i4,wt4)
      endif
 
        sig3gs = (as/2d0/pi)**2*( 
     .		    			1d0/12d0*ca**2*sF40l   
     .		    			+ nf*ca*sG40fl    
     .		    			-nf/ca*sG40fsl         
     .        	    			+ nf**2*sH40ff  ) 
 	
	sig3s = pscor4*sig3gs
    
      elseif(abs(nloop).eq.2)then
         write(6,*) 'NNLO not implemented'
         stop
         return
      endif
      
      return
      end
*
************************************************************************
*
      function sig5(i1,i2,i3,i4,i5,wtplot,var5)
      implicit real*8(a-h,o-z)
      logical plot 
      parameter(pi=3.141592653589793238d0)
      common /qcd/as,ca,cflo,cf,tr,cn 
      common/inphys/nloop,icol,njets
      common /tcuts/ymin,y0
      common /plots/plot
      common /yij5/y(5,5)

      write(6,*) 'NNLO not implemented'
      stop
      return

      if(plot)then
         call bino(1,sig5*wtplot*var5,5)
      else
         call distrib(wtdis)
         sig5=sig5/wtdis
      endif
      return
      end
*
************************************************************************
*
      function sig4s(i1,i2,i3,i4,i5,wtplot)
      implicit real*8(a-h,o-z)
      logical plot 
      common /qcd/as,ca,cflo,cf,tr,cn 
      common /tcuts/ymin,y0
      common /plots/plot
      common/inphys/nloop,icol,njets
      common /yij5/y(5,5)
      parameter(pi=3.141592653589793238d0)

      write(6,*) 'NNLO not implemented'
      stop
      return
      end
*
************************************************************************
*
      function sig3ds(i1,i2,i3,i4,i5,wtplot)
      implicit real*8(a-h,o-z)
      logical plot 
      parameter(pi=3.141592653589793238d0)
      common /qcd/as,ca,cflo,cf,tr,cn 
      common /tcuts/ymin,y0
      common /plots/plot
      common/inphys/nloop,icol,njets
      common /yij5/y(5,5)

      write(6,*) 'NNLO not implemented'
      stop
      return
      end
********************************************
