c---- the histogram management routines

************************************************************************
*                                                                      *
* write distributions                                                  *
*                                                                      *
************************************************************************
      subroutine bino(istat,wgt,npar)
      implicit real*8(a-h,o-z)
      common /pcut/ppar(4,5) 
      common/intech/iaver,imom,idist,iang,idebug
      common /jetdata/y45J,y34J,y23J,y45D,y34D,y23D,y45G,y34G,y23G
      common /cuts/ycutJ,ycutD,ycutG,Bcut,Ccut,Fcut,Tcut,Scut,em2hcut
      common /evdata/Cpar,Dpar,Spar,Apar,Planar,Tpar,
     .       Tmajor,Tminor,Opar,em2h,em2l,em2d,
     .       bmax,bmin,bsum,bdiff
      common /runinfo/itmax1,itmax2,nshot3,nshot4,nshot5(2) 
*
*
* init histograms
*
      nhis=100
      if (istat.eq.0) then
        ybinJ=0.005d0
        ybinD=0.005d0
        ybinG=0.01d0
        nbin=40

**** 1-T
        call histoi(11,0d0,0.5d0,200)
        call histoi(12,0d0,0.5d0,100)
        call histoi(13,0d0,0.5d0,50)
        call histoi(14,0d0,0.5d0,25)
        call histoi(15,0d0,0.5d0,200)
        call histoi(16,0d0,0.5d0,100)
        call histoi(17,0d0,0.5d0,50)
        call histoi(18,0d0,0.5d0,25)
**** C
        call histoi(21,0d0,1d0,400)
        call histoi(22,0d0,1d0,200)
        call histoi(23,0d0,1d0,100)
        call histoi(24,0d0,1d0,50)
        call histoi(25,0d0,1d0,400)
        call histoi(26,0d0,1d0,200)
        call histoi(27,0d0,1d0,100)
        call histoi(28,0d0,1d0,50)
**** rho=MH^2/s
        call histoi(31,0d0,0.5d0,200)
        call histoi(32,0d0,0.5d0,100)
        call histoi(33,0d0,0.5d0,50)
        call histoi(34,0d0,0.5d0,25)
        call histoi(35,0d0,0.5d0,200)
        call histoi(36,0d0,0.5d0,100)
        call histoi(37,0d0,0.5d0,50)
        call histoi(38,0d0,0.5d0,25)
**** BW
        call histoi(41,0d0,0.5d0,200)
        call histoi(42,0d0,0.5d0,100)
        call histoi(43,0d0,0.5d0,50)
        call histoi(44,0d0,0.5d0,25)
        call histoi(45,0d0,0.5d0,200)
        call histoi(46,0d0,0.5d0,100)
        call histoi(47,0d0,0.5d0,50)
        call histoi(48,0d0,0.5d0,25)
**** BT
        call histoi(1,0d0,0.5d0,200)
        call histoi(2,0d0,0.5d0,100)
        call histoi(3,0d0,0.5d0,50)
        call histoi(4,0d0,0.5d0,25)
        call histoi(5,0d0,0.5d0,200)
        call histoi(6,0d0,0.5d0,100)
        call histoi(7,0d0,0.5d0,50)
        call histoi(8,0d0,0.5d0,25)
**** -Log(y23D)
        call histoi(61,0d0,10d0,100)
        call histoi(62,0d0,10d0,50)
        call histoi(63,0d0,10d0,25)
**** -Log(y34D)
        call histoi(64,0d0,10d0,100)
        call histoi(65,0d0,10d0,50)
        call histoi(66,0d0,10d0,25)
**** -Log(y45D)
        call histoi(67,0d0,10d0,100)
        call histoi(68,0d0,10d0,50)
        call histoi(69,0d0,10d0,25)
**** sig3j (-Log(ycutD))
        call histoi(71,0d0,10d0,100)
        call histoi(72,0d0,10d0,50)
        call histoi(73,0d0,10d0,25)
**** sig4j (-Log(ycutD))
        call histoi(74,0d0,10d0,100)
        call histoi(75,0d0,10d0,50)
        call histoi(76,0d0,10d0,25)
**** sig5j (-Log(ycutD))
        call histoi(77,0d0,10d0,100)
        call histoi(78,0d0,10d0,50)
        call histoi(79,0d0,10d0,25)
**** -Log(1-T)
        call histoi(81,0d0,10d0,100)
        call histoi(82,0d0,10d0,50)
        call histoi(83,0d0,10d0,25)
**** -Log(C)
        call histoi(84,0d0,10d0,100)
        call histoi(85,0d0,10d0,50)
        call histoi(86,0d0,10d0,25)
**** -Log(BT)
        call histoi(87,0d0,10d0,100)
        call histoi(88,0d0,10d0,50)
        call histoi(89,0d0,10d0,25)
**** -Log(BW)
        call histoi(90,0d0,10d0,100)
        call histoi(91,0d0,10d0,50)
        call histoi(92,0d0,10d0,25)
**** -Log(rho)
        call histoi(93,0d0,10d0,100)
        call histoi(94,0d0,10d0,50)
        call histoi(95,0d0,10d0,25)
      endif
*
* write event into histogram
*
      if (istat.eq.1) then
      call getvar(var)
      wt=wgt/var

      if(1d0-Tpar.gt.Tcut)then				  
         dlt = -dlog(1d0-Tpar)
         call histoa(11,1d0-Tpar,wt*(1d0-Tpar))		 
         call histoa(12,1d0-Tpar,wt*(1d0-Tpar))		 
         call histoa(13,1d0-Tpar,wt*(1d0-Tpar))		 
         call histoa(14,1d0-Tpar,wt*(1d0-Tpar))
         call histoa(15,1d0-Tpar,wt)		 
         call histoa(16,1d0-Tpar,wt)		 
         call histoa(17,1d0-Tpar,wt)		 
         call histoa(18,1d0-Tpar,wt)
         call histoa(81,dlt,wt)
         call histoa(82,dlt,wt)
         call histoa(83,dlt,wt)
      endif

      if(Cpar.gt.Ccut)then
         dlc = -dlog(Cpar)
         call histoa(21,Cpar,wt*Cpar)		 
         call histoa(22,Cpar,wt*Cpar)		 
         call histoa(23,Cpar,wt*Cpar)		 
         call histoa(24,Cpar,wt*Cpar)
         call histoa(25,Cpar,wt)		 
         call histoa(26,Cpar,wt)		 
         call histoa(27,Cpar,wt)		 
         call histoa(28,Cpar,wt)
         call histoa(84,dlc,wt)
         call histoa(85,dlc,wt)
         call histoa(86,dlc,wt)
      endif

      if(em2h.gt.em2hcut)then
         dlm = -dlog(em2h)
         call histoa(31,em2h,wt*em2h)		 
         call histoa(32,em2h,wt*em2h)		 
         call histoa(33,em2h,wt*em2h)		 
         call histoa(34,em2h,wt*em2h)
         call histoa(35,em2h,wt)		 
         call histoa(36,em2h,wt)		 
         call histoa(37,em2h,wt)		 
         call histoa(38,em2h,wt)
         call histoa(93,dlm,wt)
         call histoa(94,dlm,wt)
         call histoa(95,dlm,wt)
      endif

      if(Bmax.gt.Bcut)then
         dlw = -dlog(Bmax)
         call histoa(41,Bmax,wt*Bmax)		 
         call histoa(42,Bmax,wt*Bmax)		 
         call histoa(43,Bmax,wt*Bmax)		 
         call histoa(44,Bmax,wt*Bmax)
         call histoa(45,Bmax,wt)		 
         call histoa(46,Bmax,wt)		 
         call histoa(47,Bmax,wt)		 
         call histoa(48,Bmax,wt)
         call histoa(90,dlw,wt)
         call histoa(91,dlw,wt)
         call histoa(92,dlw,wt)
      endif

      if(Bsum.gt.Bcut)then
         dlb = -dlog(Bsum)
         call histoa(1,Bsum,wt*Bsum)		 
         call histoa(2,Bsum,wt*Bsum)		 
         call histoa(3,Bsum,wt*Bsum)		 
         call histoa(4,Bsum,wt*Bsum)
         call histoa(5,Bsum,wt)		 
         call histoa(6,Bsum,wt)		 
         call histoa(7,Bsum,wt)		 
         call histoa(8,Bsum,wt)
         call histoa(87,dlb,wt)
         call histoa(88,dlb,wt)
         call histoa(89,dlb,wt)
      endif

      dly23 = 10d0
      dly34 = 10d0
      dly45 = 10d0
      if(y23D.gt.ycutD) then
         dly23 = -dlog(y23D)
         call histoa(61,dly23,wt)
         call histoa(62,dly23,wt)
         call histoa(63,dly23,wt)
      endif
      if(y34D.gt.ycutD) then
         dly34 = -dlog(y34D)
         call histoa(64,dly34,wt)
         call histoa(65,dly34,wt)
         call histoa(66,dly34,wt)
      endif
      if(y45D.gt.ycutD) then
         dly45 = -dlog(y45D)
         call histoa(67,dly45,wt)
         call histoa(68,dly45,wt)
         call histoa(69,dly45,wt)
      endif
      nbina = 100
      ybina = 10d0/real(nbina)
      nbinb = 50
      ybinb = 10d0/real(nbinb)
      nbinc = 25
      ybinc = 10d0/real(nbinc)


      do iya = 0,nbina
         dlya = real(iya)*ybina
         if (dlya.gt.dly23) then
            if (dlya.gt.dly34) then
              if (dlya.gt.dly45) call histoa(77,dlya+ybina/2d0,wt*ybina)
              if (dlya.le.dly45) call histoa(74,dlya+ybina/2d0,wt*ybina)
           else
              call histoa(71,dlya+ybina/2d0,wt*ybina)
           endif
        endif
      enddo
      do iyb = 0,nbinb
         dlyb = real(iyb)*ybinb
         if (dlyb.gt.dly23) then
            if (dlyb.gt.dly34) then
              if (dlyb.gt.dly45) call histoa(78,dlyb+ybinb/2d0,wt*ybinb)
              if (dlyb.le.dly45) call histoa(75,dlyb+ybinb/2d0,wt*ybinb)
           else
              call histoa(72,dlyb+ybinb/2d0,wt*ybinb)
           endif
        endif
      enddo
      do iyc = 0,nbinc
         dlyc = real(iyc)*ybinc
         if (dlyc.gt.dly23) then
            if (dlyc.gt.dly34) then
              if (dlyc.gt.dly45) call histoa(79,dlyc+ybinc/2d0,wt*ybinc)
              if (dlyc.le.dly45) call histoa(76,dlyc+ybinc/2d0,wt*ybinc)
           else
              call histoa(73,dlyc+ybinc/2d0,wt*ybinc)
           endif
        endif
      enddo

      endif
*
* output distributions
*
      if (istat.eq.4) then
       if(iaver.eq.0.or.iaver.eq.4)then
         write(*,*)
         write(*,*) ' (1-T) distribution (1-T)/sig dsig/d(1-T)'
         write(*,*)
         call histow(11)
         call histowf(11,11)
         write(*,*)
         write(*,*) ' (1-T) distribution (1-T)/sig dsig/d(1-T)'
         write(*,*)
         call histow(12)
         call histowf(12,12)
         write(*,*)
         write(*,*) ' (1-T) distribution (1-T)/sig dsig/d(1-T)'
         write(*,*)
         call histow(13)
         call histowf(13,13)
         write(*,*)
         write(*,*) ' (1-T) distribution (1-T)/sig dsig/d(1-T)'
         write(*,*)
         call histow(14)
         call histowf(14,14)
         write(*,*)
         write(*,*) ' (1-T) distribution 1/sig dsig/d(1-T)'
         write(*,*)
         call histow(15)
         call histowf(15,15)
         write(*,*)
         write(*,*) ' (1-T) distribution 1/sig dsig/d(1-T)'
         write(*,*)
         call histow(16)
         call histowf(16,16)
         write(*,*)
         write(*,*) ' (1-T) distribution 1/sig dsig/d(1-T)'
         write(*,*)
         call histow(17)
         call histowf(17,17)
         write(*,*)
         write(*,*) ' (1-T) distribution 1/sig dsig/d(1-T)'
         write(*,*)
         call histow(18)
         call histowf(18,18)
      endif

      if(iaver.eq.0.or.iaver.eq.2)then
         write(*,*)
         write(*,*) ' C distribution C/sig dsig/dC'
         write(*,*)
         call histow(21)
         call histowf(21,21)
         write(*,*)
         write(*,*) ' C distribution C/sig dsig/dC'
         write(*,*)
         call histow(22)
         call histowf(22,22)
         write(*,*)
         write(*,*) ' C distribution C/sig dsig/dC'
         write(*,*)
         call histow(23)
         call histowf(23,23)
         write(*,*)
         write(*,*) ' C distribution C/sig dsig/dC'
         write(*,*)
         call histow(24)
         call histowf(24,24)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/dC'
         write(*,*)
         call histow(25)
         call histowf(25,25)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/dC'
         write(*,*)
         call histow(26)
         call histowf(26,26)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/dC'
         write(*,*)
         call histow(27)
         call histowf(27,27)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/dC'
         write(*,*)
         call histow(28)
         call histowf(28,28)
      endif

      if(iaver.eq.0.or.iaver.eq.3)then
         write(*,*)
         write(*,*) 'rho(=M_H^2/s) distribution rho/sig dsig/drho'
         write(*,*)
         call histow(31)
         call histowf(31,31)
         write(*,*)
         write(*,*) 'rho distribution rho/sig dsig/drho'
         write(*,*)
         call histow(32)
         call histowf(32,32)
         write(*,*)
         write(*,*) 'rho distribution rho/sig dsig/drho'
         write(*,*)
         call histow(33)
         call histowf(33,33)
         write(*,*)
         write(*,*) 'rho distribution rho/sig dsig/drho'
         write(*,*)
         call histow(34)
         call histowf(34,34)
         write(*,*)
         write(*,*) 'rho distribution 1/sig dsig/drho'
         write(*,*)
         call histow(35)
         call histowf(35,35)
         write(*,*)
         write(*,*) 'rho distribution 1/sig dsig/drho'
         write(*,*)
         call histow(36)
         call histowf(36,36)
         write(*,*)
         write(*,*) 'rho distribution 1/sig dsig/drho'
         write(*,*)
         call histow(37)
         call histowf(37,37)
         write(*,*)
         write(*,*) 'rho distribution 1/sig dsig/drho'
         write(*,*)
         call histow(38)
         call histowf(38,38)
      endif

      if(iaver.eq.0.or.iaver.eq.1)then
         write(*,*)
         write(*,*) ' BW distribution BW/sig dsig/dBW'
         write(*,*)
         call histow(41)
         call histowf(41,41)
         write(*,*)
         write(*,*) ' BW distribution BW/sig dsig/dBW'
         write(*,*)
         call histow(42)
         call histowf(42,42)
         write(*,*)
         write(*,*) ' BW distribution BW/sig dsig/dBW'
         write(*,*)
         call histow(43)
         call histowf(43,43)
         write(*,*)
         write(*,*) ' BW distribution BW/sig dsig/dBW'
         write(*,*)
         call histow(44)
         call histowf(44,44)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/dBW'
         write(*,*)
         call histow(45)
         call histowf(45,45)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/dBW'
         write(*,*)
         call histow(46)
         call histowf(46,46)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/dBW'
         write(*,*)
         call histow(47)
         call histowf(47,47)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/dBW'
         write(*,*)
         call histow(48)
         call histowf(48,48)
      endif

      if(iaver.eq.0.or.iaver.eq.5)then
         write(*,*)
         write(*,*) ' BT distribution BT/sig dsig/dBT'
         write(*,*)
         call histow(1)
         call histowf(1,51)
         write(*,*)
         write(*,*) ' BT distribution BT/sig dsig/dBT'
         write(*,*)
         call histow(2)
         call histowf(2,52)
         write(*,*)
         write(*,*) ' BT distribution BT/sig dsig/dBT'
         write(*,*)
         call histow(3)
         call histowf(3,53)
         write(*,*)
         write(*,*) ' BT distribution BT/sig dsig/dBT'
         write(*,*)
         call histow(4)
         call histowf(4,54)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/dBT'
         write(*,*)
         call histow(5)
         call histowf(5,55)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/dBT'
         write(*,*)
         call histow(6)
         call histowf(6,56)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/dBT'
         write(*,*)
         call histow(7)
         call histowf(7,57)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/dBT'
         write(*,*)
         call histow(8)
         call histowf(8,58)
      endif
      if(iaver.eq.8)then
         write(*,*)
         write(*,*) ' -Ln(y23D) distribution 1/sig dsig/d(-Ln(y23D))'
         write(*,*)
         call histow(61)
         call histowf(61,61)
         write(*,*)
         write(*,*) ' -Ln(y23D) distribution 1/sig dsig/d(-Ln(y23D))'
         write(*,*)
         call histow(62)
         call histowf(62,62)
         write(*,*)
         write(*,*) ' -Ln(y23D) distribution 1/sig dsig/d(-Ln(y23D))'
         write(*,*)
         call histow(63)
         call histowf(63,63)
         write(*,*)
         write(*,*) ' -Ln(y34D) distribution 1/sig dsig/d(-Ln(y34D))'
         write(*,*)
         call histow(64)
         call histowf(64,64)
         write(*,*)
         write(*,*) ' -Ln(y34D) distribution 1/sig dsig/d(-Ln(y34D))'
         write(*,*)
         call histow(65)
         call histowf(65,65)
         write(*,*)
         write(*,*) ' -Ln(y34D) distribution 1/sig dsig/d(-Ln(y34D))'
         write(*,*)
         call histow(66)
         call histowf(66,66)
         write(*,*)
         write(*,*) ' -Ln(y45D) distribution 1/sig dsig/d(-Ln(y45D))'
         write(*,*)
         call histow(67)
         call histowf(67,67)
         write(*,*)
         write(*,*) ' -Ln(y45D) distribution 1/sig dsig/d(-Ln(y45D))'
         write(*,*)
         call histow(68)
         call histowf(68,68)
         write(*,*)
         write(*,*) ' -Ln(y45D) distribution 1/sig dsig/d(-Ln(y45D))'
         write(*,*)
         call histow(69)
         call histowf(69,69)
         write(*,*)
         write(*,*) ' 3-jet ratio 1/s_tot s_3j(-Ln(ycutD))'
         write(*,*)
         call histow(71)
         call histowf(71,71)
         write(*,*)
         write(*,*) ' 3-jet ratio 1/s_tot s_3j(-Ln(ycutD))'
         write(*,*)
         call histow(72)
         call histowf(72,72)
         write(*,*)
         write(*,*) ' 3-jet ratio 1/s_tot s_3j(-Ln(ycutD))'
         write(*,*)
         call histow(73)
         call histowf(73,73)
         write(*,*)
         write(*,*) ' 4-jet ratio 1/s_tot s_4j(-Ln(ycutD))'
         write(*,*)
         call histow(74)
         call histowf(74,74)
         write(*,*)
         write(*,*) ' 4-jet ratio 1/s_tot s_4j(-Ln(ycutD))'
         write(*,*)
         call histow(75)
         call histowf(75,75)
         write(*,*)
         write(*,*) ' 4-jet ratio 1/s_tot s_4j(-Ln(ycutD))'
         write(*,*)
         call histow(76)
         call histowf(76,76)
         write(*,*)
         write(*,*) ' 5-jet ratio 1/s_tot s_5j(-Ln(ycutD))'
         write(*,*)
         call histow(77)
         call histowf(77,77)
         write(*,*)
         write(*,*) ' 5-jet ratio 1/s_tot s_5j(-Ln(ycutD))'
         write(*,*)
         call histow(78)
         call histowf(78,78)
         write(*,*)
         write(*,*) ' 5-jet ratio 1/s_tot s_5j(-Ln(ycutD))'
         write(*,*)
         call histow(79)
         call histowf(79,79)
         write(*,*)
         write(*,*) ' 1-T distribution 1/sig dsig/d(-Ln(1-T))'
         write(*,*)
         call histow(81)
         call histowf(81,81)
         write(*,*)
         write(*,*) ' 1-T distribution 1/sig dsig/d(-Ln(1-T))'
         write(*,*)
         call histow(82)
         call histowf(82,82)
         write(*,*)
         write(*,*) ' 1-T distribution 1/sig dsig/d(-Ln(1-T))'
         write(*,*)
         call histow(83)
         call histowf(83,83)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/d(-Ln(C))'
         write(*,*)
         call histow(84)
         call histowf(84,84)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/d(-Ln(C))'
         write(*,*)
         call histow(85)
         call histowf(85,85)
         write(*,*)
         write(*,*) ' C distribution 1/sig dsig/d(-Ln(C))'
         write(*,*)
         call histow(86)
         call histowf(86,86)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/d(-Ln(BT))'
         write(*,*)
         call histow(87)
         call histowf(87,87)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/d(-Ln(BT))'
         write(*,*)
         call histow(88)
         call histowf(88,88)
         write(*,*)
         write(*,*) ' BT distribution 1/sig dsig/d(-Ln(BT))'
         write(*,*)
         call histow(89)
         call histowf(89,89)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/d(-Ln(BW))'
         write(*,*)
         call histow(90)
         call histowf(90,90)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/d(-Ln(BW))'
         write(*,*)
         call histow(91)
         call histowf(91,91)
         write(*,*)
         write(*,*) ' BW distribution 1/sig dsig/d(-Ln(BW))'
         write(*,*)
         call histow(92)
         call histowf(92,92)
         write(*,*)
         write(*,*) ' rho distribution 1/sig dsig/d(-Ln(rho))'
         write(*,*)
         call histow(93)
         call histowf(93,93)
         write(*,*)
         write(*,*) ' rho distribution 1/sig dsig/d(-Ln(rho))'
         write(*,*)
         call histow(94)
         call histowf(94,94)
         write(*,*)
         write(*,*) ' rho distribution 1/sig dsig/d(-Ln(rho))'
         write(*,*)
         call histow(95)
         call histowf(95,95)
      endif
      endif
*
* event errors manipulation request, pipe through
*
      if ((istat.eq.2).or.(istat.eq.3)) then
         do i=1,nhis
            call histoe(istat,i)
         enddo
      endif
*
      return
      end
       
      subroutine outfiles
      implicit real*8(a-h,o-z)
      character*20 fnames(1:100),fname
      common/intech/iaver,imom,idist,iang,idebug
      common/outfile/fname

      fnames(11) = fname(1:13)//'T1a'
      fnames(12) = fname(1:13)//'T1b'
      fnames(13) = fname(1:13)//'T1c'
      fnames(14) = fname(1:13)//'T1d'
      fnames(15) = fname(1:13)//'T2a'
      fnames(16) = fname(1:13)//'T2b'
      fnames(17) = fname(1:13)//'T2c'
      fnames(18) = fname(1:13)//'T2d'
      
      fnames(21) = fname(1:13)//'C1a'
      fnames(22) = fname(1:13)//'C1b'
      fnames(23) = fname(1:13)//'C1c'
      fnames(24) = fname(1:13)//'C1d'
      fnames(25) = fname(1:13)//'C2a'
      fnames(26) = fname(1:13)//'C2b'
      fnames(27) = fname(1:13)//'C2c'
      fnames(28) = fname(1:13)//'C2d'
      
      fnames(31) = fname(1:13)//'M1a'
      fnames(32) = fname(1:13)//'M1b'
      fnames(33) = fname(1:13)//'M1c'
      fnames(34) = fname(1:13)//'M1d'
      fnames(35) = fname(1:13)//'M2a'
      fnames(36) = fname(1:13)//'M2b'
      fnames(37) = fname(1:13)//'M2c'
      fnames(38) = fname(1:13)//'M2d'
      
      fnames(41) = fname(1:13)//'W1a'
      fnames(42) = fname(1:13)//'W1b'
      fnames(43) = fname(1:13)//'W1c'
      fnames(44) = fname(1:13)//'W1d'
      fnames(45) = fname(1:13)//'W2a'
      fnames(46) = fname(1:13)//'W2b'
      fnames(47) = fname(1:13)//'W2c'
      fnames(48) = fname(1:13)//'W2d'
      
      fnames(51) = fname(1:13)//'B1a'
      fnames(52) = fname(1:13)//'B1b'
      fnames(53) = fname(1:13)//'B1c'
      fnames(54) = fname(1:13)//'B1d'
      fnames(55) = fname(1:13)//'B2a'
      fnames(56) = fname(1:13)//'B2b'
      fnames(57) = fname(1:13)//'B2c'
      fnames(58) = fname(1:13)//'B2d'
      
      fnames(61) = fname(1:13)//'Y3a'
      fnames(62) = fname(1:13)//'Y3b'
      fnames(63) = fname(1:13)//'Y3c'
      fnames(64) = fname(1:13)//'Y4a'
      fnames(65) = fname(1:13)//'Y4b'
      fnames(66) = fname(1:13)//'Y4c'
      fnames(67) = fname(1:13)//'Y5a'
      fnames(68) = fname(1:13)//'Y5b'
      fnames(69) = fname(1:13)//'Y5c'
      
      fnames(71) = fname(1:13)//'S3a'
      fnames(72) = fname(1:13)//'S3b'
      fnames(73) = fname(1:13)//'S3c'
      fnames(74) = fname(1:13)//'S4a'
      fnames(75) = fname(1:13)//'S4b'
      fnames(76) = fname(1:13)//'S4c'
      fnames(77) = fname(1:13)//'S5a'
      fnames(78) = fname(1:13)//'S5b'
      fnames(79) = fname(1:13)//'S5c'
      
      fnames(81) = fname(1:13)//'TLa'
      fnames(82) = fname(1:13)//'TLb'
      fnames(83) = fname(1:13)//'TLc'
      fnames(84) = fname(1:13)//'CLa'
      fnames(85) = fname(1:13)//'CLb'
      fnames(86) = fname(1:13)//'CLc'
      fnames(87) = fname(1:13)//'BLa'
      fnames(88) = fname(1:13)//'BLb'
      fnames(89) = fname(1:13)//'BLc'
      fnames(90) = fname(1:13)//'WLa'
      fnames(91) = fname(1:13)//'WLb'
      fnames(92) = fname(1:13)//'WLc'
      fnames(93) = fname(1:13)//'MLa'
      fnames(94) = fname(1:13)//'MLb'
      fnames(95) = fname(1:13)//'MLc'

      if (iaver.eq.0.or.iaver.eq.4) then
         open(11,file=fnames(11))
         open(12,file=fnames(12))
         open(13,file=fnames(13))
         open(14,file=fnames(14))
         open(15,file=fnames(15))
         open(16,file=fnames(16))
         open(17,file=fnames(17))
         open(18,file=fnames(18))
      endif
      if (iaver.eq.0.or.iaver.eq.2) then
         open(21,file=fnames(21))
         open(22,file=fnames(22))
         open(23,file=fnames(23))
         open(24,file=fnames(24))
         open(25,file=fnames(25))
         open(26,file=fnames(26))
         open(27,file=fnames(27))
         open(28,file=fnames(28))
      endif
      if (iaver.eq.0.or.iaver.eq.3) then
         open(31,file=fnames(31))
         open(32,file=fnames(32))
         open(33,file=fnames(33))
         open(34,file=fnames(34))
         open(35,file=fnames(35))
         open(36,file=fnames(36))
         open(37,file=fnames(37))
         open(38,file=fnames(38))
      endif
      if (iaver.eq.0.or.iaver.eq.1) then
         open(41,file=fnames(41))
         open(42,file=fnames(42))
         open(43,file=fnames(43))
         open(44,file=fnames(44))
         open(45,file=fnames(45))
         open(46,file=fnames(46))
         open(47,file=fnames(47))
         open(48,file=fnames(48))
      endif
      if (iaver.eq.0.or.iaver.eq.5) then
         open(51,file=fnames(51))
         open(52,file=fnames(52))
         open(53,file=fnames(53))
         open(54,file=fnames(54))
         open(55,file=fnames(55))
         open(56,file=fnames(56))
         open(57,file=fnames(57))
         open(58,file=fnames(58))
      endif
      if (iaver.eq.8) then
         open(61,file=fnames(61))
         open(62,file=fnames(62))
         open(63,file=fnames(63))
         open(64,file=fnames(64))
         open(65,file=fnames(65))
         open(66,file=fnames(66))
         open(67,file=fnames(67))
         open(68,file=fnames(68))
         open(69,file=fnames(69))
         open(71,file=fnames(71))
         open(72,file=fnames(72))
         open(73,file=fnames(73))
         open(74,file=fnames(74))
         open(75,file=fnames(75))
         open(76,file=fnames(76))
         open(77,file=fnames(77))
         open(78,file=fnames(78))
         open(79,file=fnames(79))
         open(81,file=fnames(81))
         open(82,file=fnames(82))
         open(83,file=fnames(83))
         open(84,file=fnames(84))
         open(85,file=fnames(85))
         open(86,file=fnames(86))
         open(87,file=fnames(87))
         open(88,file=fnames(88))
         open(89,file=fnames(89))
         open(90,file=fnames(90))
         open(91,file=fnames(91))
         open(92,file=fnames(92))
         open(93,file=fnames(93))
         open(94,file=fnames(94))
         open(95,file=fnames(95))
      endif

      call bino(4,0d0,0)
      do i=1,8
         close(i+10)
         close(i+20)
         close(i+30)
         close(i+40)
         close(i+50)
         close(i+60)
         close(i+70)
         close(i+80)
      enddo
      close(69)
      close(79)
      close(89)
      close(90)
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      return
      end
