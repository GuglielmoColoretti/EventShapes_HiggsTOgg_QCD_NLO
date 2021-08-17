      program eerad3
      implicit real*8 (a-h,o-z)
      common/inphys/nloop,icol,njets
      common /ivegas/iwarm,iprod

      call readinit
      call writer 
      call cross(avgi,sd)
      if(nloop.eq.0)then
        write(*,11)avgi,sd 
      elseif(nloop.eq.1)then
        write(*,12)avgi,sd 
      elseif(nloop.eq.2)then
        write(*,13)avgi,sd 
      elseif(nloop.eq.-1)then
        write(*,14)avgi,sd 
      elseif(nloop.eq.-2)then 
        write(*,15)avgi,sd 
      endif	
      if (iprod.eq.1) call outfiles

   11 format(/,'        [LO]      ',g14.6,' +- ',g14.6,/)
   12 format(/,'       [NLO]      ',g14.6,' +- ',g14.6,/)
   13 format(/,'      [NNLO]      ',g14.6,' +- ',g14.6,/)
   14 format(/,'  only [NLO]      ',g14.6,' +- ',g14.6,/)
   15 format(/,' only [NNLO]      ',g14.6,' +- ',g14.6,/)

      stop
      end

      subroutine readinit
      implicit real*8(a-h,o-z)
      character*2 ianame(1:3)
      character*40 ibname(1:3)
      character*40 froot
      character*20 fname
      character*1 itag
      dimension iseeds(1:2,0:99),icolmax(0:2)
      parameter(pi=3.141592653589793238d0)
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common /tcuts/ymin,y0
      common /qcd/as,ca,cflo,cf,tr,cn 
      common/masses/rm2(1:5),shat
      common /cuts/ycutJ,ycutD,ycutG,Bcut,Ccut,Fcut,Tcut,Scut,em2hcut
      common/rseeds/i1,i2
      common /runinfo/itmax1,itmax2,nshot3,nshot4,nshot5(2) 
      common/outfile/fname
      common /ivegas/iwarm,iprod

      data iseeds/
     . 2991, 9281, 3823, 8151, 6239, 8754, 9514, 3746, 4893, 7731, 
     . 8635, 4662, 4132, 3764, 5234, 4892, 9753, 7629, 2531, 1638, 
     . 9376, 4157, 6472, 8086, 8526, 1264, 4723, 7549, 9472, 4637, 
     . 3746, 4682, 9854, 8472, 2332, 7724, 5524, 6352, 3652, 7382, 
     . 4431, 3667, 4538, 4888, 1838, 7655, 6764, 8382, 4885, 9980, 
     . 9876, 2312, 8743, 5521, 7512, 1123, 9814, 5531, 1249, 6413, 
     . 3468, 2814, 4392, 0921,  736, 2881, 4831, 8381, 2319, 4452, 
     . 3112,  978, 4791,  541, 1598, 8421, 5798, 4359, 2345, 3258, 
     . 3409,   91,  714, 8383, 2553, 7344, 9589, 4931,  941, 3881, 
     . 9145, 3885,   12, 4949, 9693, 3919, 9592, 1312, 6667, 2823, 
     . 5994, 2864, 4913, 6511, 8726, 2745, 4537, 2371, 6148, 4714, 
     . 8414,   28, 3851,  284, 3884, 2813,  847, 1822, 5871, 2442, 
     . 4771, 7411, 4728, 2771, 4881, 9831, 1236, 8831,  837, 3771, 
     . 5993, 1277, 3441, 8482,  693, 4882, 9391, 3899, 3991, 1134, 
     . 2884, 3851, 5918,  761, 9278, 9391,  101, 2391, 8124, 9391, 
     .  401, 4313, 3855, 8731, 4721, 3893,  145, 4992, 7573, 8471, 
     . 9981, 2383,  921, 3881,  388, 1384, 7628, 7237, 9492, 3884, 
     .   28, 8313, 3883, 8281, 1332, 3927, 5452, 7246, 2812, 2831, 
     . 7691, 3718, 4881,  331, 6266, 2163, 4772, 7371, 9931, 3881, 
     . 7713, 6361, 3488, 9281, 4913, 3813, 3463,   89,  918, 8381 / 

* iaver=0 : all distributions
* iaver=1 : BW
* iaver=2 : C
* iaver=3 : rho=MH^2/s
* iaver=4 : 1-T
* iaver=5 : BT
* iaver=6 : y23D
* iaver=7 : y23J
* iaver=8 : jet distributions
*
* nloop =  0: leading order only
* nloop =  1; leading order + next-to-leading order
* nloop = -1; next-to-leading order only
*
*
* icol =  0: all colours
* nloop = 0
* icol =  1:   N
* icol =  2:   NF
* nloop = +/- 1
* icol =  1:   N^2
* icol =  2:   NF*N
* icol =  3:   -NF/N
* icol =  4:   NF^2
      
      n = iargc()
      if (n.ge.2) then 
         call getarg(1,ianame(1))
         call getarg(2,ibname(1))
         nl = 1
      endif
      if (n.ge.4) then
         call getarg(3,ianame(2))
         call getarg(4,ibname(2))
         nl = 2
      endif
      inum = -1
      ifile = -1
      froot='eerad3.input'
      do i=1,nl 
         if (ianame(i).eq.'-n') then
            read(ibname(i),*) inum
         endif
         if (ianame(i).eq.'-i') then
            ifile=0
            froot = ibname(i) 
            ilen = len(froot)
         endif
      enddo
      if (inum.lt.0.or.inum.gt.99) inum = 0
      irlen = 0
      do i=1,ilen
         if (ichar(froot(i:i)).eq.32.and.irlen.eq.0) then
            irlen = i-1
         endif
      enddo
      i1 = iseeds(1,inum)
      i2 = iseeds(2,inum)
      if (inum.lt.10) then  
         fname='E0'
         write(fname(3:3),'(I1)') inum
      else
         fname='E'
         write(fname(2:3),'(I2)') inum
      endif
      open(9,file=froot)
      read(9,*) y0
      read(9,*) iaver
      read(9,*) cutvar
      read(9,*) imom
      read(9,*) iang
      read(9,*) nloop
      read(9,*) icol
      read(9,*) itag
      read(9,*) iwarm,iprod
      read(9,*) itmax1,itmax2
      read(9,*) nshot3,nshot4,nshot5a
      nshot5(1) = nshot5a
      nshot5(2) = nshot5a/5
      close(9)
      idist = 0
      if (y0.gt.1d-4) y0=1d-4
      if (y0.lt.1d-8) y0=1d-8
      iy1 = int(-log10(y0))
      iy0 = int(y0*10**(iy1))
      if (iy0.eq.0) then 
         iy1 = iy1+1
         iy0 = int(y0*10**(iy1))
      endif
      if (abs(nloop).gt.2) nloop=0
      icolmax(0) = 2
      icolmax(1) = 4
      if (icol.le.0.or.icol.gt.icolmax(abs(nloop))) icol=0
      fname = fname(1:3)//'.y'//char(iy0+48)//'d'//char(iy1+48)//'.'
      fname = fname(1:8)//'.i'//itag(1:1)//char(icol+48)//'.'

       ymin=1d0

**** the default cuts for the event shape run
      if (iaver.eq.0) then
         ycutJ=cutvar 
         ycutD=cutvar
         ycutG=cutvar 
         Bcut = cutvar
         BWcut=cutvar
         BTcut=cutvar
         Ccut=cutvar
         Fcut=cutvar
         Tcut=cutvar
         Scut=cutvar
         em2hcut=cutvar
      endif
**** the default cuts for the jet run
      if (iaver.eq.8) then
         ycutD=dexp(-10d0) 
         Bcut=ycutD
         BWcut=Bcut
         BTcut=Bcut
         Ccut=ycutD
         Fcut=ycutD
         Tcut=ycutD
         Scut=ycutD
         em2hcut=ycutD
      endif

      ome=0d0
      if(iaver.eq.2)Ccut=cutvar
      if(iaver.eq.4)Tcut=cutvar
      if(iaver.eq.3)em2hcut=cutvar
      if(iaver.eq.6)ycutD=cutvar
      if(iaver.eq.7)ycutJ=cutvar
      if(iaver.eq.1.or.iaver.eq.5) then
         Bcut=cutvar
         BWcut = Bcut
         BTcut = Bcut
      endif

      do i=1,5 
         rm2(i) = 0d0
      enddo
      shat = 1d0

      njets = 3
      idebug=0
      as=2d0*pi                
      cflo=4d0/3d0
      cf=4d0/3d0
      tr=0.5d0*5d0
      ca=3d0
      cn=ca

      return
      end


      subroutine writer
      implicit real*8 (a-h,o-z)
      character*60 slong,sno
      character*24 short
      character*9 sblank1
      character*24 sblank2
      character*34 sblank3
      character*60 starline
      character*1 star 
      character*20 fname
      common /qcd/as,ca,cflo,cf,tr,cn 
      common /tcuts/ymin,y0
      common /cuts/ycutJ,ycutD,ycutG,Bcut,Ccut,Fcut,Tcut,Scut,em2hcut
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common /runinfo/itmax1,itmax2,nshot3,nshot4,nshot5(2) 
      common /ivegas/iwarm,iprod
      common/outfile/fname
   
      starline=
     . '************************************************************'
      sblank1=' '
      sblank2=' '
      sblank3=' '
      sno=' '
      star='*'
      write(*,*)star,starline,star
      write(*,*)star,starline,star
      write(*,*)star,sno,star
      slong = ' EERAD3: event shapes and jet rates at O(alpha_s^3)'
      write(*,*)star,slong,star
      write(*,*)star,sno,star
      slong = ' A.Gehrmann-De Ridder, T.Gehrmann, N.Glover, G.Heinrich'
      write(*,*)star,slong,star
      write(*,*)star,sno,star
      slong = ' JHEP 0711 (2007) 058 [arXiv:0710.0346]'
      write(*,*)star,slong,star
      write(*,*)star,sno,star
      if(nloop.eq.0) slong= '    e+e- -> 3 jets at LO ' 
      if(nloop.eq.1) slong= '    e+e- -> 3 jets at LO + NLO ' 
      if(nloop.eq.2) slong= '    e+e- -> 3 jets at LO + NLO + NNLO' 
      if(nloop.eq.-1) slong= '    e+e- -> 3 jets at NLO only ' 
      if(nloop.eq.-2) slong= '    e+e- -> 3 jets at NNLO only ' 
      write(*,*)star,slong,star
      slong =''
      if(abs(nloop).eq.0)then
      if(icol.eq.0)slong='    all colours'
      if(icol.eq.1)slong='    colour factor:	N'
      if(icol.eq.2)slong='    colour factor:	NF'
      endif
      if(abs(nloop).eq.1)then
      if(icol.eq.0)slong='    all colours'
      if(icol.eq.1)slong='    colour factor:  N^2'
      if(icol.eq.2)slong='    colour factor:  NF*N'
      if(icol.eq.3)slong='    colour factor:  -NF/N'
      if(icol.eq.4)slong='    colour factor:  NF^2'
      endif

      write(*,*)star,slong,star
      write(*,*)star,sno,star
      write(*,*)star,starline,star

      write(*,*)star,sno,star
      slong=' input QCD parameters :'
      write(*,*)star,slong,star
      write(*,*)star,sno,star

c      short='    alpha_s(M_Z) = '
c      write(*,10)star,short,as,sblank2,star    
      short='             ca  = '
      write(*,10)star,short,ca,sblank2,star    
      short='             cf  = '
      write(*,10)star,short,cf,sblank2,star    
      short='             tr  = '
c      write(*,10)star,short,tr,sblank2,star    
c      short='             cn  = '
c      write(*,10)star,short,cn,sblank2,star    
c      short='             cflo  = '
c      write(*,10)star,short,cflo,sblank2,star    

      write(*,*)star,sno,star
      write(*,*)star,starline,star
      write(*,*)star,sno,star
      slong=' input technical parameters :'
      write(*,*)star,slong,star

      write(*,*)star,sno,star
      short='           ymin = '
      write(*,*)star,short,ymin,sblank1,star    
c      short='           y0   = '
c      write(*,*)star,short,y0,sblank1,star    

      write(*,*)star,sno,star
      write(*,*)star,starline,star
      write(*,*)star,sno,star
      if(iaver.eq.0)then
      slong = ' cuts applied : '
      write(*,*)star,slong,star
      write(*,*)star,sno,star

      short='  ycut in   Jade scheme ' 
      write(*,10)star,short,ycutJ,sblank2,star    
      short='  ycut in Durham scheme ' 
      write(*,10)star,short,ycutD,sblank2,star    
c      short='  ycut in Geneva scheme ' 
c      write(*,10)star,short,ycutG,sblank2,star    
      short='   BW,BT parameters     > ' 
      write(*,10)star,short,Bcut,sblank2,star    
      short='   C parameter        > ' 
      write(*,10)star,short,Ccut,sblank2,star    
c      short='   T major parameter  > ' 
c      write(*,10)star,short,Fcut,sblank2,star    
c      short='   S parameter        > ' 
c      write(*,10)star,short,Scut,sblank2,star    
      short='   rho=M_h^2/s  parameter   > ' 
      write(*,10)star,short,em2hcut,sblank2,star    
      short='   1-T parameter      > ' 
      write(*,10)star,short,Tcut,sblank2,star    
      elseif(iaver.eq.1)then
      slong='   weighted by BW       ' 
      write(*,*)star,slong,star
      short='   BW parameter     > ' 
      write(*,10)star,short,Bcut,sblank2,star    
      elseif(iaver.eq.2)then
      slong='   weighted by Cpar     ' 
      write(*,*)star,slong,star
      short='   C parameter       >  ' 
      write(*,10)star,short,Ccut,sblank2,star    
      elseif(iaver.eq.3)then
      slong='   weighted by rho=Mh^2/s  ' 
      write(*,*)star,slong,star
      short='  rho=M_h^2/s  parameter >   ' 
      write(*,10)star,short,em2hcut,sblank2,star    
      elseif(iaver.eq.4)then
      slong='   weighted by 1-T      ' 
      write(*,*)star,slong,star
      short='   1-T  parameter >     ' 
      write(*,10)star,short,Tcut,sblank2,star    
      elseif(iaver.eq.5)then
      slong='   weighted by BT        ' 
      write(*,*)star,slong,star
      short='   BT  parameter >       ' 
      write(*,10)star,short,Bcut,sblank2,star    
      elseif(iaver.eq.6)then
      slong='   weighted by y23D     ' 
      write(*,*)star,slong,star
      short='   ycut in Durham scheme' 
      write(*,10)star,short,ycutD,sblank2,star    
      elseif(iaver.eq.7)then
      slong='   weighted by y23J     ' 
      write(*,*)star,slong,star
      short='   ycut in  Jade scheme ' 
      write(*,10)star,short,ycutJ,sblank2,star    
      elseif(iaver.eq.8)then
      slong='   weighted by y23J     ' 
      write(*,*)star,slong,star
      slong = ' cuts applied : '
      write(*,*)star,slong,star
      write(*,*)star,sno,star
      short='   BW,BT parameters     > ' 
      write(*,10)star,short,Bcut,sblank2,star    
      short='   C parameter        > ' 
      write(*,10)star,short,Ccut,sblank2,star    
c      short='   T major parameter  > ' 
c      write(*,10)star,short,Fcut,sblank2,star    
      short='   1-T parameter      > ' 
      write(*,10)star,short,Tcut,sblank2,star    
      short='   rho=M_h^2/s  parameter   > ' 
      write(*,10)star,short,em2hcut,sblank2,star    
      short='   ycut in Durham scheme' 
      write(*,10)star,short,ycutD,sblank2,star    
      endif   
      write(*,*)star,sno,star
      write(*,*)star,starline,star
      write(*,*)star,sno,star
      slong=' Distributions:        '      
      write(*,*)star,slong,star
      if (iaver.eq.0) then
      slong='  BW, C, rho=MH^2/s, 1-T, BT '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.1) then
      slong='   BW   '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.2) then
      slong='   C   '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.3) then
      slong='   rho=MH^2/s  '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.4) then
      slong='   1-T '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.5) then
      slong='   BT  '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.6) then
      slong=' y23D  '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.8) then
      slong='   ln(1-T), ln(rho), ln(C), ln(BT), ln(BW)  '
      write(*,*)star,slong,star
      slong='   ln(y23), ln(y34),ln(y45), R_3, R_4, R_5 '
      write(*,*)star,slong,star
      endif
      if (iaver.eq.1.or.iaver.eq.2.or.iaver.eq.3.or.
     .    iaver.eq.4.or.iaver.eq.5.or.iaver.eq.6) then
      short=' Compute moment: '
      write(*,11) star,short,imom,sblank3,star
      endif
      write(*,*)star,sno,star
      write(*,*)star,starline,star
      write(*,*)star,sno,star
      short='  Output filenames: '
      write(*,12) star,short,fname(1:12),sblank2,star
      write(*,*)star,sno,star
      write(*,*)star,starline,star
      write(*,*)star,sno,star
      if (iwarm.eq.1) then
         short = ' vegas warmup steps:    '
         write(*,11)star,short,itmax1,sblank3,star
      else
         slong = ' vegas reads grid files'
         write(*,*)star,slong,star
      endif
      if (iprod.eq.1) then
         short = ' vegas production steps:'
         write(*,11)star,short,itmax2,sblank3,star
      else
         slong = 'grid run only, no production'
         write(*,*)star,slong,star
      endif
      write(*,*)star,sno,star
      write(*,*)star,starline,star
 10   format(1x,A,A,1pe12.4,A,A)
 11   format(1x,A,A,i2,A,A)
 12   format(1x,A,A,A,A,A)
      return
      end   
************************************************************************
*                                                                      *
* compute cross section (ave) with standard deviation (sd) from       *
* itmax2 iterations                                                    *
*                                                                      *
************************************************************************
      subroutine cross(ave,sd)
      implicit real*8(a-h,o-z)
      logical plot 
      COMMON /BVEG3/NDIM3,NCALL3,NPRN3
      COMMON /BVEG4/NDIM4(1),NCALL4(1),NPRN4
      COMMON /BVEG5/NDIM5(2),NCALL5(2),NPRN5
      common /plots/plot
      common /runinfo/itmax1,itmax2,nshot3,nshot4,nshot5(2) 
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common /phase/ips 
      common /ivegas/iwarm,iprod
      external sig3a,sig4a,sig5a  

      nprn3=1
      nprn4=1
      nprn5=1
      ndim3=2
      ndim4(1)=5
      ndim5(1)=8
      ndim5(2)=8

c---- number of shots for warmup run
      if (iprod.eq.1) then
         ncall3=nshot3/5d0
         ncall4(1)=nshot4/5d0
         do i=1,2
            ncall5(i)=nshot5(i)/5d0
         enddo
      else
         ncall3=nshot3
         ncall4(1)=nshot4
         do i=1,2
            ncall5(i)=nshot5(i)
         enddo
      endif
       plot=.false.

      init = 0
      if (iwarm.eq.1) then
         init = 1
      endif

*
* initialize vegas
* init=1 cold start   init=0 input grid from previous run
*
      call vegas3(init,sig3a,ave3,sd3,chi2)
        if(abs(nloop).ge.1)then
          do ips4=1,1
             call vegas4a(init,sig4a,ips4,ave4,sd4,chi2)
          enddo
       endif	
        if(abs(nloop).eq.2)then
          do ips=1,2
            if (icol.ne.6) call vegas5(init,sig5a,ips,ave5,sd5,chi2)
          enddo
        endif
*
* vegas sweeps with grid adjustments
*
      if(init.eq.1) then
         do it=1,itmax1
            call vegas3(3,sig3a,ave3,sd3,chi2)
            if(abs(nloop).ge.1)then
               rsum4 = 0d0
               sdsum4 = 0d0
               do ips4=1,1
                  call vegas4a(3,sig4a,ips4,ave4,sd4,chi2)
                  rsum4 = rsum4 + ave4
                  sdsum4 = sdsum4 + sd4**2
               enddo
               ave4 = rsum4
               sd4 = sqrt(sdsum4)
            endif	
            if(abs(nloop).eq.2)then
               rsum=0d0
               sdsum=0d0
            endif
            write(*,10)ave3,sd3  
            if(abs(nloop).ge.1)write(*,11)ave4,sd4 
            if(abs(nloop).ge.2)write(*,12)rsum,sqrt(sdsum) 
            write(*,13)ave3+ave4+rsum,sqrt(sd3**2+sd4**2+sdsum) 
         enddo
      endif
*
* main run
*
      if (iprod.eq.1) then
         ncall3=nshot3
         ncall4=nshot4
         do i=1,2
            ncall5(i)=nshot5(i)
         enddo
         plot=.true.
         call bino(0,0d0,0)
*
* initialize vegas
*
         call vegas3(2,sig3a,ave3,sd3,chi2)
         if(abs(nloop).ge.1)then
            do ips4=1,1
               call vegas4a(2,sig4a,ips4,ave4,sd4,chi2)
            enddo
         endif	
         if(abs(nloop).eq.2)then
            do ips=1,2
               if (icol.ne.6) call vegas5(2,sig5a,ips,ave5,sd5,chi2)
            enddo
         endif
*     
* vegas sweeps with frozen grid
*
         sum=0d0
         sum2=0d0
         do it=1,itmax2
            call vegas3(4,sig3a,ave3,sd3,chi2)
            if(it.eq.itmax2)then
               sum=sum+ave3
               sum2=sum2+sd3**2
            endif
            if(abs(nloop).ge.1)then
               rsum4 = 0d0
               sdsum4 = 0d0
               do ips4=1,1
                  call vegas4a(4,sig4a,ips4,ave4,sd4,chi2)
                  rsum4 = rsum4 + ave4
                  sdsum4 = sdsum4 + sd4**2
               enddo
               ave4 = rsum4
               sd4 = sqrt(sdsum4)
               if(it.eq.itmax2)then
                  sum=sum+ave4
                  sum2=sum2+sd4**2
               endif
            endif	
            if(abs(nloop).eq.2)then
               rsum=0d0
               sdsum=0d0
               if (icol.ne.6) then
                  do ips=1,2
                     call vegas5(4,sig5a,ips,ave5,sd5,chi2)
                     if(it.eq.itmax2)then
                        sum=sum+ave5
                        sum2=sum2+sd5**2
                     endif
                     rsum=rsum+ave5
                     sdsum=sdsum+sd5**2  
                  enddo
               endif            ! end if icol ne 6 
            endif
            call bino(2,0d0,0)
            if (it.lt.itmax2) then
               write(*,10)ave3,sd3  
               if(abs(nloop).ge.1)write(*,11)ave4,sd4 
               if(abs(nloop).ge.2)write(*,12)rsum,sqrt(sdsum) 
               write(*,13)ave3+ave4+rsum,sqrt(sd3**2+sd4**2+sdsum) 
            else
               write(*,14)ave3,sd3  
               if(abs(nloop).ge.1)write(*,15)ave4,sd4 
               if(abs(nloop).ge.2)write(*,16)rsum,sqrt(sdsum) 
               write(*,17)ave3+ave4+rsum,sqrt(sd3**2+sd4**2+sdsum) 
            endif
         enddo
         call bino(3,0d0,0)
      endif
 10   format(/,' sweep 3 parton    ',g14.6,' +- ',g14.6)
 11   format(/,' sweep 4 parton    ',g14.6,' +- ',g14.6)
 12   format(/,' sweep 5 parton    ',g14.6,' +- ',g14.6)
 13   format(/,' total result      ',g14.6,' +- ',g14.6,/)
 14   format(/,' final sweep 3 parton    ',g14.6,' +- ',g14.6)
 15   format(/,' final sweep 4 parton    ',g14.6,' +- ',g14.6)
 16   format(/,' final sweep 5 parton    ',g14.6,' +- ',g14.6)
 17   format(/,' final total result      ',g14.6,' +- ',g14.6,/)
*
      ave=sum  
      sd=sqrt(sum2)
      return
      end
