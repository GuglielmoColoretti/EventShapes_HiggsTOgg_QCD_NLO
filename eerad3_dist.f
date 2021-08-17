      program eerad3_dist
      implicit real*8(a-h,o-z)
      character*2 ianame
      character*40 ibname
      character*40 ifile
      character*12 froot(11)
      dimension ibindata(10),ibindata2(10),ibindataC(10)
      character*2 fextdata(10),fextdata2(10)
      common/qcdpar/rmuarr(-10:10),dlarr(-10:10),alsarr(-10:10,1:3)

      n = iargc()
      ifile = 'eerad3_dist.input'
      if (n.eq.2) then
         call getarg(1,ianame)
         call getarg(2,ibname)
         if (ianame.eq.'-i') ifile = ibname
      endif 
 
      open(9,file=ifile)
      read(9,*) iaver
      read(9,*) froot(1)
      read(9,*) froot(2)
      read(9,*) froot(3)
      read(9,*) froot(4)
      read(9,*) froot(5)
      read(9,*) froot(6)
      read(9,*) froot(7)
      read(9,*) froot(8)
      read(9,*) froot(9)
      read(9,*) froot(10)
      read(9,*) froot(11)
      read(9,*) asmz
      read(9,*) rmz
      read(9,*) roots
      read(9,*) xmu
      close(9) 


      write(6,*) 'eerad3_dist'
      write(6,*) 
      write(6,*) 'NNLO distributions from single colour factors:'
      if (iaver.eq.0) write(6,*) '  BW, C, rho=MH^2/s, 1-T, BT '
      if (iaver.eq.1) write(6,*) '  BW  '
      if (iaver.eq.2) write(6,*) '  C  '
      if (iaver.eq.3) write(6,*) '  rho = MH^2/s '
      if (iaver.eq.4) write(6,*) '  1-T '
      if (iaver.eq.5) write(6,*) '  BT '
      if (iaver.eq.6) write(6,*) '  ln(y23D) '
      if (iaver.eq.6) write(7,*) '  ln(y23J) '
      if (iaver.eq.8) 
     .          write(6,*) '  ln(1-T), ln(rho), ln(C), ln(BT), ln(BW)  '
      if (iaver.eq.8) 
     .          write(6,*) '  ln(y23), ln(y34),ln(y45), R_3, R_4, R_5 '

      write(6,*)
      write(6,101) roots
      write(6,102) asmz
      write(6,103) xmu
      write(6,*)
      
 101  format('sqrt(s)       = ',f11.5)
 102  format('alphas(MZ)    = ',f11.5)
 103  format('default scale = ',f11.5,' * roots')

      call initqcd(asmz,rmz,roots,xmu)
      
      if (iaver.le.5) then
         idata=8
         fextdata(1) = '1a'
         fextdata(2) = '1b'
         fextdata(3) = '1c'
         fextdata(4) = '1d'
         fextdata(5) = '2a'
         fextdata(6) = '2b'
         fextdata(7) = '2c'
         fextdata(8) = '2d'
         ibindata(1) = 200
         ibindata(2) = 100
         ibindata(3) = 50
         ibindata(4) = 25
         ibindata(5) = 200
         ibindata(6) = 100
         ibindata(7) = 50
         ibindata(8) = 25
      elseif(iaver.eq.8) then
         idata=3
         ibindata(1)= 100
         ibindata(2)= 50
         ibindata(3)= 25
         fextdata(1) = 'La'
         fextdata(2) = 'Lb'
         fextdata(3) = 'Lc'
      endif
      if (iaver.ge.6.and.iaver.le.8) then
         idata2 = 9
         ibindata2(1) = 100
         ibindata2(2) =  50
         ibindata2(3) =  25
         ibindata2(4) = 100
         ibindata2(5) =  50
         ibindata2(6) =  25
         ibindata2(7) = 100
         ibindata2(8) =  50
         ibindata2(9) =  25
         fextdata2(1) = '3a'
         fextdata2(2) = '3b'
         fextdata2(3) = '3c'
         fextdata2(4) = '4a'
         fextdata2(5) = '4b'
         fextdata2(6) = '4c'
         fextdata2(7) = '5a'
         fextdata2(8) = '5b'
         fextdata2(9) = '5c'
      endif


* iaver=0 : all distributions
* iaver=1 : BW
* iaver=2 : C
* iaver=3 : MH
* iaver=4 : 1-T
* iaver=5 : BT
* iaver=6 : y23D
* iaver=7 : y23J
* iaver=8 : jet distributions
*
         
      if (iaver.eq.0.or.iaver.eq.1) then
         call combinedist('W',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
      endif
      if (iaver.eq.0.or.iaver.eq.2) then
         do i=1,idata
            ibindataC(i) = 2*ibindata(i)
         enddo
         call combinedist('C',froot,idata,ibindataC(1:idata)
     .                       ,fextdata(1:idata))
      endif
      if (iaver.eq.0.or.iaver.eq.3) then
         call combinedist('M',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
      endif
      if (iaver.eq.0.or.iaver.eq.4) then
         call combinedist('T',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
      endif
      if (iaver.eq.0.or.iaver.eq.5) then
         call combinedist('B',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
      endif

      if (iaver.eq.8) then
         call combinedist('W',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
         call combinedist('C',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
         call combinedist('M',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
         call combinedist('T',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
         call combinedist('B',froot,idata,ibindata(1:idata)
     .                       ,fextdata(1:idata))
      endif

      if (iaver.ge.6.and.iaver.le.8) then
         call combinedist('Y',froot,idata2,ibindata2(1:idata2)
     .                       ,fextdata2(1:idata2))
         call combinedist('S',froot,idata2,ibindata2(1:idata2)
     .                       ,fextdata2(1:idata2))
      endif

      stop
      end


      subroutine initqcd(asmz,rmz,roots,xmu)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793238d0)
      common/lqcd/rlambda,rnf
      common/qcdpar/rmuarr(-10:10),dlarr(-10:10),alsarr(-10:10,1:3)

      rnf = 5d0
      do imu = -10,10,1
         rmuarr(imu) = xmu*roots*(2d0**(real(imu)/10d0))
         dlarr(imu) = 2d0*dlog(rmuarr(imu)/roots)
      enddo

      do iord = 1,3
         call alini(asmz,rmz,iord)
         do imu = -10,10,1
            alsarr(imu,iord) = alphas(rmuarr(imu),iord)/2d0/pi
         enddo
      enddo
      
      return
      end

      subroutine combinedist(ichar,froot,idata,ibindata,fextdata)
      implicit real*8(a-h,o-z)
      character*1 ichar
      character*12 froot(11)
      dimension ibindata(idata)
      character*2 fextdata(idata),fext
      parameter(nbins=400)
      dimension sig0(1:nbins),err0(1:nbins),sig1(1:nbins),err1(1:nbins),
     .          sig2(1:nbins),err2(1:nbins)
      common/qcdpar/rmuarr(-10:10),dlarr(-10:10),alsarr(-10:10,1:3)
      common/histo/t(1:nbins),siga(1:nbins),erra(1:nbins),
     .                  sigb(1:nbins),errb(1:nbins),
     .                  sigb1(1:nbins),errb1(1:nbins),
     .                  sigb2(1:nbins),errb2(1:nbins),
     .                  sigb3(1:nbins),errb3(1:nbins),
     .                  sigb0(1:nbins),errb0(1:nbins),
     .                  sigc1(1:nbins),errc1(1:nbins),
     .                  sigc2(1:nbins),errc2(1:nbins),
     .                  sigc3(1:nbins),errc3(1:nbins),
     .                  sigc4(1:nbins),errc4(1:nbins),
     .                  sigc5(1:nbins),errc5(1:nbins),
     .                  sigc6(1:nbins),errc6(1:nbins),
     .                  sigc0(1:nbins),errc0(1:nbins),
     .                  sigc(1:nbins),errc(1:nbins)
      common/betaf/b0,b1,b2

      alstp1 = alsarr(0,1)
      alstp2 = alsarr(0,2)
      alstp3 = alsarr(0,3)

      do k=1,idata
         fext = fextdata(k)
         ibin = ibindata(k)
         call fillhisto(froot,ichar,fext,ibin)
         open(31,file=froot(11)//'.xNNLO.'//ichar//fext(1:2))
         open(32,file=froot(11)//'.mudep.'//ichar//fext(1:2))
         open(33,file=froot(11)//'.muran.'//ichar//fext(1:2))
         open(34,file=froot(11)//'.histA.'//ichar//fext(1:2))
         open(35,file=froot(11)//'.histB.'//ichar//fext(1:2))
         open(36,file=froot(11)//'.histC.'//ichar//fext(1:2))

         do j=1,ibin
            sig0(j) = alstp1*siga(j)
            err0(j) = alstp1*erra(j)
            sig1(j) = alstp2*siga(j)+alstp2**2*sigb(j)
            err1(j) = dsqrt(alstp2**2*erra(j)**2 + alstp2**4*errb(j)**2)
            sig2(j) = alstp3*siga(j)+alstp3**2*sigb(j)+alstp3**3*sigc(j)
            err2(j) = alstp3**3*errc(j)
            
            sig0min = sig0(j)
            sig0max = sig0(j)
            sig1min = sig1(j)
            sig1max = sig1(j)
            sig2min = sig2(j)
            sig2max = sig2(j)
            
            do i = -10,10
               sig0i = alsarr(i,1)*siga(j)
               sig1i = alsarr(i,2)*siga(j)
     .              +alsarr(i,2)**2*(sigb(j)+b0*dlarr(i)*siga(j))
               sig2i = alsarr(i,3)*siga(j)
     .              +alsarr(i,3)**2*(sigb(j)+b0*dlarr(i)*siga(j))
     .              +alsarr(i,3)**3*(sigc(j)+2d0*sigb(j)*b0*dlarr(i)
     .              +siga(j)*(b0**2*dlarr(i)**2+b1*dlarr(i)))
               if (sig0i.lt.sig0min) sig0min = sig0i
               if (sig0i.gt.sig0max) sig0max = sig0i
               if (sig1i.lt.sig1min) sig1min = sig1i
               if (sig1i.gt.sig1max) sig1max = sig1i
               if (sig2i.lt.sig2min) sig2min = sig2i
               if (sig2i.gt.sig2max) sig2max = sig2i
            enddo
            therr0 = dabs(sig0max-sig0min)/2d0
            therr1 = dabs(sig1max-sig1min)/2d0
            therr2 = dabs(sig2max-sig2min)/2d0
            thmean0 = (sig0max+sig0min)/2d0
            thmean1 = (sig1max+sig1min)/2d0
            thmean2 = (sig2max+sig2min)/2d0
            span0 = 0d0
            span1 = 0d0
            span2 = 0d0
            
            if (sig0(j).gt.0d0) span0 = therr0/sig0(j)
            if (sig1(j).gt.0d0) span1 = therr1/sig1(j)
            if (sig2(j).gt.0d0) span2 = therr2/sig2(j)
            
            write(31,101) t(j),sig0(j),sig1(j),sig2(j),err2(j)
            write(32,104) t(j),thmean0,therr0,thmean1,therr1,
     .           thmean2,therr2
            write(33,103) t(j),100d0*span0,100d0*span1,100d0*span2
            write(34,105) t(j),siga(j),erra(j)
            write(35,105) t(j),sigb0(j),errb0(j)
            write(36,105) t(j),sigc0(j),errc0(j)
            
            if (sig0i.lt.sig0min) sig0min = sig0i
            if (sig0i.gt.sig0max) sig0max = sig0i
         enddo
         close(31)
         close(32)
         close(33)
         close(34)
         close(35)
         close(36)
      enddo

 101  format(3x,f11.6,1pe12.4,1pe12.4,1pe12.4,1pe12.4)
 103  format(3x,f11.6,3f12.4)
 104  format(3x,f11.6,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4)
 105  format(3x,f11.6,1pe12.4,1pe12.4)
      return
      end


      subroutine fillhisto(froot,ichar,fext,ibin)
      implicit real*8(a-h,o-z)
      character*1 ichar
      character*2 fext
      character*12 froot(11)
      parameter(nbins=400)
      logical fi(1:10)
      common/lqcd/rlambda,rnf
      common/histo/t(1:nbins),siga(1:nbins),erra(1:nbins),
     .                  sigb(1:nbins),errb(1:nbins),
     .                  sigb1(1:nbins),errb1(1:nbins),
     .                  sigb2(1:nbins),errb2(1:nbins),
     .                  sigb3(1:nbins),errb3(1:nbins),
     .                  sigb0(1:nbins),errb0(1:nbins),
     .                  sigc1(1:nbins),errc1(1:nbins),
     .                  sigc2(1:nbins),errc2(1:nbins),
     .                  sigc3(1:nbins),errc3(1:nbins),
     .                  sigc4(1:nbins),errc4(1:nbins),
     .                  sigc5(1:nbins),errc5(1:nbins),
     .                  sigc6(1:nbins),errc6(1:nbins),
     .                  sigc0(1:nbins),errc0(1:nbins),
     .                  sigc(1:nbins),errc(1:nbins)

      do i=1,10
         inquire(file=froot(i)//'.'//ichar//fext(1:2),EXIST=fi(i))
         if (.not.fi(i)) then 
            write(6,*) 'Missing file: ', froot(i)//'.'//ichar//fext(1:2)
            stop
         endif
      enddo

      open(11,file=froot(1)//'.'//ichar//fext(1:2))
      open(12,file=froot(2)//'.'//ichar//fext(1:2))
      open(13,file=froot(3)//'.'//ichar//fext(1:2))
      open(14,file=froot(4)//'.'//ichar//fext(1:2))
      open(21,file=froot(5)//'.'//ichar//fext(1:2))
      open(22,file=froot(6)//'.'//ichar//fext(1:2))
      open(23,file=froot(7)//'.'//ichar//fext(1:2))
      open(24,file=froot(8)//'.'//ichar//fext(1:2))
      open(25,file=froot(9)//'.'//ichar//fext(1:2))
      open(26,file=froot(10)//'.'//ichar//fext(1:2))
      
      do j=1,ibin
         read(11,*) t(j),siga(j),erra(j)
         read(12,*) t(j),sigb1(j),errb1(j)
         read(13,*) t(j),sigb2(j),errb2(j)
         read(14,*) t(j),sigb3(j),errb3(j)
         read(21,*) t(j),sigc1(j),errc1(j)
         read(22,*) t(j),sigc2(j),errc2(j)
         read(23,*) t(j),sigc3(j),errc3(j)
         read(24,*) t(j),sigc4(j),errc4(j)
         read(25,*) t(j),sigc5(j),errc5(j)
         read(26,*) t(j),sigc6(j),errc6(j)
         sigb0(j) = sigb1(j)+sigb2(j)+sigb3(j)
         errb0(j) = dsqrt(errb1(j)**2+errb2(j)**2+errb3(j)**2)
         sigc0(j) = sigc1(j)+sigc2(j)+sigc3(j)+sigc4(j)
     .             +sigc5(j)+sigc6(j)
         errc0(j) = dsqrt(errc1(j)**2+errc2(j)**2+errc3(j)**2
     .                   +errc4(j)**2+errc5(j)**2+errc6(j)**2)
         sigb(j) = sigb0(j) - 2d0*siga(j) 
         sigc(j) = sigc0(j) - 2d0*sigb0(j)
     .             + (-3.94282959d0+0.46118159d0*rnf)*siga(j)
         errc(j) = errc0(j)
         errb(j) = dsqrt(errb0(j)**2 + 4d0*erra(j))
      enddo

      close(11)
      close(12)
      close(13)
      close(14)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      return
      end




      subroutine alini(afix,rmu,iord)
      implicit real*8(a-h,o-z)
      common/lqcd/rlambda,rnf
      parameter(pi=3.141592653589793238d0)
c---- poor version of fixpoint algorithm
      
      ca=3d0
      cf=4d0/3d0
      tr=rnf/2d0

      b0 = 1d0/6d0*(11d0*ca-4d0*tr)
      dl0 = 2d0*pi/b0/afix
      rlambda = dexp(-dl0/2)*rmu
      al0 = alphas(rmu,iord)
      do i=1,50
      rlambda = rlambda*(afix/al0)**4
      al0 = alphas(rmu,iord)
c      write(6,*) i,rlambda,al0
      enddo
c      write(6,*) rlambda,al0

      return 
      end

      real*8 function alphas(q,iord)
      implicit real*8(a-h,o-z)
      common/lqcd/rlambda,rnf
      common/betaf/b0,b1,b2
      parameter(pi=3.141592653589793238d0)

      dl = 2d0*dlog(q/rlambda)
      dll = dlog(dl)

      ca=3d0
      cf=4d0/3d0
      tr=rnf/2d0

      b0 = 1d0/6d0*(11d0*ca-4d0*tr)
      b1 = 1d0/6d0*(17d0*ca*ca -10d0*ca*tr -6d0*cf*tr)
      b2 = 1d0/432d0*(2857d0*ca*ca*ca +108d0*cf*cf*tr -1230d0*cf*ca*tr
     .    -2830d0*ca*ca*tr +264d0*cf*tr*tr +316d0*ca*tr*tr)

      i1=0
      i2=0
      if (iord.ge.2) i1=1
      if (iord.ge.3) i2=1
      

      alphas = 2d0*pi/b0/dl*(1d0- i1* b1/b0/b0*dll/dl
     .   + i2* 1d0/b0/b0/dl/dl*((b1/b0)**2*(dll*dll-dll-1d0)+b2/b0)) 

      return 
      end

