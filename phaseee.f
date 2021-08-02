c---- phase space generators for 
c---- phase3ee:            gamma*(q) -> p1,p2,p3        
c---- phase4ee:            gamma*(q) -> p1,p2,p3,p4
c---- phase5aee,phase5bee: gamma*(q) -> p1,p2,p3,p4,p5

*
************************************************************************
*
      subroutine pick(itype,s,smin,smax,r1,wt)
      implicit real*8(a-h,o-z)
*
* itype=1  pick linearly
* itype=2  pick logarithmically
*
      if(itype.eq.1)then
        s =smin+r1*(smax-smin)
        wt=wt*(smax-smin)
      endif
      if(itype.eq.2)then
        s =smin*(smax/smin)**r1
        wt=wt*log(smax/smin)*s
      endif
      return
      end


      subroutine pickcos(costh,y0,r1,wt)
      implicit real*8(a-h,o-z)
      amax = 0.5d0*dlog((2d0+y0)/y0)
      amin = -amax
      a = amin+r1*(amax-amin)
      costh = tanh(a)
      wt = wt*(1d0-costh)*(1d0+costh)*(amax-amin)
      return
      end

      subroutine pickres(s,smin,smax,rm,rg,r1,wt)
      implicit real*8(a-h,o-z)
      zmin = atan((smin-rm**2)/rg/rm)/rg/rm
      zmax = atan((smax-rm**2)/rg/rm)/rg/rm
      z = zmin + r1*(zmax-zmin)
      wt = wt*(zmax-zmin)
      s = rm*rg*tan(rm*rg*z)+rm**2
      dv = 1d0/((s-rm**2)**2+rm**2*rg**2)
      wt = wt/dv
      return
      end



      subroutine phase3ee(x,wt,ifail)
      implicit real*8(a-h,o-z)
      dimension p1(4),p2(4),p3(4),v(4),v1(4),q(4)
      dimension p1p(4),p2p(4),p3p(4),pap(4),pbp(4),qp(4)
      dimension p1a(4),p2a(4),p3a(4)
      dimension p1ap(4),p2ap(4),p3ap(4)
      dimension p12(4),p12p(4)
      dimension pqa(4)
      dimension x(10)
      dimension bmat(4,4),ubmat(4,4),rmat(4,4),urmat(4,4),fmat(4,4)
      parameter(pi=3.141592653589793238d0)
      common/masses/rm2(1:5),shat
      common/pin/pa(4),pb(4)
      common/pmom/p(4,5)
      common/invar3/sij3(3,3),tij3(2,3)
      common/pmomang3/pang3(4,3,1)
      common/invarang3/sijang3(3,3,1)
      common /tcuts/ymin,y0

*
* phase space for 1->3 decay, q -> 1+2+3
*
* rest frame of q
* 
* massless case: integrates to P_3 = 1/256/pi**3*sab
*
      r1 = x(1)
      r2 = x(2)
      ifail =1
      wt = 0d0
      wtps = 1d0
      call p3generic(r1,r2,shat,rm2(1),rm2(2),rm2(3),p1,p2,p3,wtps)
      if (wtps.eq.0d0) return
      do i=1,4
         p(i,1) = p1(i)
         p(i,2) = p2(i)
         p(i,3) = p3(i)
         pang3(i,1,1) = p1(i)
         pang3(i,2,1) = p2(i)
         pang3(i,3,1) = p3(i)
      enddo
      call fillinv(3,pang3(1:4,1:3,1),sijang3(1:3,1:3,1))
      ifail =0
      wt = wtps

      return
      end


      subroutine fillcommon3pee
      implicit real*8(a-h,o-z)
      common/pmomang3/pang3(4,3,1)
      common/pmomang3g/pang3g(4,3,1)
      common/invarang3/sijang3(3,3,1)
      common/eventmom3p/pevt3(4,3,1)
      common/eventinv3p/sevt3(3,3,1),ievt3(1)

c---- pang3,pang3g and sija,tija are input 
c---- this subroutine fills 
c----       pevt3(component,pnum,ievent)             : 4-momenta of events
c----       sevt3(i,j,ievent),tevt3([ab],i,ievent)   : invariants of events
c----
c---- ievent = 1 fixed

      do j=1,3
         do k=1,4
            pevt3(k,j,1)    = pang3(k,j,1)
         enddo
      enddo
      do j=1,3
         do k=1,3
            sevt3(j,k,1)    = sijang3(j,k,1)
         enddo
      enddo
      ievt3(1) = 1
      return
      end


      
      subroutine phase4ee(x,wt,ifail)
      implicit real*8(a-h,o-z)
      dimension p1(4),p2(4),p3(4),p4(4),v(4),v1(4),q(4)
      dimension p1p(4),p2p(4),p3p(4),p4p(4)
      dimension p1a(4),p2a(4),p3a(4),p4a(4)
      dimension p1ap(4),p2ap(4),p3ap(4),p4ap(4)
      dimension p12(4),p12p(4)
      dimension p34(4),p34p(4),p34g(4),papg(4)
      dimension x(10)
      dimension bmat(4,4),ubmat(4,4),rmat(4,4),urmat(4,4),fmat(4,4)
      parameter(pi=3.141592653589793238d0)
      common/masses/rm2(1:5),shat
      common/pin/pa(4),pb(4)
      common/pmom/p(4,5)
      common/invar4/sij4(4,4)
      common/pmomang4/pang4(4,4,2)
      common/invarang4/sijang4(4,4,2),iacc(2)
      common /tcuts/ymin,y0
c---  corresponds to the phase space(3,4;12)
*
* phase space for 1->4 decay, q -> 1+2+3+4 with (34)->3+4
*
* (34) system is rotated around p1+p2 axis 
*
* massless case: integrates to P_4 = 1/24576/pi**5
*
      r1 = x(1)
      r2 = x(2)
      r3 = x(3)
      r4 = x(4)
      r5 = x(5)
      ifail =1
      wt = 0d0
      wtps = 1d0

      s34min = dmax1(y0*shat,(dsqrt(rm2(3))+dsqrt(rm2(4)))**2)
      s34max = dmin1((dsqrt(shat)-dsqrt(rm2(1))-dsqrt(rm2(2)))**2
     .              ,(1d0-5d0*y0)*shat)
      if (s34min.ge.s34max) return
      call pick(2,s34,s34min,s34max,r1,wtps)
      wtps = wtps/2d0/pi
      call p3generic(r2,r3,shat,rm2(1),rm2(2),s34,p1,p2,p34,wtps)
      if (wtps.eq.0d0) return
      do i=1,4
         p12(i) = p1(i)+p2(i)
      enddo

      costhmin = -1d0
      costhmax = 1d0
      phimin = 0d0
      phimax = 2d0*pi
      call pick(1,costh,costhmin,costhmax,r4,wtps)
      call pick(1,phi,phimin,phimax,r5,wtps)
      wtps = wtps*dsqrt(dlambda(s34,rm2(3),rm2(4)))/s34/32d0/pi**2      
      call makeps2com(s34,costh,phi,rm2(3),rm2(4),p3p,p4p,p3ap,p4ap)      

      call boostrest(p34,bmat)
      call unboostrest(p34,ubmat)
      do i=1,4
         v(i) = 0d0
         p12p(i) = 0d0
         do j=1,4
            p12p(i) = p12p(i) + bmat(i,j)*p12(j)
         enddo
      enddo

      call unrotatetoz(p12p,urmat)
      
      do i=1,4
         do j=1,4
            fmat(i,j) = 0d0
            do k=1,4
               fmat(i,j) = fmat(i,j)+ubmat(i,k)*urmat(k,j)
            enddo
         enddo
      enddo

      do i=1,4
         p3(i) = 0d0
         p4(i) = 0d0
         p3a(i) = 0d0
         p4a(i) = 0d0
         do j=1,4
            p3(i) = p3(i) + fmat(i,j)*p3p(j)
            p4(i) = p4(i) + fmat(i,j)*p4p(j)
            p3a(i) = p3a(i) + fmat(i,j)*p3ap(j)
            p4a(i) = p4a(i) + fmat(i,j)*p4ap(j)
         enddo
      enddo
     
      do i=1,4
         p(i,1) = p1(i)
         p(i,2) = p2(i)
         p(i,3) = p3(i)
         p(i,4) = p4(i)
         p1a(i) = p1(i)
         p2a(i) = p2(i)
         pang4(i,1,1) = p1(i)
         pang4(i,2,1) = p2(i)
         pang4(i,3,1) = p3(i)
         pang4(i,4,1) = p4(i)
         pang4(i,1,2) = p1a(i)
         pang4(i,2,2) = p2a(i)
         pang4(i,3,2) = p3a(i)
         pang4(i,4,2) = p4a(i)
      enddo
      call fillinv(4,p,sij4)
      do i=1,2
         call fillinv(4,pang4(1:4,1:4,i),sijang4(1:4,1:4,i))
      enddo
      do i=1,2
         call acceptcuts4(sijang4(1:4,1:4,i),iacc(i))
      enddo
      if (iacc(1)+iacc(2).eq.0) return
      ifail =0
      wt = wtps

      return
      end


      subroutine acceptcuts4(sij4,ipass)
      implicit real*8(a-h,o-z)
      dimension sij4(1:4,1:4)
      common/masses/rm2(1:5),shat
      common /tcuts/ymin,y0
      
      ipass = 0

      s12 = sij4(1,2)
      s13 = sij4(1,3)
      s14 = sij4(1,4)
      s23 = sij4(2,3)
      s24 = sij4(2,4)
      s34 = sij4(3,4)
      rinvmin = min(s12,s13,s14,s23,s24,s34)
      if (rinvmin.lt.y0*shat) return      
      ipass = 0
* single out on 3||4
      if (rinvmin.eq.s34) ipass=1
      return
      end


      subroutine fillcommon4pee
      implicit real*8(a-h,o-z)
      dimension i1(1:4,1:6)
      common/pmomang4/pang4(4,4,2)
      common/invarang4/sijang4(4,4,2),iacc(2)
      common/eventmom4p/pevt4(4,4,12)
      common/eventinv4p/sevt4(4,4,12),ievt4(12)

      data i1/1,2,3,4,
     .        1,3,2,4,
     .        1,4,2,3,
     .        2,3,1,4,
     .        2,4,1,3,
     .        3,4,1,2/

c---- pang4 and sija are input and correspond to an event in (3,4;12)
c---- this subroutine fills 
c----       pevt4(component,pnum,ievent)     : 4-momenta of events
c----       sevt4(i,j,ievent)                : invariants of events
c----
c---- ievent
c----  1: (3,4;12)
c----  2: (2,4;13)
c----  3: (2,3;14)
c----  4: (1,4;23)
c----  5: (1,3;24)
c----  6: (1,2;34)
c----
c----  7-12 (i,j;x) with (i,j) system rotated around x-axis 
c----
c---- full phave space coverage is obtained by summing all 6 events
c---- angular average is obtained by averaging:
c----        M2 = 1/2* ( sum_ievent=1..12   M2(ievent))
c----
      do i=1,6
         do j=1,4
            do k=1,4
               pevt4(k,i1(j,i),i)    = pang4(k,j,1)
               pevt4(k,i1(j,i),i+6)  = pang4(k,j,2)
            enddo
         enddo
      enddo

      do i=1,6
         do j=1,4
            do k=1,4
               sevt4(i1(j,i),i1(k,i),i)    = sijang4(j,k,1)
               sevt4(i1(j,i),i1(k,i),i+6)  = sijang4(j,k,2)
            enddo
         enddo
         ievt4(i)   = iacc(1)
         ievt4(i+6) = iacc(2)
      enddo

      return
      end



      subroutine phase5aee(x,wt,ifail)
      implicit real*8(a-h,o-z)
      dimension p1(4),p2(4),p3(4),p4(4),p5(4),v(4),q(4),v1(4)
      dimension p1p(4),p2p(4),p3p(4),p4p(4),p5p(4),pap(4),
     .          p3pp(4),p3app(4)
      dimension p1ap(4),p2ap(4),p3ap(4),p4ap(4),p5ap(4)
      dimension p3a(4),p4a(4),p5a(4)
      dimension p345(4),p45p(4),p45(4)
      dimension p345p(4)
      dimension p45a(4),p45ap(4),p4t(4),p5t(4),p4ta(4),p5ta(4)
      dimension p12(4),p12p(4)
      dimension x(10)
      dimension bmat(4,4),ubmat(4,4),rmat(4,4),urmat(4,4),fmat(4,4)
      parameter(pi=3.141592653589793238d0)
      common/masses/rm2(1:5),shat
      common/pin/pa(4),pb(4)
      common/pmom/p(4,5)
      common/pmomang5/pang5(4,5,4)
      common/invar5/sij5(5,5)
      common/invarang5/sijang5(5,5,4),iacc(4)
      common/tcuts/ymin,y0

* this is (3;4,5;12)
*
* phase space for 1->5 scattering, q -> 1+2+3+4+5 
* with sequential splitting (345) -> 3+(45) -> 3+4+5 
* 
* (345) system is rotated with respect to p12 axis
* (45)  system is rotated with respect to p3 axis
*
* default momenta             : pang(i,j,1)
* (45) rotated                : pang(i,j,2)
* (345) rotated               : pang(i,j,3)
* (345) rotated, (45) rotated : pang(i,j,4)
* 
* massless case: integrates to P_5 = 1/4718592/pi**7*sab^3 = 0.7016789753e-10*sab^3
*
      r1 = x(1)
      r2 = x(2)
      r3 = x(3)
      r4 = x(4)
      r5 = x(5)
      r6 = x(6)
      r7 = x(7)
      r8 = x(8)
      wt = 0d0
      ifail =1
      wtps = 1d0

      s345min = dmax1(y0*shat,
     .           (dsqrt(rm2(3))+dsqrt(rm2(4))+dsqrt(rm2(5)))**2)
      s345max = dmin1((dsqrt(shat)-dsqrt(s12))**2,(1d0-y0)*shat)
      if (s345max.le.s345min) return
      call pick(2,s345,s345min,s345max,r1,wtps)
      wtps = wtps/2d0/pi
      s45min = dmax1(y0*shat,(dsqrt(rm2(4))+dsqrt(rm2(5)))**2)
      s45max = dmin1(s345,(1d0-y0)*shat)
      if (s45max.le.s45min) return
      call pick(2,s45,s45min,s45max,r2,wtps)
      wtps = wtps/2d0/pi

      call p3generic(r3,r4,shat,rm2(1),rm2(2),s345,p1,p2,p345,wtps)
      if (wtps.eq.0d0) return
      do i=1,4
         p12(i) = p1(i)+p2(i)
      enddo

      costhmin = -1d0
      costhmax = 1d0
      phimin = 0d0
      phimax = 2d0*pi
      call pick(1,costh,costhmin,costhmax,r5,wtps)
      call pick(1,phi,phimin,phimax,r6,wtps)
      call makeps2com(s345,costh,phi,rm2(3),s45,p3p,p45p,p3ap,p45ap)      
      wtps = wtps*dsqrt(dlambda(s345,rm2(3),s45))/s345/32d0/pi**2

      call pick(1,costh,costhmin,costhmax,r7,wtps)
      call pick(1,phi,phimin,phimax,r8,wtps)
      call makeps2com(s45,costh,phi,rm2(4),rm2(5),p4p,p5p,p4ap,p5ap)      
      wtps = wtps*dsqrt(dlambda(s45,rm2(4),rm2(5)))/s45/32d0/pi**2

c---- rotate the (345) system
      call boostrest(p345,bmat)
      call unboostrest(p345,ubmat)
      do i=1,4
         v(i) = 0d0
         p12p(i) = 0d0
         do j=1,4
            p12p(i) = p12p(i) + bmat(i,j)*p12(j)
         enddo
      enddo
      call unrotatetoz(p12p,urmat)
      do i=1,4
         do j=1,4
            fmat(i,j) = 0d0
            do k=1,4
               fmat(i,j) = fmat(i,j)+ubmat(i,k)*urmat(k,j)
            enddo
         enddo
      enddo
      do i=1,4
         p3(i) = 0d0
         p45(i) = 0d0
         p3a(i) = 0d0
         p45a(i) = 0d0
         v(i) = 0d0
         do j=1,4
            p3(i) = p3(i) + fmat(i,j)*p3p(j)
            p45(i) = p45(i) + fmat(i,j)*p45p(j)
            p3a(i) = p3a(i) + fmat(i,j)*p3ap(j)
            p45a(i) = p45a(i) + fmat(i,j)*p45ap(j)
            v(i) = v(i) + fmat(i,j)*pap(j)
         enddo
      enddo

c---- rotate the (45) system
      call boostrest(p45,bmat)
      call unboostrest(p45,ubmat)
      do i=1,4
         v(i) = 0d0
         p3pp(i) = 0d0
         do j=1,4
            p3pp(i) = p3pp(i) + bmat(i,j)*p3(j)
         enddo
      enddo
      call unrotatetoz(p3pp,urmat)
      do i=1,4
         do j=1,4
            fmat(i,j) = 0d0
            do k=1,4
               fmat(i,j) = fmat(i,j)+ubmat(i,k)*urmat(k,j)
            enddo
         enddo
      enddo
      do i=1,4
         p4(i) = 0d0
         p5(i) = 0d0
         p4a(i) = 0d0
         p5a(i) = 0d0
         v(i) = 0d0
         do j=1,4
            p4(i) = p4(i) + fmat(i,j)*p4p(j)
            p5(i) = p5(i) + fmat(i,j)*p5p(j)
            p4a(i) = p4a(i) + fmat(i,j)*p4ap(j)
            p5a(i) = p5a(i) + fmat(i,j)*p5ap(j)
            v(i) = v(i) + fmat(i,j)*pap(j)
         enddo
      enddo

c---- rotate the (45ang) system
      call boostrest(p45a,bmat)
      call unboostrest(p45a,ubmat)
      do i=1,4
         v(i) = 0d0
         p3app(i) = 0d0
         do j=1,4
            p3app(i) = p3app(i) + bmat(i,j)*p3a(j)
         enddo
      enddo
      call unrotatetoz(p3app,urmat)
      do i=1,4
         do j=1,4
            fmat(i,j) = 0d0
            do k=1,4
               fmat(i,j) = fmat(i,j)+ubmat(i,k)*urmat(k,j)
            enddo
         enddo
      enddo
      do i=1,4
         p4t(i) = 0d0
         p5t(i) = 0d0
         p4ta(i) = 0d0
         p5ta(i) = 0d0
         v(i) = 0d0
         do j=1,4
            p4t(i) = p4t(i) + fmat(i,j)*p4p(j)
            p5t(i) = p5t(i) + fmat(i,j)*p5p(j)
            p4ta(i) = p4ta(i) + fmat(i,j)*p4ap(j)
            p5ta(i) = p5ta(i) + fmat(i,j)*p5ap(j)
            v(i) = v(i) + fmat(i,j)*pap(j)
         enddo
      enddo

      do i=1,4
         p(i,1) = p1(i)
         p(i,2) = p2(i)
         p(i,3) = p3(i)
         p(i,4) = p4(i)
         p(i,5) = p5(i)
         pang5(i,1,1) = p1(i)
         pang5(i,2,1) = p2(i)
         pang5(i,3,1) = p3(i)
         pang5(i,4,1) = p4(i)
         pang5(i,5,1) = p5(i)
         pang5(i,1,2) = p1(i)
         pang5(i,2,2) = p2(i)
         pang5(i,3,2) = p3(i)
         pang5(i,4,2) = p4a(i)
         pang5(i,5,2) = p5a(i)
         pang5(i,1,3) = p1(i)
         pang5(i,2,3) = p2(i)
         pang5(i,3,3) = p3a(i)
         pang5(i,4,3) = p4t(i)
         pang5(i,5,3) = p5t(i)
         pang5(i,1,4) = p1(i)
         pang5(i,2,4) = p2(i)
         pang5(i,3,4) = p3a(i)
         pang5(i,4,4) = p4ta(i)
         pang5(i,5,4) = p5ta(i)
      enddo

      do i=1,4
         call fillinv(5,pang5(1:4,1:5,i),sijang5(1:5,1:5,i))
      enddo

      do i=1,4
        call acceptcuts5a(sijang5(1:5,1:5,i),iacc(i))
      enddo
      if (iacc(1)+iacc(2)+iacc(3)+iacc(4).eq.0) return
      ifail =0
      wt = wtps

      return
      end
      

      subroutine acceptcuts5a(sij5,ipass)
      implicit real*8(a-h,o-z)
      dimension sij5(1:5,1:5)
      common/masses/rm2(1:5),shat
      common /tcuts/ymin,y0
      
      ipass = 0
      
      s12 = sij5(1,2)
      s13 = sij5(1,3)
      s14 = sij5(1,4)
      s15 = sij5(1,5)
      s23 = sij5(2,3)
      s24 = sij5(2,4)
      s25 = sij5(2,5)
      s34 = sij5(3,4)
      s35 = sij5(3,5)
      s45 = sij5(4,5)

      smin1 = shat
      smin2 = shat

      do i=1,5
         do j=i+1,5
            if (sij5(i,j).lt.smin2) then
               smin2 = sij5(i,j)
               if (smin2.lt.smin1) then
                  stmp  = smin1
                  smin1 = smin2
                  smin2 = stmp
               endif
            endif
         enddo
      enddo

      if (smin1.lt.y0*shat) return
      if (smin1.eq.s45.and.smin2.eq.s34) ipass=1
      if (smin1.eq.s45.and.smin2.eq.s35) ipass=1

      return
      end


      subroutine phase5bee(x,wt,ifail)
      implicit real*8(a-h,o-z)
      dimension p1(4),p2(4),p3(4),p4(4),p5(4),v(4),q(4),v1(4)
      dimension p1p(4),p2p(4),p3p(4),p4p(4),p5p(4)
      dimension p1a(4),p2a(4),p3a(4),p4a(4),p5a(4)
      dimension p1ap(4),p2ap(4),p3ap(4),p4ap(4),p5ap(4)
      dimension p12(4),p34(4)
      dimension x(10)
      dimension bmat(4,4),ubmat(4,4),rmat(4,4),urmat(4,4),fmat(4,4)
      parameter(pi=3.141592653589793238d0)
      common/masses/rm2(1:5),shat
      common/pin/pa(4),pb(4)
      common/pmom/p(4,5)
      common/pmomang5/pang5(4,5,4)
      common/invar5/sij5(5,5)
      common/invarang5/sijang5(5,5,4),iacc(4)
      common/tcuts/ymin,y0


* this is (1,2;5;3,4;5)
*
* phase space for 1->5 decay, q -> 1+2+3+4+5 with q=p1+p2
* with splittings (12) -> 1+2; (34)->3+4
* 
* (12) system is rotated with respect to p5 axis
* (34) system is rotated with respect to p5 axis
*
* default momenta            : pang(i,j,1)
* (12) rotated               : pang(i,j,2)
* (34) rotated               : pang(i,j,3)
* (12) rotated, (34) rotated : pang(i,j,4)
* 
* 
* massless case: integrates to P_5 = 1/4718592/pi**7*sab^3 = 0.7016789753e-10*sab^3
*
      r1 = x(1)
      r2 = x(2)
      r3 = x(3)
      r4 = x(4)
      r5 = x(5)
      r6 = x(6)
      r7 = x(7)
      r8 = x(8)
      wt = 0d0
      ifail =1
      wtps = 1d0


      s12min = dmax1(y0*shat,(dsqrt(rm2(1))+dsqrt(rm2(2)))**2)
      s12max = dmin1((dsqrt(shat)-dsqrt(rm2(3))-dsqrt(rm2(4))
     .              -dsqrt(rm2(5)))**2,(1d0-5d0*y0)*shat)
      if (s12min.ge.s12max) return
      call pick(2,s12,s12min,s12max,r1,wtps)
      wtps = wtps/2d0/pi

      s34min = dmax1(y0*shat,(dsqrt(rm2(3))+dsqrt(rm2(4)))**2)
      s34max = dmin1((dsqrt(shat)-dsqrt(s12)-dsqrt(rm2(5)))**2
     .        ,(1d0-y0)*shat)
      if (s34max.le.s34min) return
      call pick(2,s34,s34min,s34max,r2,wtps)
      wtps = wtps/2d0/pi

      call p3generic(r3,r4,shat,s12,s34,rm2(5),p12,p34,p5,wtps)
      if (wtps.eq.0d0) return

      costhmin = -1d0
      costhmax = 1d0
      phimin = 0d0
      phimax = 2d0*pi
      call pick(1,costh,costhmin,costhmax,r5,wtps)
      call pick(1,phi,phimin,phimax,r6,wtps)
      call makeps2com(s12,costh,phi,rm2(1),rm2(2),p1p,p2p,p1ap,p2ap)      
      wtps = wtps*dsqrt(dlambda(s12,rm2(1),rm2(2)))/s12/32d0/pi**2

      costhmin = -1d0
      costhmax = 1d0
      phimin = 0d0
      phimax = 2d0*pi
      call pick(1,costh,costhmin,costhmax,r7,wtps)
      call pick(1,phi,phimin,phimax,r8,wtps)
      call makeps2com(s34,costh,phi,rm2(3),rm2(4),p3p,p4p,p3ap,p4ap)      
      wtps = wtps*dsqrt(dlambda(s34,rm2(3),rm2(4)))/s34/32d0/pi**2

c---- rotate the (12) system
      call boostrest(p12,bmat)
      call unboostrest(p12,ubmat)
      do i=1,4
         v(i) = 0d0
         p5p(i) = 0d0
         do j=1,4
            p5p(i) = p5p(i) + bmat(i,j)*p5(j)
         enddo
      enddo
      call unrotatetoz(p5p,urmat)
      do i=1,4
         do j=1,4
            fmat(i,j) = 0d0
            do k=1,4
              fmat(i,j) = fmat(i,j)+ubmat(i,k)*urmat(k,j)
            enddo
         enddo
      enddo
      do i=1,4
         p1(i) = 0d0
         p2(i) = 0d0
         p1a(i) = 0d0
         p2a(i) = 0d0
         v(i) = 0d0
         do j=1,4
            p1(i) = p1(i) + fmat(i,j)*p1p(j)
            p2(i) = p2(i) + fmat(i,j)*p2p(j)
            p1a(i) = p1a(i) + fmat(i,j)*p1ap(j)
            p2a(i) = p2a(i) + fmat(i,j)*p2ap(j)
            v(i) = v(i) + fmat(i,j)*p5p(j)
         enddo
      enddo

c---- rotate the (34) system
      call boostrest(p34,bmat)
      call unboostrest(p34,ubmat)
      do i=1,4
         v(i) = 0d0
         p5p(i) = 0d0
         do j=1,4
            p5p(i) = p5p(i) + bmat(i,j)*p5(j)
         enddo
      enddo
      call unrotatetoz(p5p,urmat)
      do i=1,4
         do j=1,4
            fmat(i,j) = 0d0
            do k=1,4
               fmat(i,j) = fmat(i,j)+ubmat(i,k)*urmat(k,j)
            enddo
         enddo
      enddo
      do i=1,4
         p3(i) = 0d0
         p4(i) = 0d0
         p3a(i) = 0d0
         p4a(i) = 0d0
         v(i) = 0d0
         do j=1,4
            p3(i) = p3(i) + fmat(i,j)*p3p(j)
            p4(i) = p4(i) + fmat(i,j)*p4p(j)
            p3a(i) = p3a(i) + fmat(i,j)*p3ap(j)
            p4a(i) = p4a(i) + fmat(i,j)*p4ap(j)
            v(i) = v(i) + fmat(i,j)*p5p(j)
         enddo
      enddo

      do i=1,4
         p(i,1) = p1(i)
         p(i,2) = p2(i)
         p(i,3) = p3(i)
         p(i,4) = p4(i)
         p(i,5) = p5(i)
         pang5(i,1,1) = p1(i)
         pang5(i,2,1) = p2(i)
         pang5(i,3,1) = p3(i)
         pang5(i,4,1) = p4(i)
         pang5(i,5,1) = p5(i)
         pang5(i,1,2) = p1a(i)
         pang5(i,2,2) = p2a(i)
         pang5(i,3,2) = p3(i)
         pang5(i,4,2) = p4(i)
         pang5(i,5,2) = p5(i)
         pang5(i,1,3) = p1(i)
         pang5(i,2,3) = p2(i)
         pang5(i,3,3) = p3a(i)
         pang5(i,4,3) = p4a(i)
         pang5(i,5,3) = p5(i)
         pang5(i,1,4) = p1a(i)
         pang5(i,2,4) = p2a(i)
         pang5(i,3,4) = p3a(i)
         pang5(i,4,4) = p4a(i)
         pang5(i,5,4) = p5(i)
      enddo
      do i=1,4
         call fillinv(5,pang5(1:4,1:5,i),sijang5(1:5,1:5,i))
      enddo
      do i=1,4
        call acceptcuts5b(sijang5(1:5,1:5,i),iacc(i))
      enddo
      if (iacc(1)+iacc(2)+iacc(3)+iacc(4).eq.0) return
      ifail =0
      wt = wtps
      return
      end



      subroutine acceptcuts5b(sij5,ipass)
      implicit real*8(a-h,o-z)
      dimension sij5(1:5,1:5)
      common/masses/rm2(1:5),shat
      common /tcuts/ymin,y0
      
      ipass = 0
      s12 = sij5(1,2)
      s13 = sij5(1,3)
      s14 = sij5(1,4)
      s15 = sij5(1,5)
      s23 = sij5(2,3)
      s24 = sij5(2,4)
      s25 = sij5(2,5)
      s34 = sij5(3,4)
      s35 = sij5(3,5)
      s45 = sij5(4,5)

      smin1 = shat
      smin2 = shat

      do i=1,5
         do j=i+1,5
            if (sij5(i,j).lt.smin2) then
               smin2 = sij5(i,j)
               if (smin2.lt.smin1) then
                  stmp  = smin1
                  smin1 = smin2
                  smin2 = stmp
               endif
            endif
         enddo
      enddo

      if (smin1.lt.y0*shat) return
      if (smin1.eq.s12.and.smin2.eq.s34) ipass=1
      if (smin1.eq.s34.and.smin2.eq.s12) ipass=1

      return
      end

      

      subroutine fillcommon5apee
      implicit real*8(a-h,o-z)
      common/pmomang5/pang5(4,5,4)
      common/invarang5/sijang5(5,5,4),iacc(4)
      common/eventmom5ap/pevt5a(4,5,120)
      common/eventinv5ap/sevt5a(5,5,120),ievt5a(120)
      dimension i1(1:5,1:30)

      data i1/1,2,3,4,5,
     .        1,2,4,5,3,
     .        1,2,5,3,4,
     .        1,3,2,4,5,
     .        1,3,4,5,2,
     .        1,3,5,2,4,
     .        1,4,2,3,5,
     .        1,4,3,5,2,
     .        1,4,5,2,3,
     .        1,5,2,3,4,
     .        1,5,3,4,2,
     .        1,5,4,2,3,
     .        2,3,1,4,5,
     .        2,3,4,5,1,
     .        2,3,5,1,4,
     .        2,4,1,3,5,
     .        2,4,3,5,1,
     .        2,4,5,1,3,
     .        2,5,1,3,4,
     .        2,5,3,4,1,
     .        2,5,4,1,3,
     .        3,4,1,2,5,
     .        3,4,2,5,1,
     .        3,4,5,1,2,
     .        3,5,1,2,4,
     .        3,5,2,4,1,
     .        3,5,4,1,2,
     .        4,5,1,2,3,
     .        4,5,2,3,1,
     .        4,5,3,1,2/

c---- pang5 and sijang are input and correspond to an event in (3;45;12)
c---- this subroutine fills 
c----       pevt5a(ievent,component,pnum)            : 4-momenta of events
c----       sevt5a(i,j,ievent),tevt5a([ab],i,ievent) : invariants of events
c----
c----  1: (3;45;12)
c----  2: (4;35;12)
c----  3: (5;34;12)
c----  4: (2;45;13)
c----  5: (4;25;13)
c----  6: (5;24;13)
c----  7: (2;35;14)
c----  8: (3;25;14)
c----  9: (5;23;14)
c---- 10: (2;34;15)
c---- 11: (3;24;15)
c---- 12: (4;23;15)
c---- 13: (1;45;23)
c---- 14: (4;15;23)
c---- 15: (5;14;23)
c---- 16: (1;35;24)
c---- 17: (3;15;24)
c---- 18: (5;13;24)
c---- 19: (1;34;25)
c---- 20: (3;14;25)
c---- 21: (4;13;25)
c---- 22: (1;25;34)
c---- 23: (2;15;34)
c---- 24: (5;12;34)
c---- 25: (1;24;35)
c---- 26: (2;14;35)
c---- 27: (4;12;35)
c---- 28: (1;23;45)
c---- 29: (2;13;45)
c---- 30: (3;12;45)

c---- 31-60:  with (45) system rotated
c---- 61-90:  with (345) system rotated
c---- 91-120: with (45) and (345) system rotated

c---- full phave space coverage is obtained by summing ievent=1..30
c---- angular average is obtained by averaging:
c----        M2 = 1/4* ( sum_ievent=1..120   M2(ievent))
c----

      do i=1,30
         do j=1,5
            do k=1,4
               pevt5a(k,i1(j,i),i)     = pang5(k,j,1)
               pevt5a(k,i1(j,i),i+30)  = pang5(k,j,2)
               pevt5a(k,i1(j,i),i+60)  = pang5(k,j,3)
               pevt5a(k,i1(j,i),i+90)  = pang5(k,j,4)
            enddo
         enddo
      enddo

      do i=1,30
         do j=1,5
            do k=1,5
               sevt5a(i1(j,i),i1(k,i),i)    = sijang5(j,k,1)
               sevt5a(i1(j,i),i1(k,i),i+30) = sijang5(j,k,2)
               sevt5a(i1(j,i),i1(k,i),i+60) = sijang5(j,k,3)
               sevt5a(i1(j,i),i1(k,i),i+90) = sijang5(j,k,4)
            enddo
         enddo
         ievt5a(i)    = iacc(1)
         ievt5a(i+30) = iacc(2)
         ievt5a(i+60) = iacc(3)
         ievt5a(i+90) = iacc(4)
      enddo
      return
      end

      subroutine fillcommon5bpee
      implicit real*8(a-h,o-z)
      common/pmomang5/pang5(4,5,4)
      common/invarang5/sijang5(5,5,4),iacc(4)
      common/eventmom5bp/pevt5b(4,5,60)
      common/eventinv5bp/sevt5b(5,5,60),ievt5b(60)
      dimension i1(1:5,1:15)

      data i1/1,2,3,4,5,
     .        1,3,2,4,5,
     .        1,4,2,3,5,
     .        1,2,3,5,4,
     .        1,3,2,5,4,
     .        1,5,2,3,4,
     .        1,2,4,5,3,
     .        1,4,2,5,3,
     .        1,5,2,4,3,
     .        1,3,4,5,2,
     .        1,4,3,5,2,
     .        1,5,3,4,2,
     .        2,3,4,5,1,
     .        2,4,3,5,1,
     .        2,5,3,4,1/

c---- pang5 and sijang are input and correspond to an event in (12;34;5)
c---- this subroutine fills 
c----       pevt5b(ievent,component,pnum)  : 4-momenta of events
c----       sevt5b(i,j,ievent)             : invariants of events
c----
c----  1: (12;34;5)
c----  2: (13;24;5)
c----  3: (14;23;5)
c----  4: (12;35;4)
c----  5: (13;25;4)
c----  6: (15;23;4)
c----  7: (12;45;3)
c----  8: (14;25;3)
c----  9: (15;24;3)
c---- 10: (13;45;2)
c---- 11: (14;35;2)
c---- 12: (15;34;2)
c---- 13: (23;45;1)
c---- 14: (24;35;1)
c---- 15: (25;34;1)

c---- 16-30: with (12) system rotated
c---- 31-45: with (34) system rotated
c---- 46-60: with (12) and (34) system rotated

c---- full phave space coverage is obtained by summing ievent=1..15
c---- angular average is obtained by averaging:
c----        M2 = 1/4* ( sum_ievent=1..60   M2(ievent))
c----

      do i=1,15
         do j=1,5
            do k=1,4
               pevt5b(k,i1(j,i),i)    = pang5(k,j,1)
               pevt5b(k,i1(j,i),i+15) = pang5(k,j,2)
               pevt5b(k,i1(j,i),i+30) = pang5(k,j,3)
               pevt5b(k,i1(j,i),i+45) = pang5(k,j,4)
            enddo
         enddo
      enddo

      do i=1,15
         do j=1,5
            do k=1,5
               sevt5b(i1(j,i),i1(k,i),i)    = sijang5(j,k,1)
               sevt5b(i1(j,i),i1(k,i),i+15) = sijang5(j,k,2)
               sevt5b(i1(j,i),i1(k,i),i+30) = sijang5(j,k,3)
               sevt5b(i1(j,i),i1(k,i),i+45) = sijang5(j,k,4)
            enddo
         enddo
         ievt5b(i)    = iacc(1)
         ievt5b(i+15) = iacc(2)
         ievt5b(i+30) = iacc(3)
         ievt5b(i+45) = iacc(4)
      enddo
      return
      end

      function dot(a,b)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4)
      dot=a(4)*b(4)-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)
      return
      end

      subroutine boostrest(a,bmat)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4),c(4),v(4),bmat(4,4)
c---- computes Lorentz boost into rest frame of a
c---- requires ma \neq 0 

      rma2 = dot(a,a)
      if (rma2.le.0d0) then
         write(*,*) 'undue boost:',rma2 
         stop
      endif
      rma = dsqrt(rma2)
      gam = a(4)/rma     
      do i=1,4
         v(i) = a(i)/a(4)
      enddo
      v2 = v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
      bmat(4,4) = gam
      do i=1,3
         bmat(i,4) = -gam*v(i)
         bmat(4,i) = -gam*v(i)
         do j=1,3
            bmat(i,j) = (gam-1d0)*v(i)*v(j)/v2
         enddo
         bmat(i,i) = bmat(i,i)+1d0
      enddo
      return
      end


      subroutine unboostrest(a,bmat)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4),c(4),v(4),bmat(4,4)
c---- computes Lorentz boost from rest frame of a  
c---- of a into moving frame defined by the actual a
c---- requires ma \neq 0 

      rma2 = dot(a,a)
      if (rma2.le.0d0) then
         write(*,*) 'undue boost' 
         stop
      endif
      rma = dsqrt(rma2)
      gam = a(4)/rma     
      do i=1,4
         v(i) = -a(i)/a(4)
      enddo
      v2 = v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
      bmat(4,4) = gam
      do i=1,3
         bmat(i,4) = -gam*v(i)
         bmat(4,i) = -gam*v(i)
         do j=1,3
            bmat(i,j) = (gam-1d0)*v(i)*v(j)/v2
         enddo
         bmat(i,i) = bmat(i,i)+1d0
      enddo
      return
      end

      subroutine rotatetoz(a,rmat)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4),rmat(4,4),rmatz(4,4),rmaty(4,4)
      parameter(pi=3.141592653589793238d0)
c---- computes rotation into frame with a along z axis

      if(a(1).eq.0d0.and.a(2).eq.0d0.and.a(3).ge.0d0) then
         do i=1,4
            do j=1,4
               rmat(i,j) = 0d0
            enddo
            rmat(i,i) = 1d0
         enddo
         return
      endif

      rmatz(4,4) = 1d0
      do i=1,3
         rmatz(4,i) = 0d0
         rmatz(i,4) = 0d0
      enddo

      if(a(1).eq.0d0) then
         rmatz(1,1) = 0d0
         rmatz(1,2) = 1d0
         rmatz(2,1) = -1d0
         rmatz(2,2) = 0d0
      else
         phiz = datan(a(2)/a(1))
         cphiz = cos(phiz)
         sphiz = sin(phiz)
         rmatz(1,1) = cphiz
         rmatz(1,2) = sphiz
         rmatz(2,1) = -sphiz
         rmatz(2,2) = cphiz
      endif
      rmatz(3,3) = 1d0
      do i=1,2
         rmatz(i,3) = 0d0
         rmatz(3,i) = 0d0
      enddo

      do i=1,4
         b(i)=0d0
         do j=1,4
            b(i)=b(i) + rmatz(i,j)*a(j)
         enddo
      enddo
      
      rmaty(4,4) = 1d0
      do i=1,3
         rmaty(4,i) = 0d0
         rmaty(i,4) = 0d0
      enddo

      if(b(3).eq.0d0) then
         b00 = 1d0
         if (b(1).le.0d0) b00 = -1d0 
         rmaty(1,1) = 0d0
         rmaty(1,3) = -b00
         rmaty(3,1) = b00
         rmaty(3,3) = 0d0
      else
         phiy = datan(b(1)/b(3))
         if (b(3).lt.0d0) phiy = phiy+pi
         cphiy = cos(phiy)
         sphiy = sin(phiy)
         rmaty(1,1) = cphiy
         rmaty(1,3) = -sphiy
         rmaty(3,1) = sphiy
         rmaty(3,3) = cphiy
      endif
      rmaty(2,2) = 1d0
      do i=1,3,2
         rmaty(i,2) = 0d0
         rmaty(2,i) = 0d0
      enddo
      
      do i=1,4
         do j=1,4
            rmat(i,j) = 0d0
            do k=1,4
               rmat(i,j) = rmat(i,j) + rmaty(i,k)*rmatz(k,j)
            enddo
         enddo
      enddo

      return
      end


      subroutine unrotatetoz(a,rmat)
      implicit real*8(a-h,o-z)
      dimension a(4),b(4),rmat(4,4),rmatinv(4,4)
c---- computes rotation from a frame where a is along z 
c---- to frame where a is at its actual value

      call rotatetoz(a,rmatinv)
      do i=1,4
         do j=1,4
            rmat(i,j)=rmatinv(j,i)
         enddo
      enddo
      
      return
      end


      subroutine makeps2com(s,costh,phi,rm1sq,rm2sq,p1,p2,p1a,p2a)
      implicit real*8(a-h,o-z)
      dimension p1(4),p2(4),p1a(4),p2a(4)
      parameter(pi=3.141592653589793238d0)
      
c ---- generates momenta p1,p2 in their centre-of-mass frame
c ---- and momenta p1a,p2a rotated by pi/2 in azimuth
      
      sinth = dsqrt(1d0-costh**2)
      sp = sin(phi)
      cp = cos(phi)
      spa = sin(phi+pi/2d0)
      cpa = cos(phi+pi/2d0)

      smin = rm1sq + rm2sq
      if (s.lt.smin) then
         write(6,*) 'improper call of makeps2com'
         stop
      endif

      roots = dsqrt(s)
      
      e1 = (s + rm1sq - rm2sq)/2d0/roots
      e2 = (s + rm2sq - rm1sq)/2d0/roots
      pp = dsqrt(e1**2-rm1sq)
      
      p1(4) = e1
      p1(1) = pp*sinth*cp
      p1(2) = pp*sinth*sp
      p1(3) = pp*costh

      p1a(4) = e1
      p1a(1) = pp*sinth*cpa
      p1a(2) = pp*sinth*spa
      p1a(3) = pp*costh
      
      p2(4) = e2
      p2a(4) = e2

      do i=1,3
         p2(i) = -p1(i)
         p2a(i) = -p1a(i)
      enddo
      return
      end

      subroutine p3generic(r1,r2,s,rm1sq,rm2sq,rm3sq,p1,p2,p3,wt)
c--- generates generic 3-parton decay phase space in the rest system
c--- for arbitrary time-like external masses
c--- input: random number r1,r2
c--- input: s, external masses
c--- output: p1(4),p2(4),p3(4)
c--- processed: phase space weight wt
c--- phase space dR3 = d^3p1/(2*E1*(2*pi)^3)*[p1<->p2]*[p1<->p3]*(2*pi)^4*del(..)
      implicit real*8(a-h,o-z)
      parameter(pi=3.141592653589793238d0)
      dimension p1(4),p2(4),p3(4),q(4)
      common /tcuts/ymin,y0
      
      roots=dsqrt(s)
      q(1) = 0d0
      q(2) = 0d0
      q(3) = 0d0
      q(4) = roots
      rm1 = dsqrt(rm1sq)
      rm2 = dsqrt(rm2sq)
      rm3 = dsqrt(rm3sq)

      s23min = max(y0*s,(rm2+rm3)**2)
      s23max = (roots-rm1)**2
      if (s23max.le.s23min) then
         wtps = 0d0
         return
      endif
      call pick(1,s23,s23min,s23max,r1,wt)
      disc12 = dsqrt(dlambda(s23,s,rm1sq)*dlambda(s23,rm2sq,rm3sq))
      base12 = rm1sq+rm2sq-0.5d0/s23*(s23-s+rm1sq)*(s23+rm2sq-rm3sq)
      s12min = base12 - 0.5d0/s23*disc12
      s12max = base12 + 0.5d0/s23*disc12
      call pick(1,s12,s12min,s12max,r2,wt)
      if (s12max.le.s12min) then
         wtps = 0d0
         return
      endif
      s13 = s + rm1sq + rm2sq + rm3sq - s12 - s23
      wt = wt/128d0/pi**3/s
      
      e1 = (s+rm1sq-s23)/2d0/roots
      e2 = (s+rm2sq-s13)/2d0/roots
      e3 = (s+rm3sq-s12)/2d0/roots

      p1a = dsqrt(e1**2-rm1sq)
      p2a = dsqrt(e2**2-rm2sq)
      p3a = dsqrt(e3**2-rm3sq)
      
      ct12 = (s12-rm1sq-rm2sq-2d0*e1*e2)/(-2d0*p1a*p2a)
      if (ct12.lt.-1d0.or.ct12.gt.1d0) then
         wtps = 0d0
         return
      endif      
      st12 = dsqrt(1d0-ct12**2)

      p1(4) = e1
      p1(1) = p1a
      p1(2) = 0d0
      p1(3) = 0d0

      p2(4) = e2
      p2(1) = p2a*ct12
      p2(2) = p2a*st12
      p2(3) = 0d0
      
      do i=1,4
         p3(i) = q(i) - p1(i) - p2(i) 
      enddo
      return
      end
      
      real*8 function dlambda(a,b,c)
      implicit real*8(a-h,o-z)
      
      if (dabs(c).lt.min(dabs(a),dabs(b))) then
         dlambda = (a-b)**2 + c*c - 2d0*(a+b)*c
      elseif(dabs(b).lt.min(dabs(a),dabs(c))) then
         dlambda = (a-c)**2 + b*b - 2d0*(a+c)*b
      else
         dlambda = (b-c)**2 + a*a - 2d0*(b+c)*a
      endif
c      dlambda = a*a + b*b + c*c - 2d0*a*b - 2d0*a*c - 2d0*b*c
      return
      end


      subroutine fillinv(n,p,sij)
      implicit real*8(a-h,o-z)
      dimension p(4,1:n),sij(1:n,1:n)
      
      do i=1,n
         sij(i,i) = 0d0
         do j=i+1,n
            sij(i,j) = 2d0*dot(p(1,i),p(1,j))
            sij(j,i) = sij(i,j)
         enddo
      enddo
      return
      end
