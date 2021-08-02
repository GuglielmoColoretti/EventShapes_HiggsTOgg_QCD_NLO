c
c real dilogarithm function  0 <= x <= 1
c
      function rli2(x)
      implicit real*8(a-h,o-z)
      parameter(a1 = -0.250000000000000d0)
      parameter(a2 = -0.111111111111111d0)
      parameter(a3 = -0.010000000000000d0)
      parameter(a4 = -0.0170068027210884d0)
      parameter(a5 = -0.0194444444444444d0)
      parameter(a6 = -0.0206611570247934d0)
      parameter(a7 = -0.0214173006480699d0)
      parameter(a8 = -0.02194886637723d0)
      parameter(a9 = -0.0220893589994137d0)
      parameter(a10 = -0.0229303207720760d0)
      parameter(zeta2 =  1.644934066848226d0)
      
      if(x.gt.1.d0)then
        write(*,*)' argument greater than 1 passed to li2'
        rli2=0.d0
        return
      elseif(x.lt.0.d0)then
        write(*,*)' argument less than 0 passed to li2'
        rli2=0.d0
        return
      elseif(x.eq.1.d0)then
        rli2=zeta2
        return
      elseif(x.eq.0.d0)then
        rli2=0d0
        return
      elseif(x.gt.0.5d0)then
        y=1.d0-x
        z=-log(1.d0-y)
        z2=z*z
        rli2=-z*(1.d0+a1*z*(1.d0+a2*z*(1.d0+a3*z2*(1.d0+a4*z2*
     1 (1.d0+a5*z2*(1.d0+a6*z2*(1.d0+a7*z2*(1.d0+a8*z2*(1.d0+a9*z2*
     2 (1.d0+a10*z2))))))))))
     3 +zeta2-log(x)*log(1.d0-x)
       return
      elseif(x.le.0.5d0)then
        y=x
        z=-log(1.d0-y)
        z2=z*z
        rli2=z*(1.d0+a1*z*(1.d0+a2*z*(1.d0+a3*z2*(1.d0+a4*z2*
     1 (1.d0+a5*z2*(1.d0+a6*z2*(1.d0+a7*z2*(1.d0+a8*z2*(1.d0+a9*z2*
     2 (1.d0+a10*z2))))))))))
        return
      endif
      rli2=0d0
      return
      end
************************************************************************
*

      SUBROUTINE vegas5(ISTAT,FXN,ip,AVGI,SD,CHI2A)
C
C   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
C      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
C
      IMPLICIT REAL*8(A-H,O-Z)
      external FXN
      character char*3,ctype*2,fname*20,gridfile*16
      parameter(ipmx=2)
      COMMON/BVEG5/NDIM(ipmx),NCALL(ipmx),NPRN
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common/outfile/fname
      DIMENSION XI(ipmx,50,10),D(50,10),DI(50,10),XIN(50)
     1   ,R(50),DT(10),X(10),KG(10),IA(10)
      dimension si(ipmx),si2(ipmx),swgt(ipmx),schi(ipmx),calls(ipmx)
     1   ,DXG(ipmx),DV2G(ipmx),XND(ipmx),XJAC(ipmx)
*ng
      dimension GSI(ipmx),GSI2(ipmx),GSWGT(ipmx) 
*ng
      dimension nd(ipmx),ng(ipmx),npg(ipmx),it(ipmx),mds(ipmx)
     1   ,ndo(ipmx),ndm(ipmx)
      REAL*8 QRAN(10)
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/ 
      if(iaver.eq.0)ctype='.A'
      if(iaver.eq.1)ctype='.W'
      if(iaver.eq.2)ctype='.C'
      if(iaver.eq.3)ctype='.M'
      if(iaver.eq.4)ctype='.T'
      if(iaver.eq.5)ctype='.B'
      if(iaver.eq.6)ctype='.D'
      if(iaver.eq.7)ctype='.J'
      if(iaver.eq.8)ctype='.L'
      if(ip.eq. 1)char='v5a'
      if(ip.eq. 2)char='v5b'
      gridfile = 'E'//fname(4:13)//char//ctype
C
      if(istat.eq.0.or.istat.eq.1.or.istat.eq.2)then
c
c         initialize cumulative variables 
c
        IT(ip)=0
        SI(ip)  =0d0
        SI2(ip) =0d0
        SWGT(ip)=0d0
        SCHI(ip)=0d0
*ng
        GSI(ip)  =0d0
        GSI2(ip) =0d0
        GSWGT(ip)=0d0
*mg
        ND(ip)=NDMX
        NG(ip)=1
        MDS(ip)=1
        IF(MDS(ip).ne.0)then 
          NG(ip)=(NCALL(ip)/2d0)**(1d0/NDIM(ip))
          MDS(ip)=1
          IF((2*NG(ip)-NDMX).ge.0)then
            MDS(ip)=-1
            NPG(ip)=NG(ip)/NDMX+1
            ND(ip)=NG(ip)/NPG(ip)
            NG(ip)=NPG(ip)*ND(ip)
          endif
        endif
        K=NG(ip)**NDIM(ip)
        NPG(ip)=NCALL(ip)/K
        IF(NPG(ip).LT.2) NPG(ip)=2
        CALLS(ip)=NPG(ip)*K
        DXG(ip)=ONE/NG(ip)
        DV2G(ip)=(CALLS(ip)*DXG(ip)**NDIM(ip))**2
     .         /NPG(ip)/NPG(ip)/(NPG(ip)-ONE)
        XND(ip)=ND(ip)
        NDM(ip)=ND(ip)-1
        DXG(ip)=DXG(ip)*XND(ip)
        XJAC(ip)=ONE/CALLS(ip)
      IF(NPRN.NE.0.and.istat.ne.0) WRITE(6,200) ip,NDIM(ip),CALLS(ip) 

      endif
      if(istat.eq.0)then
C
C   read in grid   
C
        open(unit=3,file=gridfile,status='unknown')
        write(*,*)'* reading in vegas grid from ',gridfile 
        do j=1,ndim(ip)
          read(3,*) jj,(xi(ip,i,j),i=1,nd(ip))
        enddo
        close(3)
        NDO(ip)=ND(ip)
        return
      elseif(istat.eq.1)then
C
C   construct uniform grid  
C
        RC=1d0/XND(ip)
        DO  J=1,NDIM(ip)
          XI(ip,1,j)=1d0
          K=0
          XN=0d0
          DR=0d0
          I=0
4         K=K+1
          DR=DR+ONE
          XO=XN
          XN=XI(ip,K,J)
5         IF(RC.GT.DR) GO TO 4
          I=I+1
          DR=DR-RC
          XIN(I)=XN-(XN-XO)*DR
          IF(I.LT.NDM(ip)) GO TO 5
          DO  I=1,NDM(ip)
            XI(ip,I,J)=XIN(I)
          enddo
          XI(ip,ND(ip),J)=ONE
        enddo
        NDO(ip)=ND(ip)
        return
C
      elseif(istat.eq.2)then
C
C   rescale refined grid to new ND value - preserve bin density
C
        if(nd(ip).ne.ndo(ip))then
          RC=NDO(ip)/XND(ip)
          DO  J=1,NDIM(ip)
            K=0
            XN=0d0
            DR=0d0
            I=0
6           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(ip,K,J)
7           IF(RC.GT.DR) GO TO 6
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM(ip)) GO TO 7
            DO  I=1,NDM(ip)
              XI(ip,I,J)=XIN(I)
            enddo
            XI(ip,ND(ip),J)=ONE
          enddo
        endif
C
        return
C
      elseif(istat.eq.3.or.istat.eq.4)then
c
c    main integration loop
c         
        IT(ip)=IT(ip)+1
        TI =0d0
        TSI=0d0
*ng
        GTI =0d0
        GTSI=0d0
*ng
        DO J=1,NDIM(ip)
          KG(J)=1
          DO I=1,ND(ip)
           D (I,J)=0d0
           DI(I,J)=0d0
          enddo
        enddo
C
11      FB=0d0
        F2B=0d0
        K=0
12      K=K+1
        do j=1,ndim(ip)
          qran(j)=rn(1)
        enddo        
        WGT=XJAC(ip)
        DO J=1,NDIM(ip)
          XN=(KG(J)-QRAN(J))*DXG(ip)+ONE
          IA(J)=XN
          IF(IA(J).GT.1)then
            XO=XI(ip,IA(J),J)-XI(ip,IA(J)-1,J)
            RC=XI(ip,IA(J)-1,J)+(XN-IA(J))*XO
          else
            XO=XI(ip,IA(J),J)
            RC=(XN-IA(J))*XO
          endif
          X(J)=RC
          WGT=WGT*XO*XND(ip)
        enddo
C
        F=WGT
        F=F*FXN(X,WGT)
        F2=F*F
        FB=FB+F
        F2B=F2B+F2
        DO J=1,NDIM(ip)
          DI(IA(J),J)=DI(IA(J),J)+F
          IF(MDS(ip).GE.0) D(IA(J),J)=D(IA(J),J)+F2
        enddo
        IF(K.LT.NPG(ip)) GO TO 12
C
*ng

        GTI=GTI+FB
        GTSI=GTSI+F2B
*ng

        F2B=DSQRT(F2B*NPG(ip))
        F2B=(F2B-FB)*(F2B+FB)
        TI=TI+FB
        TSI=TSI+F2B
        IF(MDS(ip).lt.0) then
          DO J=1,NDIM(ip)
            D(IA(J),J)=D(IA(J),J)+F2B
          enddo
        endif
        K=NDIM(ip)
19      KG(K)=MOD(KG(K),NG(ip))+1
        IF(KG(K).NE.1) GO TO 11
        K=K-1
        IF(K.GT.0) GO TO 19
C
C   FINAL RESULTS for THIS ITERATION
C
        TSI=TSI*DV2G(ip)
        TI2=TI*TI
        WGT=TI2/TSI
        SI(ip)=SI(ip)+TI*WGT
        SI2(ip)=SI2(ip)+TI2
        SWGT(ip)=SWGT(ip)+WGT
        SCHI(ip)=SCHI(ip)+TI2*WGT
        AVGI=SI(ip)/SWGT(ip)
        SD=SWGT(ip)*IT(ip)/SI2(ip)
        CHI2A=SD*(SCHI(ip)/SWGT(ip)-AVGI*AVGI)/(IT(ip)-.999d0)
        SD=DSQRT(ONE/SD)
*ng

        GTSI=(GTSI*CALLS(ip)-GTI**2)/CALLS(ip)
        GSI(ip)=GSI(ip)+GTI/GTSI
        GSWGT(ip)=GSWGT(ip)+1d0/GTSI
        GAVGI=GSI(ip)/GSWGT(ip)
        GSD=DSQRT(ONE/GSWGT(ip))
*ng
C
        IF(NPRN.ne.0) then
          TSI=DSQRT(TSI)
          WRITE(6,201) IT(ip),ip,TI,TSI,AVGI,SD,CHI2A
*ng
          GTSI=DSQRT(GTSI)
          WRITE(6,201) IT(ip),ip,GTI,GTSI,GAVGI,GSD 
          AVGI=GAVGI
	  SD=GSD
*ng
        endif
      endif
C
C   REFINE GRID
C
      if(istat.eq.3)then
       DO J=1,NDIM(ip)
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2d0
         DT(J)=D(1,J)
         DO I=2,NDM(ip)
           D(I,J)=XO+XN
           XO=XN
           XN=D(I+1,J)
           D(I,J)=(D(I,J)+XN)/3d0
           DT(J)=DT(J)+D(I,J)
         enddo
         D(ND(ip),J)=(XN+XO)/2d0
         DT(J)=DT(J)+D(ND(ip),J)
       enddo
C
        DO 28 J=1,NDIM(ip)
        RC=0d0
        DO 24 I=1,ND(ip)
        R(I)=0d0
        IF(D(I,J).LE.0d0) GO TO 24
        XO=DT(J)/D(I,J)
        R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
24      RC=RC+R(I)
        RC=RC/XND(ip)
        K=0
        XN=0d0
        DR=XN
        I=K
25      K=K+1
        DR=DR+R(K)
        XO=XN
        XN=XI(ip,K,J)
26      IF(RC.GT.DR) GO TO 25
        I=I+1
        DR=DR-RC
        XIN(I)=XN-(XN-XO)*DR/R(K)
        IF(I.LT.NDM(ip)) GO TO 26
        DO 27 I=1,NDM(ip)
27      XI(ip,I,J)=XIN(I)
28      XI(ip,ND(ip),J)=ONE
        open(unit=2,file=gridfile,status='unknown')
        write(*,*)'* writing vegas grid to ',gridfile 
        do j=1,ndim(ip)
          write(2,*) j,(xi(ip,i,j),i=1,nd(ip))
        enddo
        close(2)
        return
      endif
C
200   FORMAT(' * Input parameters for vegas',i2,': ndim=',I2,
     1   ',  nshot= ', F10.0,'  *')
201   FORMAT(i4,'(',i2,') ',g15.7,g13.6,g15.7,g13.6,f7.2)
202   FORMAT(' DATA for AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
      END
*
*
************************************************************************
*
      SUBROUTINE vegas4a(ISTAT,FXN,ip,AVGI,SD,CHI2A)
C
C   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
C      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
C
      IMPLICIT REAL*8(A-H,O-Z)
      external FXN
      parameter(ipmx=1)
      COMMON/BVEG4/NDIM(ipmx),NCALL(ipmx),NPRN
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common/outfile/fname
      DIMENSION XI(ipmx,50,10),D(50,10),DI(50,10),XIN(50)
     1   ,R(50),DT(10),X(10),KG(10),IA(10)
      dimension si(ipmx),si2(ipmx),swgt(ipmx),schi(ipmx),calls(ipmx)
     1   ,DXG(ipmx),DV2G(ipmx),XND(ipmx),XJAC(ipmx)
*ng
      dimension GSI(ipmx),GSI2(ipmx),GSWGT(ipmx) 
*ng
      dimension nd(ipmx),ng(ipmx),npg(ipmx),it(ipmx),mds(ipmx)
     1   ,ndo(ipmx),ndm(ipmx)
      REAL*8 QRAN(10)
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/ 
      character char*3,ctype*2,fname*20,gridfile*16
      if(iaver.eq.0)ctype='.A'
      if(iaver.eq.1)ctype='.W'
      if(iaver.eq.2)ctype='.C'
      if(iaver.eq.3)ctype='.M'
      if(iaver.eq.4)ctype='.T'
      if(iaver.eq.5)ctype='.B'
      if(iaver.eq.6)ctype='.D'
      if(iaver.eq.7)ctype='.J'
      if(iaver.eq.8)ctype='.L'
      if(ip.eq. 1)char='v4a'
      gridfile = 'E'//fname(4:13)//char//ctype

C
      if(istat.eq.0.or.istat.eq.1.or.istat.eq.2)then
c
c         initialize cumulative variables 
c
        IT(ip)=0
        SI(ip)  =0d0
        SI2(ip) =0d0
        SWGT(ip)=0d0
        SCHI(ip)=0d0
*ng
        GSI(ip)  =0d0
        GSI2(ip) =0d0
        GSWGT(ip)=0d0
*mg
        ND(ip)=NDMX
        NG(ip)=1
        MDS(ip)=1
        IF(MDS(ip).ne.0)then 
          NG(ip)=(NCALL(ip)/2d0)**(1d0/NDIM(ip))
          MDS(ip)=1
          IF((2*NG(ip)-NDMX).ge.0)then
            MDS(ip)=-1
            NPG(ip)=NG(ip)/NDMX+1
            ND(ip)=NG(ip)/NPG(ip)
            NG(ip)=NPG(ip)*ND(ip)
          endif
        endif
        K=NG(ip)**NDIM(ip)
        NPG(ip)=NCALL(ip)/K
        IF(NPG(ip).LT.2) NPG(ip)=2
        CALLS(ip)=NPG(ip)*K
        DXG(ip)=ONE/NG(ip)
        DV2G(ip)=(CALLS(ip)*DXG(ip)**NDIM(ip))**2
     .         /NPG(ip)/NPG(ip)/(NPG(ip)-ONE)
        XND(ip)=ND(ip)
        NDM(ip)=ND(ip)-1
        DXG(ip)=DXG(ip)*XND(ip)
        XJAC(ip)=ONE/CALLS(ip)
      IF(NPRN.NE.0.and.istat.ne.0) WRITE(6,200) ip,NDIM(ip),CALLS(ip) 

      endif
      if(istat.eq.0)then
C
C   read in grid   
C
        open(unit=3,file=gridfile,status='unknown')
        write(*,*)'* reading in vegas grid from ',gridfile 
         do j=1,ndim(ip)
          read(3,*) jj,(xi(ip,i,j),i=1,nd(ip))
        enddo
        close(3)
        NDO(ip)=ND(ip)
        return
      elseif(istat.eq.1)then
C
C   construct uniform grid  
C
        RC=1d0/XND(ip)
        DO  J=1,NDIM(ip)
          XI(ip,1,j)=1d0
          K=0
          XN=0d0
          DR=0d0
          I=0
4         K=K+1
          DR=DR+ONE
          XO=XN
          XN=XI(ip,K,J)
5         IF(RC.GT.DR) GO TO 4
          I=I+1
          DR=DR-RC
          XIN(I)=XN-(XN-XO)*DR
          IF(I.LT.NDM(ip)) GO TO 5
          DO  I=1,NDM(ip)
            XI(ip,I,J)=XIN(I)
          enddo
          XI(ip,ND(ip),J)=ONE
        enddo
        NDO(ip)=ND(ip)
        return
C
      elseif(istat.eq.2)then
C
C   rescale refined grid to new ND value - preserve bin density
C
        if(nd(ip).ne.ndo(ip))then
          RC=NDO(ip)/XND(ip)
          DO  J=1,NDIM(ip)
            K=0
            XN=0d0
            DR=0d0
            I=0
6           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(ip,K,J)
7           IF(RC.GT.DR) GO TO 6
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM(ip)) GO TO 7
            DO  I=1,NDM(ip)
              XI(ip,I,J)=XIN(I)
            enddo
            XI(ip,ND(ip),J)=ONE
          enddo
        endif
        return
C
      elseif(istat.eq.3.or.istat.eq.4)then
c
c    main integration loop
c         
        IT(ip)=IT(ip)+1
        TI =0d0
        TSI=0d0
*ng
        GTI =0d0
        GTSI=0d0
*ng
        DO J=1,NDIM(ip)
          KG(J)=1
          DO I=1,ND(ip)
           D (I,J)=0d0
           DI(I,J)=0d0
          enddo
        enddo
C
11      FB=0d0
        F2B=0d0
        K=0
12      K=K+1
        do j=1,ndim(ip)
          qran(j)=rn(1)
        enddo        
        WGT=XJAC(ip)
        DO J=1,NDIM(ip)
          XN=(KG(J)-QRAN(J))*DXG(ip)+ONE
          IA(J)=XN
          IF(IA(J).GT.1)then
            XO=XI(ip,IA(J),J)-XI(ip,IA(J)-1,J)
            RC=XI(ip,IA(J)-1,J)+(XN-IA(J))*XO
          else
            XO=XI(ip,IA(J),J)
            RC=(XN-IA(J))*XO
          endif
          X(J)=RC
          WGT=WGT*XO*XND(ip)
        enddo
C
        F=WGT
        F=F*FXN(X,WGT)
        F2=F*F
        FB=FB+F
        F2B=F2B+F2
        DO J=1,NDIM(ip)
          DI(IA(J),J)=DI(IA(J),J)+F
          IF(MDS(ip).GE.0) D(IA(J),J)=D(IA(J),J)+F2
        enddo
        IF(K.LT.NPG(ip)) GO TO 12
C
*ng

        GTI=GTI+FB
        GTSI=GTSI+F2B
*ng
        F2B=DSQRT(F2B*NPG(ip))
        F2B=(F2B-FB)*(F2B+FB)
        TI=TI+FB
        TSI=TSI+F2B
        IF(MDS(ip).lt.0) then
          DO J=1,NDIM(ip)
            D(IA(J),J)=D(IA(J),J)+F2B
          enddo
        endif
        K=NDIM(ip)
19      KG(K)=MOD(KG(K),NG(ip))+1
        IF(KG(K).NE.1) GO TO 11
        K=K-1
        IF(K.GT.0) GO TO 19
C
C   FINAL RESULTS for THIS ITERATION
C
        TSI=TSI*DV2G(ip)
        TI2=TI*TI
        WGT=TI2/TSI
        SI(ip)=SI(ip)+TI*WGT
        SI2(ip)=SI2(ip)+TI2
        SWGT(ip)=SWGT(ip)+WGT
        SCHI(ip)=SCHI(ip)+TI2*WGT
        AVGI=SI(ip)/SWGT(ip)
        SD=SWGT(ip)*IT(ip)/SI2(ip)
        CHI2A=SD*(SCHI(ip)/SWGT(ip)-AVGI*AVGI)/(IT(ip)-.999d0)
        SD=DSQRT(ONE/SD)
*ng

        GTSI=(GTSI*CALLS(ip)-GTI**2)/CALLS(ip)
        GSI(ip)=GSI(ip)+GTI/GTSI
        GSWGT(ip)=GSWGT(ip)+1d0/GTSI
        GAVGI=GSI(ip)/GSWGT(ip)
        GSD=DSQRT(ONE/GSWGT(ip))
*ng
C
        IF(NPRN.ne.0) then
          TSI=DSQRT(TSI)
          WRITE(6,201) IT(ip),ip,TI,TSI,AVGI,SD,CHI2A
*ng
          GTSI=DSQRT(GTSI)
          WRITE(6,201) IT(ip),ip,GTI,GTSI,GAVGI,GSD 
          AVGI=GAVGI
	  SD=GSD
*ng
        endif
      endif
C
C   REFINE GRID
C
      if(istat.eq.3)then
       DO J=1,NDIM(ip)
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2d0
         DT(J)=D(1,J)
         DO I=2,NDM(ip)
           D(I,J)=XO+XN
           XO=XN
           XN=D(I+1,J)
           D(I,J)=(D(I,J)+XN)/3d0
           DT(J)=DT(J)+D(I,J)
         enddo
         D(ND(ip),J)=(XN+XO)/2d0
         DT(J)=DT(J)+D(ND(ip),J)
       enddo
C
        DO 28 J=1,NDIM(ip)
        RC=0d0
        DO 24 I=1,ND(ip)
        R(I)=0d0
        IF(D(I,J).LE.0d0) GO TO 24
        XO=DT(J)/D(I,J)
        R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
24      RC=RC+R(I)
        RC=RC/XND(ip)
        K=0
        XN=0d0
        DR=XN
        I=K
25      K=K+1
        DR=DR+R(K)
        XO=XN
        XN=XI(ip,K,J)
26      IF(RC.GT.DR) GO TO 25
        I=I+1
        DR=DR-RC
        XIN(I)=XN-(XN-XO)*DR/R(K)
        IF(I.LT.NDM(ip)) GO TO 26
        DO 27 I=1,NDM(ip)
27      XI(ip,I,J)=XIN(I)
28      XI(ip,ND(ip),J)=ONE
        open(unit=2,file=gridfile,status='unknown')
        write(*,*)'* writing vegas grid to ',gridfile 
        do j=1,ndim(ip)
          write(2,*) j,(xi(ip,i,j),i=1,nd(ip))
        enddo
        close(2)
        return
      endif
C
200   FORMAT(' * Input parameters for vegas',i2,': ndim=',I2,
     1   ',  nshot= ', F10.0,'  *')
201   FORMAT(i4,'(',i2,') ',g15.7,g13.6,g15.7,g13.6,f7.2)
202   FORMAT(' DATA for AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
      END
*
*
************************************************************************
*


      SUBROUTINE vegas3(ISTAT,FXN,AVGI,SD,CHI2A)
C
C   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
C      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
C
      IMPLICIT REAL*8(A-H,O-Z)
      external FXN
      COMMON/BVEG3/NDIM,NCALL,NPRN
      common/intech/iaver,imom,idist,iang,idebug
      common/inphys/nloop,icol,njets
      common/outfile/fname
      DIMENSION XI(50,10),D(50,10),DI(50,10),XIN(50),R(50),DT(10),X(10)
     1   ,KG(10),IA(10)
      REAL*8 QRAN(10)
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/,MDS/1/
      character char*3,ctype*2,fname*20,gridfile*16
      if(iaver.eq.0)ctype='.A'
      if(iaver.eq.1)ctype='.W'
      if(iaver.eq.2)ctype='.C'
      if(iaver.eq.3)ctype='.M'
      if(iaver.eq.4)ctype='.T'
      if(iaver.eq.5)ctype='.B'
      if(iaver.eq.6)ctype='.D'
      if(iaver.eq.7)ctype='.J'
      if(iaver.eq.8)ctype='.L'
      char='v3a'
      gridfile = 'E'//fname(4:13)//char//ctype
C
      if(istat.eq.0.or.istat.eq.1.or.istat.eq.2)then
c
c         initialize cumulative variables 
c
        IT=0
        SI  =0d0
        SI2 =0d0
        SWGT=0d0
        SCHI=0d0
*ng
        GSI   =0d0
        GSI2  =0d0
        GSWGT =0d0
*mg
        ND=NDMX
        NG=1
        IF(MDS.ne.0)then 
          NG=(NCALL/2d0)**(1d0/NDIM)
          MDS=1
          IF((2*NG-NDMX).ge.0)then
            MDS=-1
            NPG=NG/NDMX+1
            ND=NG/NPG
            NG=NPG*ND
          endif
        endif
        K=NG**NDIM
        NPG=NCALL/K
        IF(NPG.LT.2) NPG=2
        CALLS=NPG*K
        DXG=ONE/NG
        DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)
        XND=ND
        NDM=ND-1
        DXG=DXG*XND
        XJAC=ONE/CALLS
      IF(NPRN.NE.0.and.istat.ne.0) WRITE(6,200) NDIM,CALLS 
      endif
      if(istat.eq.0)then
C
C   read in grid   
C
        open(unit=3,file=gridfile,status='unknown')
        write(*,*)'* reading in vegas grid from ',gridfile 
        do j=1,ndim
          read(3,*) jj,(xi(i,j),i=1,nd)
        enddo
        close(3)
        NDO=ND
        return
      elseif(istat.eq.1)then
C
C   construct uniform grid  
C
        RC=1d0/XND
        DO  J=1,NDIM
          xi(1,j)=1d0
          K=0
          XN=0d0
          DR=0d0
          I=0
4         K=K+1
          DR=DR+ONE
          XO=XN
          XN=XI(K,J)
5         IF(RC.GT.DR) GO TO 4
          I=I+1
          DR=DR-RC
          XIN(I)=XN-(XN-XO)*DR
          IF(I.LT.NDM) GO TO 5
          DO  I=1,NDM
            XI(I,J)=XIN(I)
          enddo
          XI(ND,J)=ONE
        enddo
        NDO=ND
        return
C
      elseif(istat.eq.2)then
C
C   rescale refined grid to new ND value - preserve bin density
C

        if(nd.ne.ndo)then
          RC=NDO/XND
          DO  J=1,NDIM
            K=0
            XN=0d0
            DR=0d0
            I=0
6           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
7           IF(RC.GT.DR) GO TO 6
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM) GO TO 7
            DO  I=1,NDM
              XI(I,J)=XIN(I)
            enddo
            XI(ND,J)=ONE
          enddo
        endif
        return
C
      elseif(istat.eq.3.or.istat.eq.4)then
c
c    main integration loop
c         
        IT=IT+1
        TI =0d0
        TSI=0d0
*ng
        GTI =0d0
        GTSI=0d0
*ng
        DO J=1,NDIM
          KG(J)=1
          DO I=1,ND
           D (I,J)=0d0
           DI(I,J)=0d0
          enddo
        enddo
C
11      FB=0d0
        F2B=0d0
        K=0
12      K=K+1
        do j=1,ndim
          qran(j)=rn(1)
        enddo        
        WGT=XJAC
        DO J=1,NDIM
          XN=(KG(J)-QRAN(J))*DXG+ONE
          IA(J)=XN
          IF(IA(J).GT.1)then
            XO=XI(IA(J),J)-XI(IA(J)-1,J)
            RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
          else
            XO=XI(IA(J),J)
            RC=(XN-IA(J))*XO
          endif
          X(J)=RC
          WGT=WGT*XO*XND
        enddo
C
        F=WGT
        F=F*FXN(X,WGT)
        F2=F*F
        FB=FB+F
        F2B=F2B+F2
        DO J=1,NDIM
          DI(IA(J),J)=DI(IA(J),J)+F
          IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
        enddo
        IF(K.LT.NPG) GO TO 12
C
*ng

        GTI=GTI+FB
        GTSI=GTSI+F2B
*ng
        F2B=DSQRT(F2B*NPG)
        F2B=(F2B-FB)*(F2B+FB)
        TI=TI+FB
        TSI=TSI+F2B
        IF(MDS.lt.0) then
          DO J=1,NDIM
            D(IA(J),J)=D(IA(J),J)+F2B
          enddo
        endif
        K=NDIM
19      KG(K)=MOD(KG(K),NG)+1
        IF(KG(K).NE.1) GO TO 11
        K=K-1
        IF(K.GT.0) GO TO 19
C
C   FINAL RESULTS for THIS ITERATION
C
        TSI=TSI*DV2G
        TI2=TI*TI
        WGT=TI2/TSI
        SI=SI+TI*WGT
        SI2=SI2+TI2
        SWGT=SWGT+WGT
        SCHI=SCHI+TI2*WGT
        AVGI=SI/SWGT
        SD=SWGT*IT/SI2
        CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-.999d0)
        SD=DSQRT(ONE/SD)
*ng

        GTSI=(GTSI*CALLS-GTI**2)/CALLS
        GSI=GSI+GTI/GTSI
        GSWGT=GSWGT+1d0/GTSI
        GAVGI=GSI/GSWGT
        GSD=DSQRT(ONE/GSWGT)
*ng
C
        IF(NPRN.ne.0) then
          TSI=DSQRT(TSI)
          WRITE(6,201) IT,TI,TSI,AVGI,SD,CHI2A
*ng
          GTSI=DSQRT(GTSI)
          WRITE(6,201) IT,GTI,GTSI,GAVGI,GSD 
          AVGI=GAVGI
	  SD=GSD
*ng
        endif
        IF(NPRN.le.0) then
          DO J=1,NDIM
            WRITE(6,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
          enddo
        endif
      endif
C
C   REFINE GRID
C
      if(istat.eq.3)then
       DO J=1,NDIM
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2d0
         DT(J)=D(1,J)
         DO I=2,NDM
           D(I,J)=XO+XN
           XO=XN
           XN=D(I+1,J)
           D(I,J)=(D(I,J)+XN)/3d0
           DT(J)=DT(J)+D(I,J)
        enddo
        D(ND,J)=(XN+XO)/2d0
      DT(J)=DT(J)+D(ND,J)
      enddo
C
        DO 28 J=1,NDIM
        RC=0d0
        DO 24 I=1,ND
        R(I)=0d0
        IF(D(I,J).LE.0d0) GO TO 24
        XO=DT(J)/D(I,J)
        R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
24      RC=RC+R(I)
        RC=RC/XND
        K=0
        XN=0d0
        DR=XN
        I=K
25      K=K+1
        DR=DR+R(K)
        XO=XN
        XN=XI(K,J)
26      IF(RC.GT.DR) GO TO 25
        I=I+1
        DR=DR-RC
        XIN(I)=XN-(XN-XO)*DR/R(K)
        IF(I.LT.NDM) GO TO 26
        DO 27 I=1,NDM
27      XI(I,J)=XIN(I)
28      XI(ND,J)=ONE
        open(unit=2,file=gridfile,status='unknown')
        write(*,*)'* writing vegas grid to ',gridfile 
        do j=1,ndim
          write(2,*) j,(xi(i,j),i=1,nd)
        enddo
        close(2)
        return
      endif
C
200   FORMAT(' * Input parameters for vegas 3: ndim=',I2,
     1   ',  nshot= ', F10.0,'  *')
201   FORMAT(i4,'( 3) ',g15.7,g13.6,g15.7,g13.6,f7.2)
202   FORMAT(' DATA for AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
      END


************************************************************************
*
*
* random number generator
*
      function rn(idummy)
      real*8 rn,ran
      common/rseeds/i1,i2
      save init
      data init /1/
c      i1=3823
c      i2=8151
      if (init.eq.1) then
        init=0
*        read(*,*)i1,i2
        write(*,*) '* seeding with (',i1,',',i2,')'
        call rmarin(i1,i2)
      end if
*
  10  call ranmar(ran)
      if (ran.lt.1d-16) goto 10
      rn=ran
*
      end
*
      subroutine ranmar(rvec)
*     -----------------
* universal random number generator proposed by marsaglia and zaman
* in report fsu-scri-87-50
* in this version rvec is a double precision variable.
      implicit real*8(a-h,o-z)
      common/ raset1 / ranu(97),ranc,rancd,rancm
      common/ raset2 / iranmr,jranmr
      save /raset1/,/raset2/
      uni = ranu(iranmr) - ranu(jranmr)
      if(uni .lt. 0d0) uni = uni + 1d0
      ranu(iranmr) = uni
      iranmr = iranmr - 1
      jranmr = jranmr - 1
      if(iranmr .eq. 0) iranmr = 97
      if(jranmr .eq. 0) jranmr = 97
      ranc = ranc - rancd
      if(ranc .lt. 0d0) ranc = ranc + rancm
      uni = uni - ranc
      if(uni .lt. 0d0) uni = uni + 1d0
      rvec = uni
      end
 
      subroutine rmarin(ij,kl)
*     -----------------
* initializing routine for ranmar, must be called before generating
* any pseudorandom numbers with ranmar. the input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      implicit real*8(a-h,o-z)
      common/ raset1 / ranu(97),ranc,rancd,rancm
      common/ raset2 / iranmr,jranmr
      save /raset1/,/raset2/
* this shows correspondence between the simplified input seeds ij, kl
* and the original marsaglia-zaman seeds i,j,k,l.
* to get the standard values in the marsaglia-zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      i = mod( ij/177 , 177 ) + 2
      j = mod( ij     , 177 ) + 2
      k = mod( kl/169 , 178 ) + 1
      l = mod( kl     , 169 )
      do 300 ii = 1 , 97
        s =  0d0
        t = .5d0
        do 200 jj = 1 , 24
          m = mod( mod(i*j,179)*k , 179 )
          i = j
          j = k
          k = m
          l = mod( 53*l+1 , 169 )
          if(mod(l*m,64) .ge. 32) s = s + t
          t = .5d0*t
  200   continue
        ranu(ii) = s
  300 continue
      ranc  =   362436d0 / 16777216d0
      rancd =  7654321d0 / 16777216d0
      rancm = 16777213d0 / 16777216d0
      iranmr = 97
      jranmr = 33
      end
       double precision function rnnew(dummy)
*
*      random number function taken from knuth
*      (seminumerical algorithms).
*      method is x(n)=mod(x(n-55)-x(n-24),1/fmodul)
*      no provision yet for control over the seed number.
*
*      ranf gives one random number between 0 and 1.
*      irn55 generates 55 random numbers between 0 and 1/fmodul.
*      in55  initializes the 55 numbers and warms up the sequence.
*
       implicit double precision (a-h,o-z)
       parameter (fmodul=1.d-09)
       integer ia(55)
       save ia
       data ncall/0/
       data mcall/55/
       if( ncall.eq.0 ) then
           call in55 ( ia,234612947 )
           ncall = 1
       endif
       if ( mcall.eq.0 ) then
           call irn55(ia)
           mcall=55
       endif
       rnnew=ia(mcall)*fmodul
       mcall=mcall-1
       end

       subroutine in55(ia,ix)
       parameter (modulo=1000000000)
       integer ia(55)
       ia(55)=ix
       j=ix
       k=1
       do 10 i=1,54
       ii=mod(21*i,55)
       ia(ii)=k
       k=j-k
       if(k.lt.0)k=k+modulo
       j=ia(ii)
   10  continue
       do 20 i=1,10
       call irn55(ia)
   20  continue
       end

       subroutine irn55(ia)
       parameter (modulo=1000000000)
       integer ia(55)
       do 10 i=1,24
       j=ia(i)-ia(i+31)
       if(j.lt.0)j=j+modulo
       ia(i)=j
   10  continue
       do 20 i=25,55
       j=ia(i)-ia(i-24)
       if(j.lt.0)j=j+modulo
       ia(i)=j
   20  continue
       end

*
* The histo-handlers, makes easy interface between user routine 'bino'
* and general histogram manipulator 'ghiman'.
* 
* 'histoi' : sets up histogram 'idhis' with minimum bin value 'bmin',
*            maximum binvalue 'bmax' and number of bins 'nbin'
* 'histoa' : make specific entry in histogram 'idhis' for value 'val' and 
*            weight 'wgt'
* 'histoe' : calculate error request for histogram 'idhis', pipe through
* 'histow' : output request for histogram 'idhis', pipe through
*
      subroutine histoi(idhis,bmin,bmax,nbin)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /hispar/hmin,hwidth,ibin 
      dimension hmin(nhisto),hwidth(nhisto),ibin(nhisto)
*
* histogram initialization
*
      call ghiman(0,idhis,nbin,0d0,0d0)
      hmin(idhis)=bmin
      ibin(idhis)=nbin
      hwidth(idhis)=(bmax-bmin)/nbin
*
      return
      end
*
      subroutine histoa(idhis,val,wgt)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /hispar/hmin,hwidth,ibin 
      dimension hmin(nhisto),hwidth(nhisto),ibin(nhisto)
*
* histogram entry
*
      if(val.lt.hmin(idhis))return
      if(idhis.le.100)then 				
        iloc=1+int((val-hmin(idhis))/hwidth(idhis))
        call ghiman(1,idhis,iloc,wgt,wgt*val)
      endif
*
      return
      end
*
      subroutine histoe(istat,idhis)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /hispar/hmin,hwidth,ibin 
      dimension hmin(nhisto),hwidth(nhisto),ibin(nhisto)
*
* event errors request, pipe through with correct 'ghiman' call
*
*      write(100,*)idhis,hmin(idhis),hwidth(idhis),ibin(idhis)
      call ghiman(istat,idhis,ibin(idhis),0d0,0d0)
*
      return
      end
*
      subroutine histow(idhis)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /hispar/hmin,hwidth,ibin 
      dimension hmin(nhisto),hwidth(nhisto),ibin(nhisto)
*
* output request, pipe through with correct 'ghiman' call
*
      call ghiman(4,idhis,ibin(idhis),hmin(idhis),hwidth(idhis))
*
      return
      end

      subroutine histow1(idhis)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /hispar/hmin,hwidth,ibin 
      dimension hmin(nhisto),hwidth(nhisto),ibin(nhisto)
*
* output request, pipe through with correct 'ghiman' call
*
      call ghiman(5,idhis,ibin(idhis),hmin(idhis),hwidth(idhis))
*
      return
      end

      subroutine histowf(idhis,lun)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /hispar/hmin,hwidth,ibin 
      dimension hmin(nhisto),hwidth(nhisto),ibin(nhisto)
*
* output request, pipe through with correct 'ghiman' call
*
**  write to file, logical unit 'lun'
      call ghiman(lun,idhis,ibin(idhis),hmin(idhis),hwidth(idhis))
*
      return
      end

*
*
* General histogram manipulator, has no knowledge about specific histograms.
*
* Maximum number of histograms is given by the nhisto parameter.
* Maximum number of bins       is given by the maxbin parameter.
*
* istat = 0 : set histogram 'idhis' to value 'entry1' from bin 1 to bin 'iloc'
*             and resets sweep counter 'itmx' to zero.
* istat = 1 : add in histogram 'idhis' weight 'entry1' in  bin location 'iloc'.
* istat = 2 : accumulate event errors in histogram 'idhis' from bin 1 to bin
*             'iloc' for this sweep and increases sweepcounter 'itmx' by 1.
* istat = 3 : calculate standard deviation of the 'itmx' sweeps per bin in 
*             histogram 'idhis' from bin 1 to bin 'iloc' as monte carlo error
*             estimate.
* istat = 4 : write final output to screen in histogram 'idhis' from bin 1 to 
*             bin 'iloc' with offset 'entry1' and binwidth 'entry2'
*            (format: 
*         from bin_number = 1 to 'iloc'  
*            write 'entry1'+'entry2'*(bin_number - 0.5), bin_value, bin_error
*         endfrom) 
*             [for istat>10: write final output to file lun=istat]
* istat = 5 : write final output to screen in histogram 'idhis' from bin 1 to 
*             bin 'iloc' with offset 'entry1' and binwidth 'entry2'
*            (format: 
*         from bin_number = 1 to 'iloc'  
*            write 'entry1'+'entry2'*(bin_number - 1), bin_value, bin_error
*         endfrom)
*
      subroutine ghiman(istat,idhis,iloc,entry1,entry2)
      implicit double precision (a-h,o-z)
      parameter(nhisto=100,maxbin=400)
      common /runinfo/itmax1,itmax2,nshot3,nshot4,nshot5(2) 
      dimension bin(nhisto,4,maxbin),xbin(nhisto,4)
*
      if ((idhis.lt.1).or.(idhis.gt.nhisto)) return
      if ((iloc .lt.1).or.(iloc .gt.maxbin)) return
*
*       init histograms
*
      if (istat.eq.0) then
         do i=1,iloc
            do j=1,4
               bin(idhis,j,i)=entry1
            enddo
         enddo
         itmx=0
      endif
*
* write event into histogram
*
      if (istat.eq.1) then
         if ((iloc.ge.1).and.(iloc.le.maxbin)) then
            bin(idhis,1,iloc)=bin(idhis,1,iloc)+entry1
            bin(idhis,4,iloc)=bin(idhis,4,iloc)+1d0
            xbin(idhis,1)=xbin(idhis,1)+entry2
            xbin(idhis,4)=xbin(idhis,4)+1d0
         endif
      endif
*
* accumulate event errors
*
      if (istat.eq.2) then
         do i=1,iloc
            bin(idhis,2,i)=bin(idhis,2,i)+bin(idhis,1,i)**2
            bin(idhis,3,i)=bin(idhis,3,i)+bin(idhis,1,i)
            bin(idhis,1,i)=0d0
            xbin(idhis,2)=xbin(idhis,2)+xbin(idhis,1)**2
            xbin(idhis,3)=xbin(idhis,3)+xbin(idhis,1)
            xbin(idhis,1)=0d0
         enddo
      endif
*
* calculate event errors
*
      if (istat.eq.3) then
        do i=1,iloc
          bin(idhis,2,i)=
     .  sqrt((bin(idhis,2,i)/itmax2-(bin(idhis,3,i)/itmax2)**2)
     .        /float(itmax2-1))
        enddo
          xbin(idhis,2)=
     .  sqrt((xbin(idhis,2)/itmax2-(xbin(idhis,3)/itmax2)**2)
     .        /float(itmax2-1))
      endif
*
* output distributions
*
      if (istat.eq.4) then
         sum=0d0
         sum2=0d0
         do i=1,iloc
            y=bin(idhis,3,i)/itmax2 
            y2=bin(idhis,2,i)
	    sum=sum+y
	    sum2=sum2+y2**2
            nn=bin(idhis,4,i)
            if(idhis.le.100)then		
              x=entry1+entry2*(dfloat(i)-.5d0)
            write(6,101) x,y/entry2,y2/entry2 
            endif
         enddo
         write(6,102)sum,sqrt(sum2)
         write(6,103)xbin(idhis,3)/itmax2,xbin(idhis,2)
      endif
 101     format(3x,f11.6,1pe12.4,1pe12.4)
 102     format(' sum ',1pe12.4,1pe12.4)
 103     format(' <y> ',1pe12.4,1pe12.4)
      if (istat.ge.11) then
         sum=0d0
         sum2=0d0
         do i=1,iloc
            y=bin(idhis,3,i)/itmax2 
            y2=bin(idhis,2,i)
	    sum=sum+y
	    sum2=sum2+y2**2
            nn=bin(idhis,4,i)
            if(idhis.le.100)then			
              x=entry1+entry2*(dfloat(i)-.5d0)
            write(istat,101) x,y/entry2,y2/entry2 
            endif
         enddo
c         write(11,102)sum,sqrt(sum2)
c         write(11,103)xbin(idhis,3)/itmax2,xbin(idhis,2)
      endif
      if (istat.eq.5) then
         sum=0d0
         sum2=0d0
         do i=1,iloc
            y=bin(idhis,3,i)/itmax2 
            y2=bin(idhis,2,i)
	    sum=sum+y
	    sum2=sum2+y2**2
            nn=bin(idhis,4,i)
            x=entry1+entry2*(dfloat(i)-1d0)
            write(6,101) x,y/entry2,y2/entry2 
         enddo
         write(6,102)sum,sqrt(sum2)
         write(6,103)xbin(idhis,3)/itmax2,xbin(idhis,2)
      endif
*     
      return
      end
*
