C **** Program PME *****************************************************
C *                                                                    *
C *  Sample program for PME Method    Ver.1.0 3/25/02     by Shige.    *
C *                                                                    *
C * CUTEG(cutoff of reciprocal vector), GAM(beta in ref),              *
C * NM1,2,3(mesh size), IBSOD(order of b-spline)                       *
C *  should be chosen carefully                                        *
C **********************************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
C  Ref.1  Essmann et al. J.Chem.Phys 103(19), 8577(1995)
C  Ref.2  Deserno et al. J.Chem.Phys 109(18), 7678(1998)
C
C
      PROGRAM PMESH
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (NM=2300)
C
      COMMON /NUMBE/ NATOM
      COMMON /VALUQ/ QX(NM),QY(NM),QZ(NM),DQX(NM),DQY(NM),DQZ(NM)
      COMMON /UNIV1/ RLV(3,3),RLVI(3,3),RLVD(3,3)
C
      COMMON /PMEPM/ IBSOD,NM1,NM2,NM3,NK1,NK2,NK3
      COMMON /PMEP2/ RCUTPML,RCUTPMS
      COMMON /GAMMA/ GAM
C
C PME PARAMETER
C
c
c  order of B-spline (even number)
      IBSOD=8
c
c  mech size
      NM1=64
      NM2=64
      NM3=64
c
c  n1,n2,n3 in Ref.1 eq.(4.6)
      NK1=1
      NK2=1
      NK3=1
c
c  cutoff radius of list vector
      RCUTPML=14d-10
c
c  cutoff radius for Direct term
      RCUTPMS=9d-10
C
c  make table for FFT
      CALL TABLE(NM1,NM2,NM3)
c
c  some inputs
      CALL INPUTG
c
c  make reciprocal vector set
      CALL GLORDER
c
c  make list vector
      CALL LSVEC
c
c  pre process of PME
      CALL PREPME
c
c  calculate PME energy, force, and stress
      CALL PME(NATOM,QX,QY,QZ,RLV,RLVI)
C
      STOP 'NORMAL END'
      END
C
C **** GLORDER *********************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE GLORDER
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (IGG3=100000)
      COMMON /UNIV1/ RLV(3,3),RLVI(3,3),RLVD(3,3)
      COMMON /RECV1/ G1LVX,G1LVY,G1LVZ
      COMMON /RECV2/ G2LVX,G2LVY,G2LVZ
      COMMON /RECV3/ G3LVX,G3LVY,G3LVZ
      COMMON /CUTEW/ CUTEG,CUTEL
      COMMON /NUMGL/ NUML,NUMG
      COMMON /GGGGG/ GG1(IGG3),GG2(IGG3),GG3(IGG3)
      COMMON /FFTPM/ NG1(IGG3),NG2(IGG3),NG3(IGG3),NNG(IGG3)
      COMMON /PMEPM/ IBSOD,NM1,NM2,NM3,NK1,NK2,NK3
C
c make reciprocal vector set
      NUMG=0
      DO 20 IGX=-NM1/2,NM1/2-1
      DO 20 IGY=-NM2/2,NM2/2-1
      DO 20 IGZ=-NM3/2,NM3/2-1
        DIX=IGX
        DIY=IGY
        DIZ=IGZ
        DKGG1=DIX*G1LVX+DIY*G2LVX+DIZ*G3LVX
        DKGG2=DIX*G1LVY+DIY*G2LVY+DIZ*G3LVY
        DKGG3=DIX*G1LVZ+DIY*G2LVZ+DIZ*G3LVZ
        DKG2=DKGG1*DKGG1+DKGG2*DKGG2+DKGG3*DKGG3
        IF(DKG2.LE.CUTEG**2) THEN
          NUMG=NUMG+1
          NG1(NUMG)=IGX
          NG2(NUMG)=IGY
          NG3(NUMG)=IGZ
        END IF 
  20  CONTINUE
      WRITE(*,*)'CUTEG = ',CUTEG
      WRITE(*,*)'NUMG = ',NUMG
c
c numbering of reciprocal vector
      DO 40 I=1,NUMG
        IF(NG1(I).LT.0) THEN
          NNG1=NG1(I)+NM1 
        ELSE
          NNG1=NG1(I) 
        END IF
        IF(NG2(I).LT.0) THEN
          NNG2=NG2(I)+NM2 
        ELSE
          NNG2=NG2(I) 
        END IF
        IF(NG3(I).LT.0) THEN
          NNG3=NG3(I)+NM3 
        ELSE
          NNG3=NG3(I) 
        END IF
        NNG(I)=NM2*NM3*NNG1+NM3*NNG2+NNG3+1 
  40  CONTINUE
C
c reorder of reciprocal vector (not neccesory except for G=0)
      DO 50 I=1,NUMG 
      DO 50 J=I,NUMG
        DKGG1=NG1(I)*G1LVX+NG2(I)*G2LVX+NG3(I)*G3LVX
        DKGG2=NG1(I)*G1LVY+NG2(I)*G2LVY+NG3(I)*G3LVY
        DKGG3=NG1(I)*G1LVZ+NG2(I)*G2LVZ+NG3(I)*G3LVZ
        DKG2I=DKGG1*DKGG1+DKGG2*DKGG2+DKGG3*DKGG3
        DKGG1=NG1(J)*G1LVX+NG2(J)*G2LVX+NG3(J)*G3LVX
        DKGG2=NG1(J)*G1LVY+NG2(J)*G2LVY+NG3(J)*G3LVY
        DKGG3=NG1(J)*G1LVZ+NG2(J)*G2LVZ+NG3(J)*G3LVZ
        DKG2J=DKGG1*DKGG1+DKGG2*DKGG2+DKGG3*DKGG3
        IF(DKG2I.GT.DKG2J) THEN
          NTEMP1=NG1(I)
          NTEMP2=NG2(I)
          NTEMP3=NG3(I)
          NNTEMP=NNG(I)
          NG1(I)=NG1(J)
          NG2(I)=NG2(J)
          NG3(I)=NG3(J)
          NNG(I)=NNG(J)
          NG1(J)=NTEMP1
          NG2(J)=NTEMP2
          NG3(J)=NTEMP3
          NNG(J)=NNTEMP
        ELSE IF(DKG2I.EQ.DKG2J) THEN
          IF(NG1(I).GT.NG1(J)) THEN
            NTEMP1=NG1(I)
            NTEMP2=NG2(I)
            NTEMP3=NG3(I)
            NNTEMP=NNG(I)
            NG1(I)=NG1(J)
            NG2(I)=NG2(J)
            NG3(I)=NG3(J)
            NNG(I)=NNG(J)
            NG1(J)=NTEMP1
            NG2(J)=NTEMP2
            NG3(J)=NTEMP3
            NNG(J)=NNTEMP
          ELSE IF(NG1(I).EQ.NG1(J)) THEN
             IF(NG2(I).GT.NG2(J)) THEN
               NTEMP1=NG1(I)
               NTEMP2=NG2(I)
               NTEMP3=NG3(I)
               NNTEMP=NNG(I)
               NG1(I)=NG1(J)
               NG2(I)=NG2(J)
               NG3(I)=NG3(J)
               NNG(I)=NNG(J)
               NG1(J)=NTEMP1
               NG2(J)=NTEMP2
               NG3(J)=NTEMP3
               NNG(J)=NNTEMP
             ELSE IF(NG2(I).EQ.NG2(J)) THEN
               IF(NG3(I).GT.NG3(J)) THEN
                 NTEMP1=NG1(I)
                 NTEMP2=NG2(I)
                 NTEMP3=NG3(I)
                 NNTEMP=NNG(I)
                 NG1(I)=NG1(J)
                 NG2(I)=NG2(J)
                 NG3(I)=NG3(J)
                 NNG(I)=NNG(J)
                 NG1(J)=NTEMP1
                 NG2(J)=NTEMP2
                 NG3(J)=NTEMP3
                 NNG(J)=NNTEMP
               END IF
             END IF
          END IF
        END IF
   50 CONTINUE
C 
c storing reciprocal vector
      DO 60 I=1,NUMG
        GG1(I)=DBLE(NG1(I))*G1LVX+DBLE(NG2(I))*G2LVX
     &           +DBLE(NG3(I))*G3LVX
        GG2(I)=DBLE(NG1(I))*G1LVY+DBLE(NG2(I))*G2LVY
     &           +DBLE(NG3(I))*G3LVY
        GG3(I)=DBLE(NG1(I))*G1LVZ+DBLE(NG2(I))*G2LVZ
     &           +DBLE(NG3(I))*G3LVZ
   60 CONTINUE
C
      NUML=1
C
      RETURN
      END
C
C **** PREPME   ********************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE PREPME
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (NM=2300)
      PARAMETER (NA=3,NAN=(NA*(NA+1))/2)
      PARAMETER (LLM3=NM*NM)
      PARAMETER (LM3=NM*NM)
      PARAMETER (PI=3.1415926535897932384626D0)
      PARAMETER (NGX=64,NGY=64,NGZ=64)
      PARAMETER (NXYZ=NGX*NGY*NGZ)
      PARAMETER (IGG3=100000)
C
      COMMON /NUMBE/ NATOM
      COMMON /ATOM / WM(NA),NKA
      COMMON /IDENT/ KIND(NM)
      COMMON /NUMGL/ NUML,NUMG
      COMMON /GGGGG/ GG1(IGG3),GG2(IGG3),GG3(IGG3)
      COMMON /BDBDB/ BDASHR(IGG3),BDASHI(IGG3)
      COMMON /RECIP/ OMC
C
C *** EWALD
C
      COMMON /GAMMA/ GAM
      COMMON /NUMVV/ DNZ(NA),DNUMV(NAN)
      COMMON /PRECA/ EWASUM3,EWASUM22
      COMMON /PMEPM/ IBSOD,NM1,NM2,NM3,NK1,NK2,NK3
      COMMON /FFTPM/ NG1(IGG3),NG2(IGG3),NG3(IGG3),NNG(IGG3)
      COMMON /PMMAT/ CC(IGG3),QMATR(NXYZ),QMATI(NXYZ)
C
C *** 1 EV = 1.60219E-19 J
C     DATA ELV,ELV2/1.60219D-19,2.567012796D-38/
      DATA ELV,ELV2/1.60219D-19,2.3071132D-28/
      DATA ELV2A /2.3071132D-18/
C
      R2RTPI=2.0D0/DSQRT(PI)
      RTPI=1.D0/DSQRT(PI)
      RGAM=1.0D0/(GAM*GAM)
c
c constant term of ewald method(never change in simulation)
      EWASUM3=0.0D0
      DO IA=1,NATOM
        DNUMVA=DNZ(KIND(IA))
        EWASUM3=EWASUM3+DNUMVA*DNUMVA*GAM*RTPI*ELV2
      ENDDO
C
c
c B(G)*C(G) ref.1 eq(3.9) and eq(4.8)
      DO I=1,NUMG
        CC(I)=0D0
      ENDDO
      DO IGGG=1,NUMG
c
c ref.1 eq(4.4) B(G)
       DRSUM1=0D0
       DISUM1=0D0
       DRSUM2=0D0
       DISUM2=0D0
       DRSUM3=0D0
       DISUM3=0D0
       DO K=0,IBSOD-2
         CALL MMN(IBSOD,DBLE(K+1),DMNK1)
         DRSUM1=DRSUM1+DMNK1*DCOS(2D0*PI*DBLE(NG1(IGGG)*K)/DBLE(NM1))
         DISUM1=DISUM1+DMNK1*DSIN(2D0*PI*DBLE(NG1(IGGG)*K)/DBLE(NM1))
         DRSUM2=DRSUM2+DMNK1*DCOS(2D0*PI*DBLE(NG2(IGGG)*K)/DBLE(NM2))
         DISUM2=DISUM2+DMNK1*DSIN(2D0*PI*DBLE(NG2(IGGG)*K)/DBLE(NM2))
         DRSUM3=DRSUM3+DMNK1*DCOS(2D0*PI*DBLE(NG3(IGGG)*K)/DBLE(NM3))
         DISUM3=DISUM3+DMNK1*DSIN(2D0*PI*DBLE(NG3(IGGG)*K)/DBLE(NM3))
       ENDDO
       IF(IABS(NG1(IGGG))*2.EQ.NM1) THEN
         DRSUM1=0d0
         DISUM1=0d0
       ENDIF
       IF(IABS(NG2(IGGG))*2.EQ.NM2) THEN
         DRSUM2=0d0
         DISUM2=0d0
       ENDIF
       IF(IABS(NG3(IGGG))*2.EQ.NM3) THEN
         DRSUM3=0d0
         DISUM3=0d0
       ENDIF
       DRR1=DCOS(2D0*PI*DBLE((IBSOD-1)*NG1(IGGG))/DBLE(NM1))
       DII1=DSIN(2D0*PI*DBLE((IBSOD-1)*NG1(IGGG))/DBLE(NM1))
       DRR2=DCOS(2D0*PI*DBLE((IBSOD-1)*NG2(IGGG))/DBLE(NM2))
       DII2=DSIN(2D0*PI*DBLE((IBSOD-1)*NG2(IGGG))/DBLE(NM2))
       DRR3=DCOS(2D0*PI*DBLE((IBSOD-1)*NG3(IGGG))/DBLE(NM3))
       DII3=DSIN(2D0*PI*DBLE((IBSOD-1)*NG3(IGGG))/DBLE(NM3))
       BBM1=DRSUM1*DRSUM1+DISUM1*DISUM1
       BBM2=DRSUM2*DRSUM2+DISUM2*DISUM2
       BBM3=DRSUM3*DRSUM3+DISUM3*DISUM3
ccccccccccccccc
       SSR=DRR1*DRSUM1+DII1*DISUM1
       SSI=DII1*DRSUM1-DRR1*DISUM1
       PPR=DRR2*DRSUM2+DII2*DISUM2
       PPI=DII2*DRSUM2-DRR2*DISUM2
       QQR=DRR3*DRSUM3+DII3*DISUM3
       QQI=DII3*DRSUM3-DRR3*DISUM3
       BR=SSR*PPR*QQR-SSI*PPI*QQR-SSR*PPI*QQI-SSI*PPR*QQI
       BI=SSR*PPI*QQR+SSI*PPR*QQR-SSR*PPR*QQI-SSI*PPI*QQI 
       IF(DABS(BBM1*BBM2*BBM3).GT.1D-15) THEN
         BR=BR/(BBM1*BBM2*BBM3)
         BI=BI/(BBM1*BBM2*BBM3)
       ELSE
         BR=0D0
         BI=0D0
       ENDIF
       BDASHR(IGGG)=BR
       BDASHI(IGGG)=BI
      ENDDO
c
c C(G)
      DO IGGG=2,NUMG
       G2=GG1(IGGG)*GG1(IGGG)+GG2(IGGG)*GG2(IGGG)+GG3(IGGG)*GG3(IGGG)
       IF(0.25D0*RGAM*G2.GT.40D0) THEN
         RG2G2=0.D0
       ELSE 
         RG2G2=4D0*PI*OMC*DEXP(-0.25D0*RGAM*G2)/G2
       END IF 
       CC(IGGG)=RG2G2
      ENDDO
C
      RETURN
      END
C **** PME    **********************************************************
C *                                                                    *
C *      THIS IS MAIN ROUTINE FOR PME                                  *
C *                                                                    *
C **********************************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE PME(NATOM,QX,QY,QZ,RLV,RLVI)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (NM=2300)
      PARAMETER (NA=3,NAN=(NA*(NA+1))/2)
      PARAMETER (LLM3=NM*NM)
      PARAMETER (LM3=NM*NM)
      PARAMETER (PI=3.1415926535897932384626D0)
      PARAMETER (NGX=64,NGY=64,NGZ=64)
      PARAMETER (NXYZ=NGX*NGY*NGZ)
      PARAMETER (IGG3=100000)
C
      COMMON /IDENT/ KIND(NM)
      COMMON /ATOM / WM(NA),NKA
      COMMON /LVEC3/ LL3I(LM3),LL3J(LM3),KKMAX3,KMAX3
      COMMON /PERBC/ IPBCX,IPBCY,IPBCZ
      COMMON /RECV1/ G1LVX,G1LVY,G1LVZ
      COMMON /RECV2/ G2LVX,G2LVY,G2LVZ
      COMMON /RECV3/ G3LVX,G3LVY,G3LVZ
      COMMON /RECIP/ OMC
      COMMON /NUMGL/ NUML,NUMG
      COMMON /GGGGG/ GG1(IGG3),GG2(IGG3),GG3(IGG3)
      COMMON /BDBDB/ BDASHR(IGG3),BDASHI(IGG3)
C
      COMMON /GAMMA/ GAM
      COMMON /NUMVV/ DNZ(NA),DNUMV(NAN)
      COMMON /PRECA/ EWASUM3,EWASUM22
      COMMON /PMEPM/ IBSOD,NM1,NM2,NM3,NK1,NK2,NK3
      COMMON /FFTPM/ NG1(IGG3),NG2(IGG3),NG3(IGG3),NNG(IGG3)
      COMMON /PMMAT/ CC(IGG3),QMATR(NXYZ),QMATI(NXYZ)
C
      DIMENSION QX(NM),QY(NM),QZ(NM)
      DIMENSION DDQX(NM),DDQY(NM),DDQZ(NM)
C
      COMMON /PMEP2/ RCUTPML,RCUTPMS
C
      DIMENSION LKINT3(LLM3),Q3(LLM3)
      DIMENSION QIJ3X(LLM3),QIJ3Y(LLM3),QIJ3Z(LLM3)
      DIMENSION RLV(3,3),RLVD(3,3),RLVDD(3,3),RLVI(3,3)
      DIMENSION L3I(LLM3),L3J(LLM3)
C
C *** 1 EV = 1.60219E-19 J
      DATA ELV,ELV2/1.60219D-19,2.3071132D-28/
      DATA ELV2A /2.3071132D-18/
C
      R2RTPI=2.0D0/DSQRT(PI)
      RTPI=1.D0/DSQRT(PI)
      RGAM=1.0D0/(GAM*GAM)
C
C *** SET CHILD LIST VECTOR
C
      K=0
      DO 20 KK=1,KKMAX3
        I=LL3I(KK)
        J=LL3J(KK)
        QXIJ=QX(I)-QX(J)
        QYIJ=QY(I)-QY(J)
        QZIJ=QZ(I)-QZ(J)
C
        IF(IPBCX.EQ.0) THEN
          IF(QXIJ.LT.-0.5) QXIJ=QXIJ+1D0
          IF(QXIJ.GE.0.5) QXIJ=QXIJ-1D0
        ENDIF
C
        IF(IPBCY.EQ.0) THEN
          IF(QYIJ.LT.-0.5) QYIJ=QYIJ+1D0
          IF(QYIJ.GE.0.5) QYIJ=QYIJ-1D0
        ENDIF
C
        IF(IPBCZ.EQ.0) THEN
          IF(QZIJ.LT.-0.5) QZIJ=QZIJ+1D0
          IF(QZIJ.GE.0.5) QZIJ=QZIJ-1D0
        ENDIF
C
        QXIJT=QXIJ*RLV(1,1)+QYIJ*RLV(1,2)+QZIJ*RLV(1,3)
        QYIJT=QXIJ*RLV(2,1)+QYIJ*RLV(2,2)+QZIJ*RLV(2,3)
        QZIJT=QXIJ*RLV(3,1)+QYIJ*RLV(3,2)+QZIJ*RLV(3,3)
        QXIJ=QXIJT
        QYIJ=QYIJT
        QZIJ=QZIJT
C
        QQQ=QXIJ*QXIJ+QYIJ*QYIJ+QZIJ*QZIJ
        KINT=((NKA+NKA-KIND(I))*(KIND(I)-1))/2+KIND(J)
C
C **** SET POTENTIAL CUT OFF
C
        IF(QQQ.LT.RCUTPMS**2) THEN
          K=K+1
          QIJ3X(K)=QXIJ
          QIJ3Y(K)=QYIJ
          QIJ3Z(K)=QZIJ
          Q3(K)=DSQRT(QQQ)
          L3I(K)=I
          L3J(K)=J
          LKINT3(K)=KINT
        END IF
   20 CONTINUE
      KMAX3=K
      write(*,*)'kmax3=',kmax3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C ***** EWALD'S ELECTROSTATIC ENERGY
C
        EP=0D0
        DO I=1,NM1*NM2*NM3
          QMATR(I)=0D0
          QMATI(I)=0D0
        ENDDO
        DO N=1,NATOM
c
c charge of atom N
          QI=DNZ(KIND(N))
c
c scaled fractional coordinate
          UN1=QX(N)*DBLE(NM1)
          UN2=QY(N)*DBLE(NM2)
          UN3=QZ(N)*DBLE(NM3)
c
c calculate ref.1 eq(4.6)
          DO IKX=-NM1/2,NM1/2-1
          DO IKY=-NM2/2,NM2/2-1
          DO IKZ=-NM3/2,NM3/2-1
            IF(IKX.LT.0) THEN
              IK1=IKX+NM1 
            ELSE
              IK1=IKX 
            END IF
            IF(IKY.LT.0) THEN
              IK2=IKY+NM2 
            ELSE
              IK2=IKY 
            END IF
            IF(IKZ.LT.0) THEN
              IK3=IKZ+NM3 
            ELSE
              IK3=IKZ 
            END IF
            KKK=NM2*NM3*IK1+NM3*IK2+IK3+1
            DMN1=0D0
            DO N1=-NK1,NK1
              U=UN1-DBLE(IK1)-DBLE(N1*NM1)
              CALL MMN(IBSOD,U,VAL)
              DMN1=DMN1+VAL
            ENDDO
            DMN2=0D0
            DO N2=-NK2,NK2
              U=UN2-DBLE(IK2)-DBLE(N2*NM2)
              CALL MMN(IBSOD,U,VAL)
              DMN2=DMN2+VAL
            ENDDO
            DMN3=0D0
            DO N3=-NK3,NK3
              U=UN3-DBLE(IK3)-DBLE(N3*NM3)
              CALL MMN(IBSOD,U,VAL)
              DMN3=DMN3+VAL
            ENDDO
            QMATR(KKK)=QMATR(KKK)+QI*DMN1*DMN2*DMN3
          ENDDO
          ENDDO
          ENDDO
        ENDDO
c
c  Q(s) --> Q(-G)
        CALL FFT3D(1D0,NM1,NM2,NM3,QMATR,QMATI)
c
c  V=sum_G B(G)C(G)[Q(G)]**2
        EREC=0D0
        SXX=0D0
        SYY=0D0
        SZZ=0D0
        SXY=0D0
        SYZ=0D0
        SZX=0D0
        DO IGGG=2,NUMG
          GGG=GG1(IGGG)*GG1(IGGG)
     &       +GG2(IGGG)*GG2(IGGG)
     &       +GG3(IGGG)*GG3(IGGG)
          CO=2d0/GGG+0.5d0*RGAM
          CXX=1d0-GG1(IGGG)*GG1(IGGG)*CO
          CYY=1d0-GG2(IGGG)*GG2(IGGG)*CO
          CZZ=1d0-GG3(IGGG)*GG3(IGGG)*CO
          CXY=GG1(IGGG)*GG2(IGGG)*CO
          CYZ=GG2(IGGG)*GG3(IGGG)*CO
          CZX=GG3(IGGG)*GG1(IGGG)*CO
          BB2=BDASHR(IGGG)*BDASHR(IGGG)+BDASHI(IGGG)*BDASHI(IGGG)
          ERECC=BB2*CC(IGGG)
     &          *(QMATR(NNG(IGGG))**2+QMATI(NNG(IGGG))**2)
c energy
          EREC=EREC+ERECC
c stress
          SXX=SXX+CXX*EREC
          SYY=SYY+CYY*EREC
          SZZ=SZZ+CZZ*EREC
          SXY=SXY+CXY*EREC
          SYZ=SYZ+CYZ*EREC
          SZX=SZX+CZX*EREC
c
        ENDDO
        EP=EREC*0.5D0*ELV2
        SXX=SXX*0.5d0*ELV2*OMC
        SYY=SYY*0.5d0*ELV2*OMC
        SZZ=SZZ*0.5d0*ELV2*OMC
        SXY=SXY*0.5d0*ELV2*OMC
        SYZ=SYZ*0.5d0*ELV2*OMC
        SZX=SZX*0.5d0*ELV2*OMC
        add1=ep
c
c local term of Ewald method
         add21=0d0
        DO 53 K=1,KKMAX3
          I=L3I(K)
          J=L3J(K)
          KI=LKINT3(K)
          DNUMVA=DNUMV(KI)
          DRT=Q3(K)
          DTAX=QIJ3X(K)
          DTAY=QIJ3Y(K)
          DTAZ=QIJ3Z(K)
          IF(DRT.LT.1d-20) GOTO 52
c energy
          ADD=DNUMVA*(1.0D0-DERF(DRT*GAM))/DRT*ELV2
          EP=EP+ADD
          add21=add21+add
c stress
          CO=(1.0D0-DERF(DRT*GAM))/(DRT*DRT*DRT)
     &      +(2d0*GAM*RTPI*DEXP(-GAM*GAM*DRT*DRT))/DRT*DRT
          SXX=SXX+CO*DTAX*DTAX*DNUMVA*ELV2
          SYY=SYY+CO*DTAY*DTAY*DNUMVA*ELV2
          SZZ=SZZ+CO*DTAZ*DTAZ*DNUMVA*ELV2
          SXY=SXY+CO*DTAX*DTAY*DNUMVA*ELV2
          SYZ=SYZ+CO*DTAY*DTAZ*DNUMVA*ELV2
          SZX=SZX+CO*DTAZ*DTAX*DNUMVA*ELV2
   52   CONTINUE
   53   CONTINUE
c
        EP=EP-EWASUM3
c
c EP:total electrostatic energy!!
c S??: stress tensor
c
        write(*,*)gam,ep,add1,add21,ewasum3
        write(*,*)'stress',sxx,syy,szz,sxy,syz,szx
C
C ***** EWALD'S ELECTROSTATIC FORCE
C
c
c ddq?(i) force vector acting on atom i
        DO I=1,NATOM
          DDQX(I)=0D0
          DDQY(I)=0D0
          DDQZ(I)=0D0
        ENDDO
c
        DO I=1,NATOM
          QI=DNZ(KIND(I))
          QXI=QX(I)*RLV(1,1)+QY(I)*RLV(1,2)+QZ(I)*RLV(1,3)
          QYI=QX(I)*RLV(2,1)+QY(I)*RLV(2,2)+QZ(I)*RLV(2,3)
          QZI=QX(I)*RLV(3,1)+QY(I)*RLV(3,2)+QZ(I)*RLV(3,3)
          FXI=0.0D0
          FYI=0.0D0
          FZI=0.0D0
          DO IGGG=2,NUMG
            GT=GG1(IGGG)*QXI+GG2(IGGG)*QYI+GG3(IGGG)*QZI
            DSI=DSIN(GT)
            DCO=DCOS(GT)
            CBQR=CC(IGGG)*(BDASHR(IGGG)*QMATR(NNG(IGGG))
     &                    +BDASHI(IGGG)*QMATI(NNG(IGGG)))
            CBQI=CC(IGGG)*(BDASHR(IGGG)*QMATI(NNG(IGGG))
     &                    -BDASHI(IGGG)*QMATR(NNG(IGGG)))
            PROD=CBQR*DSI+CBQI*DCO
            FXIT=GG1(IGGG)*QI*PROD*ELV2
            FYIT=GG2(IGGG)*QI*PROD*ELV2
            FZIT=GG3(IGGG)*QI*PROD*ELV2
            FXI=FXI+FXIT
            FYI=FYI+FYIT
            FZI=FZI+FZIT
          ENDDO
          DDQX(I)=DDQX(I)+FXI
          DDQY(I)=DDQY(I)+FYI
          DDQZ(I)=DDQZ(I)+FZI
        ENDDO
c
c force from local term
        DO 33 K=1,KKMAX3
          I=LL3I(K)
          J=LL3J(K)
          KI=LKINT3(K)
          DNUMVA=DNUMV(KI)
          DTAX=QIJ3X(K)
          DTAY=QIJ3Y(K)
          DTAZ=QIJ3Z(K)
          DRTX=DTAX
          DRTY=DTAY
          DRTZ=DTAZ
          DRT2=DRTX*DRTX+DRTY*DRTY+DRTZ*DRTZ
          DRT=DSQRT(DRT2)
          if(drt.lt.1d-20) goto 32
          DRT3=1.0D0/(DRT*DRT2*1D30)
          IF(DRT2*GAM*GAM.GT.40D0) THEN
            DER=0.D0  
          ELSE
            DER=-R2RTPI*DEXP(-DRT2*GAM*GAM)
          END IF
          ER=1.0D0-DERF(DRT*GAM)
          GDE=ER-GAM*DER*DRT
          PROD=DNUMVA*GDE*DRT3*ELV2*1D30
          FXIJ=DRTX*PROD
          FYIJ=DRTY*PROD
          FZIJ=DRTZ*PROD
          DDQX(I)=DDQX(I)+FXIJ
          DDQY(I)=DDQY(I)+FYIJ
          DDQZ(I)=DDQZ(I)+FZIJ
          DDQX(J)=DDQX(J)-FXIJ
          DDQY(J)=DDQY(J)-FYIJ
          DDQZ(J)=DDQZ(J)-FZIJ
   32   CONTINUE
   33   CONTINUE
        do i=1,natom
        write(150,*)i,ddqx(i),ddqy(i),ddqz(i)
        enddo
c
      RETURN
      END
C **** INPUTG **********************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE INPUTG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (PI=3.1415926535897932384626D0)
      PARAMETER (EPS=1.0D-23)
      PARAMETER (NA=3,NAN=(NA*(NA+1))/2)
      PARAMETER (NM=2300)
      PARAMETER (NB=10)
C
      COMMON /ATOM / WM(NA),NKA
      COMMON /BOUFN/ NBOUNF,NF(NB),MBX(NB),MBY(NB),MBZ(NB)
      COMMON /BOUFV/ BFX(NB),BFY(NB),BFZ(NB)
      COMMON /BOUQN/ NBOUNQ,NQ(NB),NBX(NB),NBY(NB),NBZ(NB)
      COMMON /BOUQV/ BQX(NB),BQY(NB),BQZ(NB)
      COMMON /DDIF1/ DDQX1(NM),DDQY1(NM),DDQZ1(NM),DDS1,RLVDD1(3,3)
      COMMON /DDIF2/ DDQX2(NM),DDQY2(NM),DDQZ2(NM),DDS2,RLVDD2(3,3)
      COMMON /DDIF3/ DDQX3(NM),DDQY3(NM),DDQZ3(NM),DDS3,RLVDD3(3,3)
      COMMON /ENERG/ EP,EK,EZ,ES,ET,EW,EWW,EH
      COMMON /HEATB/ TONDO,QW,IPROB
      COMMON /STCON/ SIGXX,SIGYY,SIGZZ,SIGXY,SIGYZ,SIGZX,WW
      COMMON /IDENT/ KIND(NM)
c      COMMON /INITQ/ QX0(NM),QY0(NM),QZ0(NM)
      COMMON /INITE/ S0,Z0
      COMMON /NORML/ RNOR,ENOR,TNOR
      COMMON /NUMBE/ NATOM
      COMMON /PERBC/ IPBCX,IPBCY,IPBCZ
      COMMON /STEP / NEOUT,NGOUT,NSTEP
      COMMON /TIMEC/ T,TSTART,TFINIS,DT0,DT1,DT2
      COMMON /VALUQ/ QX(NM),QY(NM),QZ(NM),DQX(NM),DQY(NM),DQZ(NM)
      COMMON /UNIV1/ RLV(3,3),RLVI(3,3),RLVD(3,3)
      COMMON /RECV1/ G1LVX,G1LVY,G1LVZ
      COMMON /RECV2/ G2LVX,G2LVY,G2LVZ
      COMMON /RECV3/ G3LVX,G3LVY,G3LVZ
      COMMON /RECIP/ OMC
      COMMON /HEATA/ TONDODEC
      COMMON /IPRO2/ IRORM,IMEAN
      COMMON /EXTRF/ FEXX,FEXY,FEXZ
C
C *** INPUT UNIT CELL PARAMETER
C
      READ(3,*) RLV(1,1),RLV(2,1),RLV(3,1)
      READ(3,*) RLV(1,2),RLV(2,2),RLV(3,2)
      READ(3,*) RLV(1,3),RLV(2,3),RLV(3,3)
      DO I=1,3
        DO J=1,3
          RLVI(I,J)=RLV(I,J)
        ENDDO
      ENDDO
C
C *** CALCULATE INVERSE MATRIX
C
      CALL INVERS(3,RLVI)
C
C *** CALCULATE RECIPROCAL VECTOR OF UNIT CELL
C
      OMC=1.0D0/(RLV(1,1)*(RLV(2,2)*RLV(3,3)-RLV(3,2)*RLV(2,3))
     &          +RLV(2,1)*(RLV(3,2)*RLV(1,3)-RLV(1,2)*RLV(3,3))
     &          +RLV(3,1)*(RLV(1,2)*RLV(2,3)-RLV(2,2)*RLV(1,3)))
      OMC0=OMC
      G1LVX=2.0D0*PI*OMC*(RLV(2,2)*RLV(3,3)-RLV(3,2)*RLV(2,3))
      G1LVY=2.0D0*PI*OMC*(RLV(3,2)*RLV(1,3)-RLV(1,2)*RLV(3,3))
      G1LVZ=2.0D0*PI*OMC*(RLV(1,2)*RLV(2,3)-RLV(2,2)*RLV(1,3))
      G2LVX=2.0D0*PI*OMC*(RLV(2,3)*RLV(3,1)-RLV(3,3)*RLV(2,1))
      G2LVY=2.0D0*PI*OMC*(RLV(3,3)*RLV(1,1)-RLV(1,3)*RLV(3,1))
      G2LVZ=2.0D0*PI*OMC*(RLV(1,3)*RLV(2,1)-RLV(2,3)*RLV(1,1))
      G3LVX=2.0D0*PI*OMC*(RLV(2,1)*RLV(3,2)-RLV(3,1)*RLV(2,2))
      G3LVY=2.0D0*PI*OMC*(RLV(3,1)*RLV(1,2)-RLV(1,1)*RLV(3,2))
      G3LVZ=2.0D0*PI*OMC*(RLV(1,1)*RLV(2,2)-RLV(2,1)*RLV(1,2))
C
C *** NUMBER OF ATOM TYPE
C
      READ(3,*) NKA
      NKAN=(NKA*(NKA+1))/2
C
C *** THE NUMBER OF ATOMS
C
      READ(3,*) NATOM
C
C *** COORDINATIES & MOMENTUM
C
      DO 50 I=1,NATOM
        READ(1,*) ITEMP,KIND(I),QX(I),QY(I),QZ(I)
        READ(1,*) DQX(I),DQY(I),DQZ(I)
   50 CONTINUE
      RETURN
      END
C **** LSVEC ***********************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE LSVEC
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (NM=2300)
      PARAMETER (LM3=NM*NM)
C
      COMMON /LVEC3/ LL3I(LM3),LL3J(LM3),KKMAX3,KMAX3
      COMMON /NUMBE/ NATOM
      COMMON /IDENT/ KIND(NM)
      COMMON /VALUQ/ QX(NM),QY(NM),QZ(NM),DQX(NM),DQY(NM),DQZ(NM)
      COMMON /UNIV1/ RLV(3,3),RLVI(3,3),RLVD(3,3)
      COMMON /PERBC/ IPBCX,IPBCY,IPBCZ
C
      COMMON /PMEP2/ RCUTPML,RCUTPMS
C
      KK3=0
      DO 10 I=1,NATOM
      DO 10 J=1,NATOM
C
        QXIJ=QX(I)-QX(J)
        QYIJ=QY(I)-QY(J)
        QZIJ=QZ(I)-QZ(J)
C
        IF(IPBCX.EQ.0) THEN
          IF(QXIJ.LT.-0.5) QXIJ=QXIJ+1D0
          IF(QXIJ.GE.0.5) QXIJ=QXIJ-1D0
        ENDIF
C
        IF(IPBCY.EQ.0) THEN
          IF(QYIJ.LT.-0.5) QYIJ=QYIJ+1D0
          IF(QYIJ.GE.0.5) QYIJ=QYIJ-1D0
        ENDIF
C
        IF(IPBCZ.EQ.0) THEN
          IF(QZIJ.LT.-0.5) QZIJ=QZIJ+1D0
          IF(QZIJ.GE.0.5) QZIJ=QZIJ-1D0
        ENDIF
C
        QXIJT=QXIJ*RLV(1,1)+QYIJ*RLV(1,2)+QZIJ*RLV(1,3)
        QYIJT=QXIJ*RLV(2,1)+QYIJ*RLV(2,2)+QZIJ*RLV(2,3)
        QZIJT=QXIJ*RLV(3,1)+QYIJ*RLV(3,2)+QZIJ*RLV(3,3)
        QXIJ=QXIJT
        QYIJ=QYIJT
        QZIJ=QZIJT
C
        QQ=QXIJ*QXIJ+QYIJ*QYIJ+QZIJ*QZIJ
C
        IF(QQ.LT.RCUTPML**2) THEN
          IF(J.GT.I) THEN
            KK3=KK3+1
            LL3I(KK3)=I
            LL3J(KK3)=J
          ENDIF
        ENDIF
   10 CONTINUE
C
      IF(KK3.GT.LM3) THEN
        WRITE(*,*) 'KK3,LM3=',KK3,LM3
        STOP 'TOO LARGE SYSTEM'
      ENDIF
      KKMAX3=KK3
      WRITE(*,*) 'kkmax,lm=',KKMAX3,LM3
C
      RETURN
      END
C **** INVERS **********************************************************
C **********************************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE INVERS(N,AA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AA(3,3),A(3,6),B(3,6)
C
      DO 1 I=1,N
        DO 1 J=1,N
  1     A(I,J)=AA(I,J)
      NP1=N+1
      M=N+N
C
C     GENERATE AUGMENTING IDENTITY MATRIX
C
      DO 8 J=NP1,M
        DO 8 I=1,N
        IF(J-N-I) 7,6,7
  6     A(I,J)=1.0
        GOTO 8
  7     A(I,J)=0.0
  8   CONTINUE
C
  9   IF(A(1,1)) 15,10,15
 10   K=M-N
      DO 13 I=2,K
        IF(A(I,1)) 11,13,11
 11     DO 12 J=1,M
          TEMP=A(I,J)
          A(I,J)=A(1,J)
          A(1,J)=TEMP
 12   CONTINUE
      GOTO 15
 13   CONTINUE
      WRITE(*,*)'MATRIX IS SINGULAR'
      STOP
C
 15   DO 16 J=2,M
      DO 16 I=2,N
 16     B(I-1,J-1)=A(I,J)-A(1,J)*A(I,1)/A(1,1)
      DO 17 J=2,M
 17     B(N,J-1)=A(1,J)/A(1,1)
      M=M-1
      DO 18 J=1,M
        DO 18 I=1,N
 18     A(I,J)=B(I,J)
      IF(M-N) 9,19,9
 19   CONTINUE
      DO 20 I=1,N
        DO 20 J=1,N
        AA(I,J)=A(I,J)
 20   CONTINUE
      RETURN
      END
C
C **** DDERF ***********************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      FUNCTION DDERF(ARG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IF(ARG.LT.0.0D0) THEN
        DDERF=0.0D0
      ELSE
        T=1.0D0/(1.0D0+0.3275911D0*ARG)
        TD=-0.3275911D0/(1.0D0+0.3275911D0*ARG)
     &                 /(1.0D0+0.3275911D0*ARG)
        T2=T*T
        T3=T2*T
        POLY=T*0.254829592D0+T2*(-0.284496736D0)+T3*1.421413741D0
     &         +T2*T2*(-1.453152027D0)+T2*T3*1.061405429D0
        POLYD=(0.254829592D0
     &           +2.0D0*T*(-0.284496736D0)
     &           +3.0D0*T2*1.421413741D0
     &           +4.0D0*T*T2*(-1.453152027D0)
     &           +5.0D0*T2*T2*1.061405429D0)*TD
        DDERF=(2.0D0*POLY*ARG-POLYD)*DEXP(-ARG*ARG)
      ENDIF
      RETURN
      END
C **** DERF *************************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      FUNCTION DERF(ARG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IF(ARG.LT.0.0D0) THEN
        DERF=0.0D0
      ELSE
        T=1.0D0/(1.0D0+0.3275911D0*ARG)
        T2=T*T
        T3=T2*T
        POLY=T*0.254829592D0+T2*(-0.284496736D0)+T3*1.421413741D0
     &         +T2*T2*(-1.453152027D0)+T2*T3*1.061405429D0
        DERF=1.0D0-POLY*DEXP(-ARG*ARG)
      ENDIF
      RETURN
      END
C **********************************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      BLOCK DATA INITIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NA=3,NAN=(NA*(NA+1))/2)
C
      COMMON /NUMVV/ DNZ(NA),DNUMV(NAN)
      COMMON /GAMMA/ GAM
      COMMON /CUTEW/ CUTEG,CUTEL
C
      DATA DNZ/-1.104,1.472,0E+00/
      DATA DNUMV/1.218816,-1.625088,2.166784,0D0,0D0,0D0/
      DATA GAM/0.30D+10/
      DATA CUTEG,CUTEL/5e10,3/
C
      END
C **********************************************************************
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      SUBROUTINE MMN(N,U,VAL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ANS(10,10)
      IF(U.GE.0D0.AND.U.LE.DBLE(N))THEN
      DO IU=1,N-1
        IF(U.GE.DBLE(IU-1).AND.U.LE.DBLE(IU-1)+2D0) THEN
          ANS(IU,2)=1D0-DABS(U-DBLE(IU))
        ELSE
          ANS(IU,2)=0D0
        ENDIF
      ENDDO
      DO I=3,N
        DO J=1,N-I+1
        ANS(J,I)=((U-DBLE(J-1))/DBLE(I-1))*ANS(J,I-1)
     &   +((DBLE(I)-U+DBLE(J-1))/DBLE(I-1))*ANS(J+1,I-1)
        ENDDO
      ENDDO
      VAL=ANS(1,N)
      ELSE
      VAL=0D0
      ENDIF
      RETURN
      END
C***********************************************************
C  3-Dimensional FFT(Fast-Fourier-Transform)  program
C
C  TABLE.F RFFT.F CFFT.F VFFT.F FFT.F
C
      SUBROUTINE FFT3D(QSS,NR,NC,NV,XX,YY)
C***********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IMS3=262144)
      COMMON /IFLG1/ IRF,IFF
      COMMON /IFLG2/ QS
      COMMON /IDAT4/ IGC,IGR,IGV
      COMMON /IDAT1/ X(IMS3),Y(IMS3) 
      COMMON /IDAT5/ NMAX
      DIMENSION XX(IMS3),YY(IMS3)
C
      QS=QSS
C
      DO 10 I=1,NR*NC*NV
        X(I)=XX(I)
        Y(I)=YY(I)
   10 CONTINUE 
C
      CALL RFFT(NR,NC,NV,NMAX)
      CALL CFFT(NR,NC,NV,NMAX)
      CALL VFFT(NR,NC,NV,NMAX)
C
      DO 20 I=1,NR*NC*NV
        XX(I)=X(I)
        YY(I)=Y(I)
   20 CONTINUE 
C    
      RETURN
      END
C
C***********************************************************
      SUBROUTINE TABLE(NR,NC,NV)
C***********************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (PI=3.1415926535897932384626D0)
      PARAMETER (IMS=64)
      PARAMETER (IMS3=262144)
      COMMON /IFLG1/ IRF,IFF
      COMMON /IFLG2/ QS
      COMMON /IDAT3/ C(IMS/2),S(IMS/2),IK(IMS)
      COMMON /IDAT4/ IGC,IGR,IGV
      COMMON /IDAT5/ NMAX
C
      IRF=1
      IGR=IDNINT(DLOG(DBLE(NR))/DLOG(2D0))
      IGC=IDNINT(DLOG(DBLE(NC))/DLOG(2D0))
      IGV=IDNINT(DLOG(DBLE(NV))/DLOG(2D0))
C
      IF(NR.GE.NC) THEN
        IF(NV.GE.NR) THEN
          ND=NV
          IG=IGV
          NMAX=NV
        ELSE
          ND=NR
          IG=IGR
          NMAX=NR
        ENDIF
      ELSE
        IF(NV.GE.NC) THEN
          ND=NV
          IG=IGV
          NMAX=NV
        ELSE
          ND=NC
          IG=IGC
          NMAX=NC
        ENDIF
      ENDIF
      Q=PI*2.0D0/DBLE(ND)
      DO 10 IT=1,ND/2
        C(IT)=DCOS(Q*DBLE(IT-1))
        S(IT)=DSIN(Q*DBLE(IT-1))
  10  CONTINUE
      IK(1)=0
      MN=ND/2
      MC=1
      DO 20 IT1=1,IG
        DO 30 IT2=0,MC-1
          IK(IT2+MC+1)=IK(IT2+1)+MN
  30    CONTINUE
        MN=MN/2
        MC=2*MC
  20  CONTINUE
      RETURN
      END
C
C***********************************************************
      SUBROUTINE FFT(IG,ND,NA)
C***********************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (IMS=64)
      COMMON /IFLG2/ QS
      COMMON /IDAT2/ XF(IMS),YF(IMS)
      COMMON /IDAT3/ C(IMS/2),S(IMS/2),IK(IMS)
C
      L2=1
      DO 10 LF=1,IG
        N2=ND/2/L2
        DO 20 MF=1,L2
        DO 20 NF=0,N2-1
          IX=NF*L2*NA+1
          KA=NF+2*N2*(MF-1)+1
          KB=KA+N2
          TR=XF(KA)-XF(KB)
          TJ=YF(KA)-YF(KB)
          XF(KA)=XF(KA)+XF(KB)
          YF(KA)=YF(KA)+YF(KB)
          XF(KB)=TR*C(IX)+TJ*QS*S(IX)
          YF(KB)=TJ*C(IX)-TR*QS*S(IX)
  20    CONTINUE
        L2=L2*2
  10  CONTINUE
      RETURN
      END
C
C**************************************************************
      SUBROUTINE RFFT(NR,NC,NV,NMAX)
C**************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (IMS=64)
      PARAMETER (IMS3=262144)
      COMMON /IFLG1/ IRF,IFF
      COMMON /IFLG2/ QS
      COMMON /IDAT1/ X(IMS3),Y(IMS3)
      COMMON /IDAT2/ XF(IMS),YF(IMS)
      COMMON /IDAT3/ C(IMS/2),S(IMS/2),IK(IMS)
      COMMON /IDAT4/ IGC,IGR,IGV
C
      IG=IGR
      ND=NR
      NA=NMAX/NR
C
      DO 10 IC=1,NC
      DO 10 IV=1,NV
        DO 20 IR=1,NR
          IXN=(IR-1)*NV*NC+(IC-1)*NV+IV
          XF(IR)=X(IXN)
          YF(IR)=Y(IXN)
  20    CONTINUE
        CALL FFT(IG,ND,NA)
        DO 30 IR=1,NR
          IXN=(IR-1)*NV*NC+(IC-1)*NV+IV
          X(IXN)=XF(IK((IR-1)*NA+1)+1)
          Y(IXN)=YF(IK((IR-1)*NA+1)+1)
  30    CONTINUE
  10  CONTINUE
      RETURN
      END
C
C************************************************************
      SUBROUTINE CFFT(NR,NC,NV,NMAX)
C************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (IMS=64)
      PARAMETER (IMS3=262144)
      COMMON /IFLG1/ IRF,IFF
      COMMON /IFLG2/ QS
      COMMON /IDAT1/ X(IMS3),Y(IMS3)
      COMMON /IDAT2/ XF(IMS),YF(IMS)
      COMMON /IDAT3/ C(IMS/2),S(IMS/2),IK(IMS)
      COMMON /IDAT4/ IGC,IGR,IGV
C
      IG=IGC
      ND=NC
      NA=NMAX/NC
C
      DO 10 IR=1,NR
      DO 10 IV=1,NV
        DO 20 IC=1,NC
          IXN=(IR-1)*NV*NC+(IC-1)*NV+IV
          XF(IC)=X(IXN)
          YF(IC)=Y(IXN)
  20    CONTINUE
        CALL FFT(IG,ND,NA)
        DO 30 IC=1,NC
          IXN=(IR-1)*NV*NC+(IC-1)*NV+IV
          X(IXN)=XF(IK((IC-1)*NA+1)+1)
          Y(IXN)=YF(IK((IC-1)*NA+1)+1)
  30    CONTINUE
  10  CONTINUE
      RETURN
      END
C
C***********************************************************
      SUBROUTINE VFFT(NR,NC,NV,NMAX)
C***********************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (IMS=64)
      PARAMETER (IMS3=262144)
      COMMON /IFLG1/ IRF,IFF
      COMMON /IFLG2/ QS
      COMMON /IDAT1/ X(IMS3),Y(IMS3)
      COMMON /IDAT2/ XF(IMS),YF(IMS)
      COMMON /IDAT3/ C(IMS/2),S(IMS/2),IK(IMS)
      COMMON /IDAT4/ IGC,IGR,IGV
C
      IG=IGV
      ND=NV
      NA=NMAX/NV
C
      DO 10 IR=1,NR
      DO 10 IC=1,NC
        DO 20 IV=1,NV
          IXN=(IR-1)*NV*NC+(IC-1)*NV+IV
          XF(IV)=X(IXN)
          YF(IV)=Y(IXN)
  20    CONTINUE
        CALL FFT(IG,ND,NA)
        DO 30 IV=1,NV
          IXN=(IR-1)*NV*NC+(IC-1)*NV+IV
          X(IXN)=XF(IK((IV-1)*NA+1)+1)
          Y(IXN)=YF(IK((IV-1)*NA+1)+1)
  30    CONTINUE
  10  CONTINUE
      RETURN
      END
C
