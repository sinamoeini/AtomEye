C --------------------------------------------------------------------
c     Test Shannon's information theory: S = -\sum_i p_i * log(p_i)
c     by doing binary search on (0-1) till accuracy \epsilon.
c     The average number of inquiries should be S_total / S_search
c     with S_total = -log(\epsilon) = -1/\epsilon*\epsilon*log(\epsilon), 
c     S_search = -\eta*log(\eta)-(1-\eta)*log(1-\eta),
c     where \eta is the ratio of devision.
c     This formula is also derivable using a recursion relation 
c     and linear scaling between epsilon and the searching range.
C --------------------------------------------------------------------
      
      PROGRAM Search
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (EPSILON=0.00001, MESH_ETA=200, NTRIAL=50000)
      external INTEGER time
      data LP /25/

      OPEN (UNIT=LP,STATUS='UNKNOWN',FORM='FORMATTED',FILE='dat.out')
      print *,' \eta   average searching time   theoretical prediction'
      do i=1,mesh_eta
      eta = 0.5d0/mesh_eta*i
      sum = 0
      do j=1,NTRIAL
      sum = sum+NUMBER_OF_SEARCH(rand(),eta,EPSILON)
      enddo
      sum = sum/NTRIAL
      theory = log(EPSILON)/(eta*log(eta)+(1-eta)*log(1-eta))
      write (*,'(1x, F6.4,7x,F6.2,19x,F6.2)') eta, sum, theory
      write (LP,'(F5.3,1x,F6.2,1x,F6.2)') eta, sum, theory
      enddo
      close(lp)
      
      STOP
      END
      

      INTEGER FUNCTION NUMBER_OF_SEARCH (X,ETA,EPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION LOW,MID
      
      LOW = 0.
      HIGH = 1.
      NUMBER_OF_SEARCH = 0
      
 19   IF (HIGH-LOW.GT.EPS) THEN
      NUMBER_OF_SEARCH = NUMBER_OF_SEARCH+1
      MID = LOW + (HIGH-LOW)*ETA
      IF (X.GE.MID) THEN 
      LOW = MID
      ELSE
      HIGH = MID
      ENDIF
      ELSE
      RETURN
      ENDIF
      GOTO 19
      
      END
      
      
      
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
           idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum

      ran1=min(AM*iy,RNMX)
      return
      END
