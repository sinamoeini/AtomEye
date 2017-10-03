c     ----------------------------------
c     driver program to test ewald.c   
c      
c               Li Ju, June 18, 1997
c     ----------------------------------

      program driver
      implicit double precision (a-h,o-z)
      parameter (mma=256, EPS=1D-10, mmd=3*mma)
      double precision charge(mma),s1(mma),s2(mma),s3(mma),fx(mma),
     a  fy(mma), fz(mma), a(3,3), h(3,3), stress(3,3), tau(mma,3), 
     a  ntype(mma), kpt(3), t(mma,3), mass(mma)
      complex*16 phi(mmd,mmd),ap(mmd*(mmd+1)/2),add,
     a               zauz(2*mmd),z(mmd,mmd)
      double precision auz(6*mmd), w(mmd)
      character *70 buf,name
      data lp_ec /20/
      external init_ewald, ewald, exit_ewald, ewald_dynamical
      external simple_ewald
      external ewald_potential, total_ewald_energy
      double precision function ewald_potential, total_ewald_energy
      dimension pote_matrix(mma*mma)
      
      PI = dacos(-1.d0)
      mass(1) = 28.
      mass(2) = 12.
      
      read (*,'(a50)') buf
c     "Total number of phases to calculate E_Coulomb:"
      read *, nphase

      do 642 iphase=1,nphase

      read (*,'(/,A70,/,A70)') buf,buf
C     "Name:"
      read (*,'(A50)') name

      buf = name(1:index(name,' ')-1) // '_Ec.out'
      open (unit=lp_ec, status='unknown', form='formatted',
     a      file=buf)
      read (*,'(a50)') buf
c     "Lattice vectors:"
      read *,((a(i,j),j=1,3),i=1,3)
      read (*,'(a50)') buf
c     "number of atoms in the unit cell:"
      read *,npa
      read (*,'(a50)') buf
c     "type      x         y          z"
      read *, (ntype(n),(t(n,i),i=1,3),n=1,npa)
      read (*,'(a50)') buf
c     "Number of volumes:"
      read *, nvolume
      read (*,'(a50)') buf
c     "their respective constant of lengths (a) and charges (e):"

      do 431 ivolume = 1, nvolume
      read *,aa,(charge(n),n=1,npa)

      do i=1,3
      do j=1,3
      h(i,j) = a(i,j)*aa
      enddo
      do n=1,npa
      tau(n,i) = t(n,i)*aa
      enddo
      enddo
      
c     Make sure that the total charge sum to 0:
      sumcharge = 0.d0
      do 12 i = 1, npa
 12   sumcharge = sumcharge + charge(i) 
      do 13 i = 1, npa
 13   charge(i) = charge(i) - sumcharge/npa
      
      if (ivolume.eq.1) then 
      print *, 'calculating E_Coulomb for ',name(1:index(name,' ')-1), 
     a         ' with EPS = ',EPS
      print *, '*********************************************'
      call init_ewald(npa,h,eps)
      endif

c     get the volume and the reduced coordinates of atoms
      g11=h(2,2)*h(3,3)-h(2,3)*h(3,2)
      g22=h(3,3)*h(1,1)-h(3,1)*h(1,3)
      g33=h(1,1)*h(2,2)-h(1,2)*h(2,1)
      g12=h(2,3)*h(3,1)-h(2,1)*h(3,3)
      g23=h(3,1)*h(1,2)-h(3,2)*h(1,1)
      g31=h(1,2)*h(2,3)-h(1,3)*h(2,2)
      g13=h(2,1)*h(3,2)-h(3,1)*h(2,2)
      g21=h(3,2)*h(1,3)-h(1,2)*h(3,3)
      g32=h(1,3)*h(2,1)-h(2,3)*h(1,1)
      volume=h(1,1)*g11+h(1,2)*g12+h(1,3)*g13
      
      do i=1,npa
      s1(i) = (g11*tau(i,1)+g12*tau(i,2)+g13*tau(i,3))/volume
      s2(i) = (g21*tau(i,1)+g22*tau(i,2)+g23*tau(i,3))/volume
      s3(i) = (g31*tau(i,1)+g32*tau(i,2)+g33*tau(i,3))/volume
      enddo
      
      volume=abs(volume)
      
      call ewald(charge,s1,s2,s3,h,pote,fx,fy,fz,stress)
      print *, ' '
      print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *, 'AA = ', AA, 'pote = ', pote

C     test simple_ewald.c
      call simple_ewald(h(1,1),h(1,2),h(1,3),h(2,1),h(2,2),h(2,3),
     A h(3,1),h(3,2),h(3,3),npa,s1,s2,s3,eps,pote_matrix)
      if (npa.le.4) print *, 'coeff. matrix = ', 
     A (pote_matrix(i),i=1,npa*npa)
      simple_pote = total_ewald_energy(1.d0,npa,charge,pote_matrix)
      print *, 'simple_pote = ', simple_pote
      if (abs((simple_pote-pote)/pote).gt.eps) 
     A stop 'simple_pote != pote'
      print *, 'potential field = ', 
     A (ewald_potential(n,1.d0,npa,charge,pote_matrix),n=1,npa)
      
c     sanity check No.1: Virial Theorem
      err = (stress(1,1)+stress(2,2)+stress(3,3))*volume-pote
      if (abs(err).gt.eps) then 
      print *, 'Failed sanity check No.1 (Virial Theorem) by ', err
      do i=1,npa
      print *, 'Force(x,y,z) on atom ',i,': ',fx(i),fy(i),fz(i)
      enddo
      print *, ' ' 
      endif

      if (npa.eq.2) then 
c     sanity check No.2: Blackman's sum rule; 
c     <Theory of Lattice Dynamics> Page 230, (6.3.37) and Page 225, (6.3.18)
c     kpt(1) = 0.432
c     kpt(2) = 1.654
c     kpt(3) = -0.67
c     sanity check No.3: LO-TO splitting;
c     <Electronic Structure and Properties of Solids> Page 219, (9-22)
      kpt(1) = 0.
      kpt(2) = 0.
      kpt(3) = 0.
      call ewald_dynamical (kpt,phi,mmd)
      lenai=0
      do 1070 jj=1,3*npa
      jp=(jj-1)/3+1
      do 1070 ii=1,jj
      ip=(ii-1)/3+1
      lenai=lenai+1
      ap(lenai)=phi(ii,jj)/sqrt(mass(ntype(ip))*mass(ntype(jp)))
 1070 continue
#ifdef _AIX      
c     ESSL driver
      CALL ZHPEV(21, AP, W, Z, MMD, 3*npa, AUZ, 6*MMD)
#else
c     LAPACK driver
      INFO = 0
      CALL ZHPEV('V', 'U', 3*npa, AP, W, Z, MMD, ZAUZ, AUZ, INFO)
#endif
      print *, 'The force constant matrix = '
      write (*,422) 
     a ((real(phi(i,j)),imag(phi(i,j)),j=1,3*npa),i=1,3*npa)
 422  format(6(6('[',f7.4,',',f7.4,']'),/))
      print *, 'The eigenvalues are'
      write (*,'(8(2x,f10.7))') (w(i),i=1,3*npa)
      sum = 0.
      do i=1,3*npa
      sum = sum+w(i)
      enddo
      print *, 'Their sum = ', sum
      reduced = 1. / ( 1./mass(1) + 1./mass(2) ) 
      print *, 'LO-TO splitting should = ',
     a         4*PI*charge(1)**2/reduced/2.*npa/volume
      print *, 'The actual splitting = ', w(3*npa)-w(1)
      endif
      
c     calculate the nearest neighbour distance to the first atom
      rmin2 = 1000. 
      do 654 n=2,npa
      dx = tau(n,1)-tau(1,1)
      dy = tau(n,2)-tau(1,2)
      dz = tau(n,3)-tau(1,3)
      do 654 i=-1,1
      do 654 j=-1,1
      do 654 k=-1,1
      rx = dx + i*h(1,1) + j*h(2,1) + k*h(3,1)
      ry = dy + i*h(1,2) + j*h(2,2) + k*h(3,2)
      rz = dz + i*h(1,3) + j*h(2,3) + k*h(3,3)
      r2 = rx*rx+ry*ry+rz*rz
      rmin2 = min(rmin2, r2)
 654  continue
      
      print *, 'Volume/atom(A^3)    Ec/atom(eV)   Madelung constant'
      write (*,'(3(2x,f11.7,4x),/)') volume/npa, pote/npa,
     a              -pote/npa*2*sqrt(rmin2)/charge(1)**2
      write (lp_ec, '(3(1x,f11.7,4x))') volume/npa, pote/npa,
     a              -pote/npa*2*sqrt(rmin2)/charge(1)**2
      print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print *,' '
 431  continue
      call exit_ewald
      
      print *, '*********************************************'
      print *,' '
      close (lp_ec)
      read (*,'(/,a50)') buf
c     "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
 642  continue

      stop
      end


