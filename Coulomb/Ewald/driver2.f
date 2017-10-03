c     ------------------------------------------------
c     driver program to test ewald.c:
c     here we see how good energy is conserved in MD
c      
c               Li Ju, June 18, 1997
c     ------------------------------------------------

      program driver2
      implicit double precision (a-h,o-z)
      parameter (mma=12, EPS=1D-8)
      dimension charge(mma),sx(mma),sy(mma),sz(mma),fx(mma),fy(mma),
     a  fz(mma), a(3,3), h(3,3), stress(3,3), tau(mma,3), t(mma,3),
     a ntype(mma),sx1(mma),sy1(mma),sz1(mma),fsx(mma),fsy(mma),fsz(mma),
     a  sx2(mma), sy2(mma), sz2(mma), sx3(mma), sy3(mma), sz3(mma),
     a  sx4(mma), sy4(mma), sz4(mma), sx5(mma), sy5(mma), sz5(mma)
      double precision kine
      external ewald_potential, total_ewald_energy
      double precision function ewald_potential, total_ewald_energy
      dimension pote_matrix(mma*mma)
            
      delta = 0.001
      cc = delta*delta/2.d0
      np = 2

c     set parameters in predictor-corrector method
      f02=3.d0/16.d0
      f12=251.d0/360.d0
      f32=11.d0/18.d0
      f42=1.d0/6.d0
      f52=1.d0/60.d0
      
      h(1,1) = 0.5
      h(1,2) = 0.5
      h(1,3) = 0.
      h(2,1) = 0.5
      h(2,2) = 0.
      h(2,3) = 0.5
      h(3,1) = 0.
      h(3,2) = 0.5
      h(3,3) = 0.5
      
      sx(1) = 0.
      sy(1) = 0.
      sz(1) = 0.
      charge(1) = 1.
      
      sx(2) = 0.3
      sy(2) = 0.3
      sz(2) = 0.3    
      charge(2) = -1.

      do 380 i = 1,np
      sx1(i) = 0.d0
      sy1(i) = 0.d0
      sz1(i) = 0.d0
      sx2(i) = 0.d0
      sy2(i) = 0.d0
      sz2(i) = 0.d0
      sx3(i) = 0.d0
      sy3(i) = 0.d0
      sz3(i) = 0.d0
      sx4(i) = 0.d0
      sy4(i) = 0.d0
      sz4(i) = 0.d0
      sx5(i) = 0.d0
      sy5(i) = 0.d0
      sz5(i) = 0.d0
380   continue
      
c     get the inverse matrix
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
      
      call init_ewald(np,h,eps)

      kp = 0
 20   continue
      
c     predictor

      do 123 i = 1, np
      sx(i)  = sx(i)+sx1(i)+sx2(i)+sx3(i)+sx4(i)+sx5(i)
      sy(i)  = sy(i)+sy1(i)+sy2(i)+sy3(i)+sy4(i)+sy5(i)
      sz(i)  = sz(i)+sz1(i)+sz2(i)+sz3(i)+sz4(i)+sz5(i)
      sx1(i) = sx1(i)+2.d0*sx2(i)+3.d0*sx3(i)+4.d0*sx4(i)+5.d0*sx5(i)
      sy1(i) = sy1(i)+2.d0*sy2(i)+3.d0*sy3(i)+4.d0*sy4(i)+5.d0*sy5(i)
      sz1(i) = sz1(i)+2.d0*sz2(i)+3.d0*sz3(i)+4.d0*sz4(i)+5.d0*sz5(i)
      sx2(i) = sx2(i)+3.d0*sx3(i)+6.d0*sx4(i)+10.d0*sx5(i)
      sy2(i) = sy2(i)+3.d0*sy3(i)+6.d0*sy4(i)+10.d0*sy5(i)
      sz2(i) = sz2(i)+3.d0*sz3(i)+6.d0*sz4(i)+10.d0*sz5(i)
      sx3(i) = sx3(i)+4.d0*sx4(i)+10.d0*sx5(i)
      sy3(i) = sy3(i)+4.d0*sy4(i)+10.d0*sy5(i)
      sz3(i) = sz3(i)+4.d0*sz4(i)+10.d0*sz5(i)
      sx4(i) = sx4(i)+5.d0*sx5(i)
      sy4(i) = sy4(i)+5.d0*sy5(i)
      sz4(i) = sz4(i)+5.d0*sz5(i)
 123  continue
      
      call ewald (charge,sx,sy,sz,h,pote,fx,fy,fz,stress)
C     test simple_ewald.c
      call simple_ewald (h(1,1),h(1,2),h(1,3),h(2,1),h(2,2),h(2,3),
     A h(3,1),h(3,2),h(3,3),np,sx,sy,sz,eps,pote_matrix)
      simple_pote = total_ewald_energy(1.d0,np,charge,pote_matrix)
      print *, 'pote=',pote, 'simple_pote = ', simple_pote
      if (abs((simple_pote-pote)/pote).gt.eps) then 
      print *, 'sx(1)=', sx(1), 'sy(1)=', sy(1),'sz(1)=', sz(1)
      print *, 'sx(2)=', sx(2), 'sy(2)=', sy(2),'sz(2)=', sz(2)
      print *, 'pote = ',pote, 'kine = ', kine, 'tote = ', pote+kine      
      stop 'simple_pote != pote'
      endif
      
      do 210 i = 1,np
      fsx(i) = (g11*fx(i)+g12*fy(i)+g13*fz(i))/volume
      fsy(i) = (g21*fx(i)+g22*fy(i)+g23*fz(i))/volume
      fsz(i) = (g31*fx(i)+g32*fy(i)+g33*fz(i))/volume
 210  continue

      kine = 0.
      do 330 i=1,np
      sxerr=sx2(i)-cc*fsx(i)
      syerr=sy2(i)-cc*fsy(i)
      szerr=sz2(i)-cc*fsz(i)

      sx(i)  = sx(i) -sxerr*f02
      sx1(i) = sx1(i)-sxerr*f12
      sx2(i) = sx2(i)-sxerr
      sx3(i) = sx3(i)-sxerr*f32
      sx4(i) = sx4(i)-sxerr*f42
      sx5(i) = sx5(i)-sxerr*f52

      sy(i)  = sy(i) -syerr*f02
      sy1(i) = sy1(i)-syerr*f12
      sy2(i) = sy2(i)-syerr
      sy3(i) = sy3(i)-syerr*f32
      sy4(i) = sy4(i)-syerr*f42
      sy5(i) = sy5(i)-syerr*f52
      
      sz(i)  = sz(i) -szerr*f02
      sz1(i) = sz1(i)-szerr*f12
      sz2(i) = sz2(i)-szerr
      sz3(i) = sz3(i)-szerr*f32
      sz4(i) = sz4(i)-szerr*f42
      sz5(i) = sz5(i)-szerr*f52

      x1 = sx1(i)*h(1,1) + sy1(i)*h(2,1) + sz1(i)*h(3,1)
      y1 = sx1(i)*h(1,2) + sy1(i)*h(2,2) + sz1(i)*h(3,2)
      z1 = sx1(i)*h(1,3) + sy1(i)*h(2,3) + sz1(i)*h(3,3)
      kine = kine + (x1**2+y1**2+z1**2)/2./delta/delta
 330  continue
      
      if (mod(kp,20).eq.0) then 
      print *, 'at step ', kp
      print *, 'sx(1)=', sx(1), 'sy(1)=', sy(1),'sz(1)=', sz(1)
      print *, 'sx(2)=', sx(2), 'sy(2)=', sy(2),'sz(2)=', sz(2)
      print *, 'pote = ',pote, 'kine = ', kine, 'tote = ', pote+kine
      print *, ' '
      endif
      
      kp = kp+1
      goto 20

      stop
      end




