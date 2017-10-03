C     generate an LO config with very small k:
C     f77 LO.f; a.out; more conLO; ./driver1.IRIX64 < conLO > 1
      program LO
      implicit double precision (a-h,o-z)
      parameter (KPeriod=32,mma=4,mmp=mma*KPeriod)
      dimension charge(mmp)
      data lp_con /20/
      open (unit=lp_con, status='unknown', form='formatted',
     a      file='conLO')
      
      write(lp_con,'(
     A "Total number of phases to calculate E_Coulomb:",/,"2",)')
      
C     reference state
      write(lp_con,'(/,
     A "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",/,
     A "Name:",/,
     A "3C",/,
     A "Lattice vectors:",/,
     A "0.5 	0.5	0. ",/,
     A "0.5	0.	0.5",/,
     A "0.	0.5	0.5",/,
     A "Number of atoms in the unit cell:",/,
     A "2",/,
     A "Type   (X	Y	Z):",/,
     A "1	0.	0.	0.",/,
     A "2	0.25	0.25	0.25",/,
     A "Number of volumes:",/,
     A "1",/,
     A "Their respective constant of lengths (A) and charge (e):",/,
     A "4.        3.         5.",/,/,
     A "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')

      write(lp_con,'(/,
     A "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",/,
     A "Name:",/,
     A "3C_LO",/,
     A "Lattice vectors:")')
      
      write(lp_con,'(3f14.9)') 0.5,0.5,0.d0
      write(lp_con,'(3f14.9)')-0.5,0.5,0.d0
      write(lp_con,'(3f14.9)') 0.d0,0.d0,1.*KPeriod
      
      write(lp_con,'(
     A "Number of atoms in the unit cell:",/,
     A i6)') mmp
      write (lp_con,'("Type   (X	Y	Z):")')
      amplitude = 0.001
C     sinusoidal wave:
      do i = 0,KPeriod-1
      write (lp_con,'(i2,3f16.12)') 1, 0.d0,0.d0,
     A i     -amplitude*dsin(2*acos(-1.d0)/KPeriod*i)
      write (lp_con,'(i2,3f16.12)') 2, 0.25,0.25,
     A i+0.25+amplitude*dsin(2*acos(-1.d0)/KPeriod*(i+0.25))
      write (lp_con,'(i2,3f16.12)') 1, 0.d0,0.50,
     A i+0.50-amplitude*dsin(2*acos(-1.d0)/KPeriod*(i+0.50))
      write (lp_con,'(i2,3f16.12)') 2, 0.25,0.75,
     A i+0.75+amplitude*dsin(2*acos(-1.d0)/KPeriod*(i+0.75))
      enddo
      
      ip = 1
      do n = 1,KPeriod
      charge(ip) = 3.
      charge(ip+1) = 5.
      charge(ip+2) = 3.
      charge(ip+3) = 5.
      ip = ip+4
      enddo
      
      write (lp_con,'("Number of volumes:",/,
     A "1",/,
     A "Their respective constant of lengths (A) and charge (e):",/,
     A "4. ",10000f3)') (charge(ip),ip=1,mmp)
      
      write (lp_con,'(/,
     A "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
      close(lp_con)
      
      stop
      end
      
