C     Test if 2D dipolar solid has formation energy
      
      PROGRAM Converge
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAX_Y=50000)
      DIMENSION E_LINE_LINE(MAX_Y-1)
      CHARACTER *70 BUF
      
      READ (*,'(A70)') buf
c     'aspect ratio:'
      read *, aspect_ratio
      READ (*,'(/,A70)') buf
c     'mesh points in x-direction:'
      read *, mesh_x
c     make it even
      mesh_x = mesh_x/2 * 2 
      mesh_y = aspect_ratio * mesh_x

c     calculate the self-energy of a dipolar line
      E_LINE_SELF = 0.
      do idx = 1, mesh_x-1
      idx2 = idx * idx
      E_LINE_SELF = E_LINE_SELF - 1.d0 / idx2 * (mesh_x-idx)
      enddo
      
c     calculate line-line interaction energy
      do 432 idy = 1, mesh_y-1
      idy2 = idy * idy
      E_LINE_LINE(idy) = 0.
      do 433 ix1 = 1, mesh_x/2
      do 433 ix2 = 1, mesh_x
      idx2 = (ix2-ix1) * (ix2-ix1)
      ir2 = idx2 + idy2 
      E_LINE_LINE(idy) = E_LINE_LINE(idy) + dble(idy2-idx2)/ir2/ir2
 433  continue
      E_LINE_LINE(idy) = E_LINE_LINE(idy) * 2.
 432  continue
      
      total = mesh_y * E_LINE_SELF
c     1D summation 
      do 543 idy = 1, mesh_y-1
 543  total = total + E_LINE_LINE(idy) * (mesh_y-idy)

      print *,'self energy of line =', E_LINE_SELF/mesh_x
      print *,total/mesh_x/mesh_y
      
      stop
      end

