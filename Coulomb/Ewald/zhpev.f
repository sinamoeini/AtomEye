* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module ZHPEV from package LAPACK.
* Retrieved from NETLIB on Sun Jun 14 17:57:20 1998.
* ======================================================================
      SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix in packed storage.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the Hermitian matrix
*          A, packed columnwise in a linear array.  The j-th column of A
*          is stored in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*
*          On exit, AP is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
*          and first superdiagonal of the tridiagonal matrix T overwrite
*          the corresponding elements of A, and if UPLO = 'L', the
*          diagonal and first subdiagonal of T overwrite the
*          corresponding elements of A.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (max(1, 2*N-1))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK,
     $                   ISCALE
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHP
      EXTERNAL           LSAME, DLAMCH, ZLANHP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZDSCAL, ZHPTRD, ZSTEQR,
     $                   ZUPGTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) )
     $          THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = AP( 1 )
         RWORK( 1 ) = 1
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = ZLANHP( 'M', UPLO, N, AP, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL ZDSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      END IF
*
*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      CALL ZHPTRD( UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ),
     $             IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     ZUPGTR to generate the orthogonal matrix, then call ZSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         INDWRK = INDTAU + N
         CALL ZUPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ,
     $                WORK( INDWRK ), IINFO )
         INDRWK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ,
     $                RWORK( INDRWK ), INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      RETURN
*
*     End of ZHPEV
*
      END
      SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         AP( * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHPTRD reduces a complex Hermitian matrix A stored in packed form to
*  real symmetric tridiagonal form T by a unitary similarity
*  transformation: Q**H * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the Hermitian matrix
*          A, packed columnwise in a linear array.  The j-th column of A
*          is stored in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the unitary
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the unitary matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) COMPLEX*16 array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
*  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
*  overwriting A(i+2:n,i), and tau is stored in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   HALF = ( 0.5D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, I1, I1I1, II
      COMPLEX*16         ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHPMV, ZHPR2, ZLARFG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        I1 is the index in AP of A(1,I+1).
*
         I1 = N*( N-1 ) / 2 + 1
         AP( I1+N-1 ) = DBLE( AP( I1+N-1 ) )
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            ALPHA = AP( I1+I-1 )
            CALL ZLARFG( I, ALPHA, AP( I1 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               AP( I1+I-1 ) = ONE
*
*              Compute  y := tau * A * v  storing y in TAU(1:i)
*
               CALL ZHPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU,
     $                     1 )
*
*              Compute  w := y - 1/2 * tau * (y'*v) * v
*
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, AP( I1 ), 1 )
               CALL ZAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL ZHPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )
*
            END IF
            AP( I1+I-1 ) = E( I )
            D( I+1 ) = AP( I1+I )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = AP( 1 )
      ELSE
*
*        Reduce the lower triangle of A. II is the index in AP of
*        A(i,i) and I1I1 is the index of A(i+1,i+1).
*
         II = 1
         AP( 1 ) = DBLE( AP( 1 ) )
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            ALPHA = AP( II+1 )
            CALL ZLARFG( N-I, ALPHA, AP( II+2 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               AP( II+1 ) = ONE
*
*              Compute  y := tau * A * v  storing y in TAU(i:n-1)
*
               CALL ZHPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1,
     $                     ZERO, TAU( I ), 1 )
*
*              Compute  w := y - 1/2 * tau * (y'*v) * v
*
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, AP( II+1 ),
     $                 1 )
               CALL ZAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL ZHPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1,
     $                     AP( I1I1 ) )
*
            END IF
            AP( II+1 ) = E( I )
            D( I ) = AP( II )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = AP( II )
      END IF
*
      RETURN
*
*     End of ZHPTRD
*
      END
      DOUBLE COMPLEX   FUNCTION ZLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      COMPLEX*16         X, Y
*     ..
*
*  Purpose
*  =======
*
*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX*16
*  Y       (input) COMPLEX*16
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
*     ..
*     .. Executable Statements ..
*
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR,
     $             ZI )
      ZLADIV = DCMPLX( ZR, ZI )
*
      RETURN
*
*     End of ZLADIV
*
      END
      DOUBLE PRECISION FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         AP( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLANHP  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  complex hermitian matrix A,  supplied in packed form.
*
*  Description
*  ===========
*
*  ZLANHP returns the value
*
*     ZLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in ZLANHP as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          hermitian matrix A is supplied.
*          = 'U':  Upper triangular part of A is supplied
*          = 'L':  Lower triangular part of A is supplied
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHP is
*          set to zero.
*
*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The upper or lower triangle of the hermitian matrix A, packed
*          columnwise in a linear array.  The j-th column of A is stored
*          in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*          Note that the  imaginary parts of the diagonal elements need
*          not be set and are assumed to be zero.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            K = 0
            DO 20 J = 1, N
               DO 10 I = K + 1, K + J - 1
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   10          CONTINUE
               K = K + J
               VALUE = MAX( VALUE, ABS( DBLE( AP( K ) ) ) )
   20       CONTINUE
         ELSE
            K = 1
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( AP( K ) ) ) )
               DO 30 I = K + 1, K + N - J
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   30          CONTINUE
               K = K + N - J + 1
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is hermitian).
*
         VALUE = ZERO
         K = 1
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( AP( K ) ) )
               K = K + 1
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( AP( K ) ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         K = 2
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, AP( K ), 1, SCALE, SUM )
               K = K + J
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
  120       CONTINUE
         END IF
         SUM = 2*SUM
         K = 1
         DO 130 I = 1, N
            IF( DBLE( AP( K ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( AP( K ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
            IF( LSAME( UPLO, 'U' ) ) THEN
               K = K + I + 1
            ELSE
               K = K + N - I + 1
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      ZLANHP = VALUE
      RETURN
*
*     End of ZLANHP
*
      END
      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARF applies a complex elementary reflector H to a complex M-by-N
*  matrix C, from either the left or the right. H is represented in the
*  form
*
*        H = I - tau * v * v'
*
*  where tau is a complex scalar and v is a complex vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
*  tau.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) COMPLEX*16 array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) COMPLEX*16
*          The value tau in the representation of H.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL ZGEMV( 'Conjugate transpose', M, N, ONE, C, LDC, V,
     $                  INCV, ZERO, WORK, 1 )
*
*           C := C - v * w'
*
            CALL ZGERC( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL ZGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL ZGERC( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of ZLARF
*
      END
      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFG generates a complex elementary reflector H of order n, such
*  that
*
*        H' * ( alpha ) = ( beta ),   H' * H = I.
*             (   x   )   (   0  )
*
*  where alpha and beta are scalars, with beta real, and x is an
*  (n-1)-element complex vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a complex scalar and v is a complex (n-1)-element
*  vector. Note that H is not hermitian.
*
*  If the elements of x are all zero and alpha is real, then tau = 0
*  and H is taken to be the unit matrix.
*
*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) COMPLEX*16
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) COMPLEX*16
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
*
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END
      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASET initializes a 2-D array A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set. The lower triangle
*                      is unchanged.
*          = 'L':      Lower triangular part is set. The upper triangle
*                      is unchanged.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          On entry, M specifies the number of rows of A.
*
*  N       (input) INTEGER
*          On entry, N specifies the number of columns of A.
*
*  ALPHA   (input) COMPLEX*16
*          All the offdiagonal array elements are set to ALPHA.
*
*  BETA    (input) COMPLEX*16
*          All the diagonal array elements are set to BETA.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
*                   A(i,i) = BETA , 1 <= i <= min(m,n)
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASET
*
      END
      SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n complex matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) DOUBLE PRECISION arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP
      COMPLEX*16         TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZLASR
*
      END
      SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASSQ returns the values scl and ssq such that
*
*     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
*  assumed to be at least unity and the value of ssq will then satisfy
*
*     1.0 .le. ssq .le. ( sumsq + 2*n ).
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*            i
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  SCALE and SUMSQ are overwritten by scl and ssq respectively.
*
*  The routine makes only one pass through the vector X.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION
*          The vector x as described above.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with the value  scl .
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with the value  ssq .
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASSQ
*
      END
      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band complex Hermitian matrix can also
*  be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  Hermitian matrix.  On entry, Z must contain the
*                  unitary matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the unitary
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original Hermitian matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is unitarily similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA,
     $                   ZLASET, ZLASR, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL ZLASET( 'F', N, N, CZERO, CONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL ZLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL ZLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL ZLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL ZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL DLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
*
*     End of ZSTEQR
*
      END
      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
*  which is defined as the last n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by ZGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQLF in the last k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQLF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL ZLARF( 'L', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL ZSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2L
*
      END
      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by ZGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m by n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQRF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'L', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2R
*
      END
      SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDQ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUPGTR generates a complex unitary matrix Q which is defined as the
*  product of n-1 elementary reflectors H(i) of order n, as returned by
*  ZHPTRD using packed storage:
*
*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangular packed storage used in previous
*                 call to ZHPTRD;
*          = 'L': Lower triangular packed storage used in previous
*                 call to ZHPTRD.
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The vectors which define the elementary reflectors, as
*          returned by ZHPTRD.
*
*  TAU     (input) COMPLEX*16 array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZHPTRD.
*
*  Q       (output) COMPLEX*16 array, dimension (LDQ,N)
*          The N-by-N unitary matrix Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N-1)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IJ, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNG2L, ZUNG2R
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUPGTR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to ZHPTRD with UPLO = 'U'
*
*        Unpack the vectors which define the elementary reflectors and
*        set the last row and column of Q equal to those of the unit
*        matrix
*
         IJ = 2
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   10       CONTINUE
            IJ = IJ + 2
            Q( N, J ) = CZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            Q( I, N ) = CZERO
   30    CONTINUE
         Q( N, N ) = CONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL ZUNG2L( N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to ZHPTRD with UPLO = 'L'.
*
*        Unpack the vectors which define the elementary reflectors and
*        set the first row and column of Q equal to those of the unit
*        matrix
*
         Q( 1, 1 ) = CONE
         DO 40 I = 2, N
            Q( I, 1 ) = CZERO
   40    CONTINUE
         IJ = 3
         DO 60 J = 2, N
            Q( 1, J ) = CZERO
            DO 50 I = J + 1, N
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   50       CONTINUE
            IJ = IJ + 2
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL ZUNG2R( N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK,
     $                   IINFO )
         END IF
      END IF
      RETURN
*
*     End of ZUPGTR
*
      END
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   E, F
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
*
      RETURN
*
*     End of DLADIV
*
      END
      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
*     ..
*
*  Purpose
*  =======
*
*  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, and RT2
*  is the eigenvalue of smaller absolute value.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) and (2,1) elements of the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of DLAE2
*
      END
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) DOUBLE PRECISION
*  SN1     (output) DOUBLE PRECISION
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of DLAEV2
*
      END
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANST  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real symmetric tridiagonal matrix A.
*
*  Description
*  ===========
*
*  DLANST returns the value
*
*     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANST as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
*          set to zero.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of A.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) sub-diagonal or super-diagonal elements of A.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $              ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+
     $                 ABS( E( I-1 ) ) )
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      DLANST = ANORM
      RETURN
*
*     End of DLANST
*
      END
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
*     ..
*
*  Purpose
*  =======
*
*  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*  Z       (input) DOUBLE PRECISION
*          X, Y and Z specify the values x, y and z.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         DLAPY3 = ZERO
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY3
*
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of DLARTG
*
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END
      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
*     ..
*
*  Purpose
*  =======
*
*  Sort the numbers in D in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the array to be sorted.
*          On exit, D has been sorted into increasing order
*          (D(1) <= ... <= D(N) ) or into decreasing order
*          (D(1) >= ... >= D(N) ), depending on ID.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
      SUBROUTINE DSTERF( N, D, E, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm failed to find all of the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDM1, LENDP1,
     $                   LENDSV, LM1, LSV, M, MM1, NM1, NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 60 M = L, LENDM1
               TST = ABS( E( M ) )
               IF( TST.LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
*
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         MM1 = M - 1
         DO 80 I = MM1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 110 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )
               IF( TST.LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $            GO TO 120
  110       CONTINUE
         END IF
*
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         LM1 = L - 1
         DO 130 I = M, LM1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( LM1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
*
      END IF
*
*     Undo scaling if necessary
*
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 160 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  160    CONTINUE
         RETURN
      END IF
      GO TO 10
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
*
      RETURN
*
*     End of DSTERF
*
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module DSCAL from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 14:50:56 1998.
* ======================================================================
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
            
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module ZHPMV from package BLAS2.
* Retrieved from NETLIB on Mon Jun 15 14:52:24 1998.
* ======================================================================
      SUBROUTINE ZHPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX*16         AP( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHPMV  performs the matrix-vector operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n hermitian matrix, supplied in packed form.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  AP     - COMPLEX*16       array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
*           and a( 2, 2 ) respectively, and so on.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
*           and a( 3, 1 ) respectively, and so on.
*           Note that the imaginary parts of the diagonal elements need
*           not be set and are assumed to be zero.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DBLE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 6
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZHPMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when AP contains the upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               K     = KK
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + DCONJG( AP( K ) )*X( I )
                  K      = K      + 1
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*DBLE( AP( KK + J - 1 ) )
     $                         + ALPHA*TEMP2
               KK     = KK     + J
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, K = KK, KK + J - 2
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + DCONJG( AP( K ) )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*DBLE( AP( KK + J - 1 ) )
     $                           + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + J
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when AP contains the lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J ) + TEMP1*DBLE( AP( KK ) )
               K      = KK     + 1
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + DCONJG( AP( K ) )*X( I )
                  K      = K      + 1
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
               KK     = KK     + ( N - J + 1 )
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY ) + TEMP1*DBLE( AP( KK ) )
               IX      = JX
               IY      = JY
               DO 110, K = KK + 1, KK + N - J
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + DCONJG( AP( K ) )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + ( N - J + 1 )
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZHPMV .
*
      END
      

      SUBROUTINE SMESSG(NUNIT,IP,NMESS) 
C     DEFINE THE TEXT OF ERROR MESSAGES.
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      LOGICAL EX,OP,INQ
      INTEGER NUNIT(3)
      CHARACTER *256 MESS(09) 
      CHARACTER *256 M
      DATA MESS(01)/
     .' THE COMPILER IS GENERATING BAD CODE FOR IN-LINE DOT PRODUCTS OR
     .IS INCORRECTLY EVALUATING THE ARITHMETIC EXPRESSIONS J*((J+1)*J)/2
     . - (J+1)*J*(J-1)/3, J=1 THRU 32.'/
      DATA MESS(02)/
     .' ABNORMAL OR EARLY END-OF-FILE WHILE READING NAME OF FILE THAT CO
     .NTAINS THE NAMES OF THE SUBPROGRAMS AND THE SUMMARY FILES.'/
      DATA MESS(03)/
     .' THE ABOVE FILE NAME MUST BE PRESENT ON THE SYSTEM.  IT IS NOT.
     .THIS FILE CONTAINS THE NAMES OF THE SUBPROGRAMS AND THE SUMMARY FI
     .LES.'/
      DATA MESS(04)/
     .' ABNORMAL OR EARLY END-OF-FILE WHILE READING NAMES OF SUBPROGRAMS
     . FROM THE ABOVE FILE NAME.'/
      DATA MESS(05)/
     .' ABNORNAL OR EARLY END-OF-FILE WHILE READING NAMES OF FILES FOR S
     .UMMARY OUTPUT.'/
      DATA MESS(06)/
     .' ENTER NAME AND UNIT NUMBER OF FILE CONTAINING NAMES OF SUBPROGRA
     .MS AND SUMMARY FILES.  ONE ITEM PER LINE, PLEASE.'/
      DATA MESS(07)/
     .' THE SNAP-SHOT FILE OF ACTIVE TESTS CANNOT BE OPENED WITH ''NEW''
     . STATUS OR IT CANNOT BE DELETED.  THIS FILE SHOULD NOT BE PRESENT
     .ON THE SYSTEM.'/
      DATA MESS(08)/
     .' THE SUMMARY FILE OF ACTIVE TESTS CANNOT BE OPENED WITH ''UNKOWN'
     .' STATUS. THIS FILE SHOULD NOT BE PRESENT ON THE SYSTEM.'/
      M = MESS(NMESS)
      NL = 256
      NS = 72
      INQ = .TRUE.
      DO 10 I = NL,1,-1
         IF (ICHAR(M(I:I)).NE.ICHAR(' ')) GO TO 20
   10 CONTINUE
      NL = 0
      GO TO 30
*
   20 NL = I
C     FOUND NS = POINTER TO LAST NONBLANK IN MESSAGE.
   30 CONTINUE
C     NOW OUTPUT THE MESSAGE.  PARSE IT SO THAT UP TO NS CHARS. PER LINE
C     PRINT, BUT DO NOT BREAK WORDS ACCROSS LINES.
      IS = 1
   40 CONTINUE
      IE = MIN(NL,IS+NS)
      IF (IS.GE.IE) GO TO 70
   50 CONTINUE
      IF (ICHAR(M(IE:IE)).EQ.ICHAR(' ') .OR. NL-IS.LT.NS) GO TO 60
      IE = IE - 1
      IF (IE.GT.IS) GO TO 50
   60 CONTINUE
      IF (INQ) THEN 
          INQUIRE (UNIT=NUNIT(IP),EXIST=EX,OPENED=OP)
      END IF
C     IF THE INTENDED UNIT IS NOT OPENED, SEND OUTPUT TO
C     STANDARD OUTPUT SO IT WILL BE SEEN.
      IF ( .NOT. OP .OR. .NOT. EX .OR. NUNIT(IP).EQ.0) THEN 
          IF (IE.EQ.NL) THEN
              WRITE (*,'(A,/)') M(IS:IE)
*
          ELSE
              WRITE (*,'(A)') M(IS:IE)
          END IF
*
          INQ = .FALSE.
*
      ELSE
          LUNIT = NUNIT(IP)
          WRITE (LUNIT,'(A)') M(IS:IE)
      END IF
*
      IS = IE
      GO TO 40
*
   70 CONTINUE
      RETURN
      END 
*
      SUBROUTINE SCHCK1(ISNUM,SNAME,IG,DOPE,NUNIT,AVIGR,FATAL)
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C     THIS PROGRAM HAS TWO PARTS.  THE FIRST PART CHECKS TO SEE
C     IF ANY DATA GETS CHANGED WHEN NONE SHOULD.  FOR EXAMPLE WHEN
C     USING AN INVALID OPTION OR NONPOSITVE PROBLEM DIMENSIONS.
C     THE SECOND PART CHECKS THAT THE RESULTS ARE REASONABLY ACCURATE.
C     DIMENSION AND PROBLEM SIZE DATA.. 
      INTEGER INC(04),IDIM(08),NUNIT(2) 
      REAL ALF(04),BET(04),SDIFF
      LOGICAL ISAME(13),LSE,FATAL,SAME,NCHNG,RESET
      CHARACTER *128 DOPE(2)
      CHARACTER *6 SNAME
      CHARACTER *3 ICH
      CHARACTER *1 ICHS,ICI
      INTEGER LA,LV 
      PARAMETER (LA=4096,LV=4096,LMN=2048)
      REAL A(LA),AS(LA),X(LV),XS(LV)
      REAL Y(LV),YS(LV),YT(LMN),XT(LMN) 
      REAL ALPHA,ALS,BETA,BLS,T,TRANSL,XN
      PARAMETER (ZERO=0.E0,HALF=.5E0,ONE=1.E0)
      INTEGER LIJU_ZERO(3)
      DATA LIJU_ZERO/3*0/
      COMMON /ARRAYS/AR,AS,X,XS,Y,YS,YT,XT
      EXTERNAL SDIFF
*
      DATA ALF/-1.E0,2.E0,.3E0,1.E0/
      DATA BET/-1.E0,0.E0,.3E0,1.E0/
      DATA INC/-2,-1,1,2/
      DATA IDIM/1,2,4,8,64,128,2048,0/
      DATA ICH/'NT/'/
      FATAL = .FALSE.
C     CHECK GENERAL MATRIX-VECTOR PRODUCT, Y = ALPHA*A*X+BETA*Y, NO.1-2.
      IF (ISNUM.LT.0) GO TO 220
      NC = 0
      RESET = .TRUE.
      AVIGR = ZERO
      IX = 0
   10 IX = IX + 1
      IF (IX.GT.4) GO TO 200
      INCX = INC(IX)
      ALPHA = ALF(IX)
      IY = 0
   20 IY = IY + 1
      IF (IY.GT.4) GO TO 190
      INCY = INC(IY)
      BETA = BET(IY)
      MM = 0
   30 MM = MM + 1
      IF (MM.GT.8) GO TO 180
      M = IDIM(MM)
      NN = 0
   40 NN = NN + 1
      IF (NN.GT.8) GO TO 170
      N = IDIM(NN)
      IC = 0
   50 IC = IC + 1
      IF (IC.GT.3) GO TO 160
      IF (FATAL) GO TO 210
C     SET DEFAULT BANDWIDTH SO PRINTING WILL BE OK.
      KL = MAX(0,M-1)
      KU = MAX(0,N-1)
C     DEFINE THE NUMBER OF ARGUMENTS AND THE Y ARGUMENT NUMBER.
      IF (ISNUM.EQ.1) THEN
          LDA = MAX(M,1)
          NARGS = 11
          IYARG = 10
*
      ELSE IF (ISNUM.EQ.2) THEN
          NARGS = 13
          IYARG = 12
C     DEFINE BANDWIDTH OF MATRIX FOR TEST OF SGBMV.
          KL = MAX(0,MIN(M-1,M/2))
          KU = MAX(0,MIN(N-1,N/2))
          LDA = MAX(KL+KU+1,M)
      END IF
*
      ICI = ICH(IC:IC)
      IF (ICHAR(ICI).EQ.ICHAR('T')) THEN
          ML = N
          NL = M
          INCCA = 1 
          INCRA = LDA
*
      ELSE
          ML = M
          NL = N
          INCCA = LDA
          INCRA = 1 
      END IF
*
C     IF NOT ENOUGH STORAGE, SKIP THIS CASE. (AVOID EXPLICT LDA*N).
      IF (SQRT(REAL(N))*SQRT(REAL(LDA)).GT.SQRT(REAL(LA))) GO TO 50
C     DO (PREPARE NOTES FOR THIS TEST)
C
C     PRINT A MESSAGE THAT GIVES DEBUGGING INFORMATION.  THIS
C     MESSAGE SAYS..
C     IN SUBPROGRAM XXXXX TESTS WERE ACTIVE WITH
C     OPTION = 'A'
C      M = IIII,     N = IIII,
C     INCX = IIII,  INCY = IIII,
C     KL = IIII,    KU = IIII.
C     THE MAIN IDEA HERE IS THAT A SERIOUS BUG THAT OCCURS DURING THE 
C     TESTING OF THESE SUBPROGRAMS MAY LOSE PROGRAM CONTROL.  THIS
C     AUXILLIARY FILE CONTAINS THE DIMENSIONS THAT RESULTED IN THE LOSS
C     OF CONTROL.  HENCE IT PROVIDES THE IMPLEMENTOR WITH MORE COMPLETE
C     INFORMATION ABOUT WHERE TO START TRACKING DOWN THE BUG.
      IF (NUNIT(1).GT.0) THEN 
C     IF UNIT IS NOT AVAILABLE WITH 'NEW' STATUS, OPEN WITH 
C     'OLD' AND THEN DELETE IT.
          ISTAT = 1 
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.1) GO TO 60
C     GET RID OF ANY OLD FILE.
          CLOSE (UNIT=NUNIT(1),STATUS='DELETE',ERR=60)
   60     CONTINUE
          ISTAT = 2 
C     CREATE A NEW FILE FOR THE NEXT TEST.
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.0) GO TO 80
          NMESS = 7 
C     DO (PRINT A MESSAGE)
C     PRINT AN ERROR MESSAGE ABOUT OPENING THE NAME FILE.
          CALL SMESSG(LIJU_ZERO,1,NMESS)
          FATAL = .TRUE.
          GO TO 210 
*
   80     CONTINUE
          WRITE (NUNIT(1),9001) SNAME,ICI,M,N,INCX,INCY,KL,KU
C     CLOSE THE FILE SO USEFUL STATUS INFORMATION IS SEALED.
          CLOSE (UNIT=NUNIT(1))
      END IF
C     DO (DEFINE A SET OF PROBLEM DATA) 
      ASSIGN 90 TO IGO3
      GO TO 340
*
   90 CONTINUE
C     DO (CALL SUBROUTINE)
      ASSIGN 100 TO IGO1
      GO TO 280
*
  100 CONTINUE
      IF (M.LE.0 .OR. N.LE.0 .OR. ICHAR(ICI).EQ.ICHAR('/')) THEN
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 110 TO IGO2
          GO TO 240 
*
  110     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 120 I = 1,NARGS
             SAME = SAME .AND. ISAME(I) 
             IF ( .NOT. ISAME(I)) THEN
                 WRITE (NUNIT(2),9011) SNAME,I,ICI,M,N,INCX,INCY,KL,KU
             END IF 
*
  120     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 210
*
          END IF
*
      ELSE
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 130 TO IGO2
          GO TO 240 
*
  130     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 140 I = 1,NARGS
             NCHNG = (I.EQ.IYARG .OR. ISAME(I))
             SAME = SAME .AND. NCHNG
             IF ( .NOT. NCHNG) THEN
                 WRITE (NUNIT(2),9021) SNAME,I,ICI,M,N,INCX,INCY,KL,KU
             END IF 
*
  140     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 210
*
          END IF
*
          NC = NC + 1
C     DO (COMPUTE A CORRECT RESULT)
          ASSIGN 150 TO IGO4
          GO TO 370 
*
  150     CONTINUE
C     IF GOT REALLY BAD ANSWER, PRINT NOTE AND GO BACK.
          IF (FATAL) GO TO 200
*
      END IF
*
      GO TO 50
*
  160 CONTINUE
      GO TO 40
*
  170 CONTINUE
      GO TO 30
*
  180 CONTINUE
      GO TO 20
*
  190 CONTINUE
      GO TO 10
*
  200 CONTINUE
C     REPORT ON ACCURACY OF DATA.
      WRITE (NUNIT(2),9031) ISNUM,SNAME,AVIGR,IG
      GO TO 230
*
  210 CONTINUE
      WRITE (NUNIT(2),9041) ISNUM,SNAME 
      GO TO 230
*
  220 CONTINUE
      WRITE (NUNIT(2),9051) - ISNUM,SNAME
  230 CONTINUE
      RETURN
*
  240 CONTINUE
C     PROCEDURE (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
      IF (ISNUM.EQ.1) THEN
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = MS .EQ. M
          ISAME(3) = NS .EQ. N
          ISAME(4) = ALS .EQ. ALPHA
          ISAME(5) = .TRUE.
          IF (M.GT.0 .AND. N.GT.0) ISAME(5) = LSE(AS,A,M,N,LDA)
          ISAME(6) = LDAS .EQ. LDA
          ISAME(7) = .TRUE.
          IF (NL.GT.0 .AND. INCX.NE.0) ISAME(7) = LSE(XS,X,1,NL,
     .        ABS(INCX))
          ISAME(8) = INCXS .EQ. INCX
          ISAME(9) = BLS .EQ. BETA
          ISAME(10) = .TRUE.
          IF (ML.GT.0 .AND. INCY.NE.0) ISAME(10) = LSE(YS,Y,1,ML,
     .        ABS(INCY))
          ISAME(11) = INCYS .EQ. INCY
*
      ELSE IF (ISNUM.EQ.2) THEN
C     COMPARE THE MATRIX IN THE SGBMV DATA STRUCTURE WITH
C     THE SAVED COPY.
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = MS .EQ. M
          ISAME(3) = NS .EQ. N
          ISAME(4) = KLS .EQ. KL
          ISAME(5) = KUS .EQ. KU
          ISAME(6) = ALS .EQ. ALPHA
          ISAME(7) = .TRUE.
          IF (N.GT.0 .AND. M.GT.0) THEN 
              DO 260 J = 1,N
                 DO 250 I = MAX(1,J-KU),MIN(M,J+KL)
                    IF (AS(1+ (I-1)+ (J-1)*LDA).NE.
     .                  A(1+ (KU+I-J)+ (J-1)*LDA)) THEN
                        ISAME(7) = .FALSE.
                        GO TO 270
*
                    END IF
*
  250            CONTINUE
  260         CONTINUE
  270         CONTINUE
          END IF
*
          ISAME(8) = LDAS .EQ. LDA
          ISAME(9) = .TRUE.
          IF (NL.GT.0 .AND. INCX.NE.0) ISAME(9) = LSE(XS,X,1,NL,
     .        ABS(INCX))
          ISAME(10) = INCXS .EQ. INCX
          ISAME(11) = BLS .EQ. BETA
          ISAME(12) = .TRUE.
          IF (ML.GT.0 .AND. INCY.NE.0) ISAME(12) = LSE(YS,Y,1,ML,
     .        ABS(INCY))
          ISAME(13) = INCYS .EQ. INCY
      END IF
*
      GO TO IGO2
*
  280 CONTINUE
C     PROCEDURE (CALL SUBROUTINE)
C     SAVE EVERY DATUM BEFORE THE CALL. 
      ICHS = ICI
      MS = M
      NS = N
      KLS = KL
      KUS = KU
      ALS = ALPHA
      DO 290 I = 1,LDA*N
         AS(I) = A(I)
  290 CONTINUE
      LDAS = LDA
C     SAVE COPY OF THE X AND Y VECTORS. 
      IBX = 1
      IF (INCX.LT.0) IBX = 1 + (1-NL)*INCX
      DO 300 J = 1,NL
         XS(IBX+ (J-1)*INCX) = X(IBX+ (J-1)*INCX) 
  300 CONTINUE
      INCXS = INCX
      BLS = BETA
      IBY = 1
      IF (INCY.LT.0) IBY = 1 + (1-ML)*INCY
      DO 310 I = 1,ML
         YS(IBY+ (I-1)*INCY) = Y(IBY+ (I-1)*INCY) 
  310 CONTINUE
      INCYS = INCY
      IF (ISNUM.EQ.1) THEN
          CALL SGEMV(ICI,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
      ELSE IF (ISNUM.EQ.2) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH SGBMV.
          DO 330 J = 1,N
             DO 320 I = MAX(1,J-KU),MIN(M,J+KL)
                A(1+ (KU+I-J)+ (J-1)*LDA) = AS(1+ (I-1)+ (J-1)*LDA)
  320        CONTINUE
  330     CONTINUE
          CALL SGBMV(ICI,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      END IF
*
      GO TO IGO1
*
  340 CONTINUE
C     PROCEDURE (DEFINE A SET OF PROBLEM DATA)
C     DO NOTHING IF BOTH DIMENSIONS ARE NOT POSITIVE.
      IF (M.LE.0 .OR. N.LE.0) GO TO IGO3
      TRANSL = ZERO 
      CALL SMAKE(A,M,N,LDA,RESET,TRANSL)
C     TRIM AWAY ELEMENTS OUTSIDE THE BANDWIDTH FOR SGBMV.
      IF (ISNUM.EQ.2) THEN
          DO 360 J = 1,N
             DO 350 I = 1,M
                T = A(1+ (I-1)+ (J-1)*LDA)
                IF (J.GT.I .AND. J-I.GT.KU) T = ZERO
                IF (I.GT.J .AND. I-J.GT.KL) T = ZERO
                A(1+ (I-1)+ (J-1)*LDA) = T
  350        CONTINUE
  360     CONTINUE
      END IF
*
      TRANSL = 500.E0
      RESET = .FALSE.
      CALL SMAKE(X,1,NL,MAX(1,ABS(INCX)),RESET,TRANSL)
      IF (NL.GT.1 .AND. INCX.EQ.1) X(NL/2) = ZERO 
      TRANSL = ZERO 
      CALL SMAKE(Y,1,ML,MAX(1,ABS(INCY)),RESET,TRANSL)
      GO TO IGO3
*
  370 CONTINUE
C     PROCEDURE (COMPUTE A CORRECT RESULT)
C     COMPUTE THE CONDITION NUMBER TO USE AS GAUGE FOR ACCURATE RESULTS.
C     THIS IS RETURNED IN XT(*).
C     COMPUTE THE APPROXIMATE CORRECT RESULT.
C     THIS IS RETURNED IN YT(*).
      IF (INCY.LT.0) THEN
          IBY = (1-ML)*INCY + 1
*
      ELSE
          IBY = 1
      END IF
*
      DO 390 I = 1,ML
         YT(I) = BETA*YS(IBY+ (I-1)*INCY)
         XT(I) = YS(IBY+ (I-1)*INCY)**2 
         IF (INCX.LT.0) THEN
             IBX = (1-NL)*INCX + 1
*
         ELSE
             IBX = 1
         END IF
*
         DO 380 J = 1,NL
            YT(I) = YT(I) + AS(1+ (I-1)*INCRA+ (J-1)*INCCA)*ALPHA*
     .              XS(IBX+ (J-1)*INCX) 
            XT(I) = XT(I) + AS(1+ (I-1)*INCRA+ (J-1)*INCCA)**2
  380    CONTINUE
         XT(I) = SQRT(XT(I))
  390 CONTINUE
      XN = BETA**2
      DO 400 J = 1,NL
         XN = XN + XS(IBX+ (J-1)*INCX)**2
  400 CONTINUE
      XN = SQRT(XN) 
C     COMPUTE THE GAUGES FOR THE RESULTS.
      DO 410 I = 1,ML
         XT(I) = XT(I)*XN
  410 CONTINUE
C     COMPUTE THE DIFFERENCES. THEY SHOULD BE SMALL FOR CORRECT RESULTS.
      DO 420 I = 1,ML
         YT(I) = YT(I) - Y(IBY+ (I-1)*INCY)
  420 CONTINUE
C     COMPUTE THE GRADE OF THIS RESULT. 
      IGR = 0
      T = ONE
  430 CONTINUE
C     THIS TEST ALLOWS UP TO A LOSS OF FULL PRECISION BEFORE QUITTING.
      IF (IGR.GE.IG) GO TO 460
      DO 440 I = 1,ML
         IF (SDIFF(T*ABS(YT(I))+XT(I),XT(I)).EQ.ZERO) GO TO 440
         T = T*HALF 
         IGR = IGR + 1
         GO TO 430
*
  440 CONTINUE
C     IF THE LOOP COMPLETES, ALL VALUES ARE 'SMALL.'  THE VALUE IGR/IG
C     IS THE GRADE ASSIGNED.  THE VALUE OF IGR IS MAXIMIZED OVER ALL THE
C     PROBLEMS.
  450 CONTINUE
      AVIGR = MAX(AVIGR,REAL(IGR))
      GO TO IGO4
*
  460 CONTINUE
      FATAL = .TRUE.
      GO TO 450
*
*     LAST EXECUTABLE LINE OF SCHCK1
 9001 FORMAT (' IN SUBPROGRAM ',A,/,' TESTS ACTIVE WITH OPTION = ',A,/,
     .  '  M =',I4,', N = ',I4,/,' INCX = ',I2,', INCY = ',I2,/,' KL =',
     .  I4,', KU =',I4)
 9011 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WITH INVALID INPUT.',/,' OPTION = ',A,', M =',I4,
     .  ', N = ',I4,/,' INCX = ',I2,', INCY = ',I2,/,' KL =',I4,
     .  ', KU =',I4)
 9021 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WHILE COMPUTING',/,' OPTION = ',A,', M =',I4,
     .  ', N = ',I4,/,' INCX = ',I2,', INCY = ',I2,/,' KL =',I4,
     .  ', KU =',I4)
 9031 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'RECEIVED A LOSS GRADE OF ', 
     .  F5.2,' OUT OF ',I3)
 9041 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'FAILED.')
 9051 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'NOT TESTED.')
      END 
      SUBROUTINE SCHCK2(ISNUM,SNAME,IG,DOPE,NUNIT,AVIGR,FATAL)
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     TEST SSYMV, 03, SSBMV, 04, AND SSPMV, 05.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C     THIS PROGRAM HAS TWO PARTS.  THE FIRST PART CHECKS TO SEE
C     IF ANY DATA GETS CHANGED WHEN NONE SHOULD.  FOR EXAMPLE WHEN
C     USING AN INVALID OPTION OR NONPOSITVE PROBLEM DIMENSIONS.
C     THE SECOND PART CHECKS THAT THE RESULTS ARE REASONABLY ACCURATE.
C     DIMENSION AND PROBLEM SIZE DATA.. 
      INTEGER INC(04),IDIM(06),NUNIT(2) 
      REAL ALF(04),BET(04)
      LOGICAL ISAME(13),LSE,FATAL,SAME,NCHNG,RESET
      CHARACTER *128 DOPE(2)
      CHARACTER *6 SNAME
      CHARACTER *3 ICH
      CHARACTER *1 ICHS,ICI
      INTEGER LA,LV 
      PARAMETER (LA=4096,LV=4096,LMN=2048)
      REAL ALPHA,ALS,BETA,BLS,T,TRANSL,XN
      REAL A(LA),AS(LA),X(LV),XS(LV)
      REAL Y(LV),YS(LV),YT(LMN),XT(LMN) 
      PARAMETER (ZERO=0.E0,HALF=.5E0,ONE=1.E0)
      COMMON /ARRAYS/AR,AS,X,XS,Y,YS,YT,XT
      EXTERNAL SDIFF
      INTEGER LIJU_ZERO(3)
      DATA LIJU_ZERO/3*0/
*
      DATA ALF/-1.E0,2.E0,.3E0,1.E0/
      DATA BET/-1.E0,0.E0,.3E0,1.E0/
      DATA INC/-2,-1,1,2/
      DATA IDIM/1,2,4,8,64,0/ 
      DATA ICH/'LU/'/
      FATAL = .FALSE.
C     CHECK SYMMETRIC MATRIX-VECTOR PRODUCT, Y = ALPHA*A*X+BETA*Y, 3-5.
      IF (ISNUM.LT.0) GO TO 200
      NC = 0
      RESET = .TRUE.
      AVIGR = ZERO
      IX = 0
   10 IX = IX + 1
      IF (IX.GT.4) GO TO 180
      INCX = INC(IX)
      ALPHA = ALF(IX)
      IY = 0
   20 IY = IY + 1
      IF (IY.GT.4) GO TO 170
      INCY = INC(IY)
      BETA = BET(IY)
      NN = 0
   30 NN = NN + 1
      IF (NN.GT.6) GO TO 160
      N = IDIM(NN)
      IC = 0
   40 IC = IC + 1
      IF (IC.GT.3) GO TO 150
      IF (FATAL) GO TO 190
      ICI = ICH(IC:IC)
C     DEFINE DEFAULT VALUE OF K SO PRINTING IS OK.
      K = MAX(0,N-1)
C     DEFINE THE NUMBER OF ARGUMENTS AND THE Y ARGUMENT NUMBER.
      LDA = MAX(N,1)
      IF (ISNUM.EQ.3) THEN
          NARGS = 10
          IYARG = 09
*
      ELSE IF (ISNUM.EQ.4) THEN
          NARGS = 11
          IYARG = 10
C     DEFINE BANDWIDTH OF MATRIX FOR TEST OF SSBMV.
          K = INT(SQRT(REAL(N))+HALF) - 1
*
      ELSE IF (ISNUM.EQ.5) THEN
          NARGS = 9 
          IYARG = 8 
      END IF
C     DO (PREPARE NOTES FOR THIS TEST)
C
C     PRINT A MESSAGE THAT GIVES DEBUGGING INFORMATION.  THIS
C     MESSAGE SAYS..
C     IN SUBPROGRAM XXXXX TESTS WERE ACTIVE WITH
C     OPTION = 'A'
C     N = IIII,
C     INCX = IIII,  INCY = IIII,
C     K = IIII.
C     THE MAIN IDEA HERE IS THAT A SERIOUS BUG THAT OCCURS DURING THE 
C     TESTING OF THESE SUBPROGRAMS MAY LOSE PROGRAM CONTROL.  THIS
C     AUXILLIARY FILE CONTAINS THE DIMENSIONS THAT RESULTED IN THE LOSS
C     OF CONTROL.  HENCE IT PROVIDES THE IMPLEMENTOR WITH MORE COMPLETE
C     INFORMATION ABOUT WHERE TO START TRACKING DOWN THE BUG.
      IF (NUNIT(1).GT.0) THEN 
C     IF UNIT IS NOT AVAILABLE WITH 'NEW' STATUS, OPEN WITH 
C     'OLD' AND THEN DELETE IT.
          ISTAT = 1 
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.1) GO TO 50
C     GET RID OF ANY OLD FILE.
          CLOSE (UNIT=NUNIT(1),STATUS='DELETE',ERR=50)
   50     CONTINUE
          ISTAT = 2 
C     CREATE A NEW FILE FOR THE NEXT TEST.
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.0) GO TO 70
   60     CONTINUE
          NMESS = 7 
C     DO (PRINT A MESSAGE)
C     PRINT AN ERROR MESSAGE ABOUT OPENING THE NAME FILE.
          CALL SMESSG(LIJU_ZERO,1,NMESS)
          FATAL = .TRUE.
          GO TO 190 
*
   70     CONTINUE
          WRITE (NUNIT(1),9001) SNAME,ICI,N,INCX,INCY,K
C     CLOSE THE FILE SO USEFUL STATUS INFORMATION IS SEALED.
          CLOSE (UNIT=NUNIT(1))
      END IF
C     DO (DEFINE A SET OF PROBLEM DATA) 
      ASSIGN 80 TO IGO3
      GO TO 370
*
   80 CONTINUE
C     DO (CALL SUBROUTINE)
      ASSIGN 90 TO IGO1
      GO TO 290
*
   90 CONTINUE
      IF (N.LE.0 .OR. ICHAR(ICI).EQ.ICHAR('/')) THEN
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 100 TO IGO2
          GO TO 220 
*
  100     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 110 I = 1,NARGS
             SAME = SAME .AND. ISAME(I) 
             IF ( .NOT. ISAME(I)) THEN
                 WRITE (NUNIT(2),9011) SNAME,I,ICI,N,INCX,INCY,K
             END IF 
*
  110     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 190
*
          END IF
*
      ELSE
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 120 TO IGO2
          GO TO 220 
*
  120     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 130 I = 1,NARGS
             NCHNG = (I.EQ.IYARG .OR. ISAME(I))
             SAME = SAME .AND. NCHNG
             IF ( .NOT. NCHNG) THEN
                 WRITE (NUNIT(2),9021) SNAME,I,ICI,N,INCX,INCY,K
             END IF 
*
  130     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 190
*
          END IF
*
          NC = NC + 1
C     DO (COMPUTE A CORRECT RESULT)
          ASSIGN 140 TO IGO4
          GO TO 420 
*
  140     CONTINUE
C     IF GOT REALLY BAD ANSWER, PRINT NOTE AND GO BACK.
          IF (FATAL) GO TO 180
*
      END IF
*
      GO TO 40
*
  150 CONTINUE
      GO TO 30
*
  160 CONTINUE
      GO TO 20
*
  170 CONTINUE
      GO TO 10
*
  180 CONTINUE
C     REPORT ON ACCURACY OF DATA.
      WRITE (NUNIT(2),9031) ISNUM,SNAME,AVIGR,IG
      GO TO 210
*
  190 CONTINUE
      WRITE (NUNIT(2),9041) ISNUM,SNAME 
      GO TO 210
*
  200 CONTINUE
      WRITE (NUNIT(2),9051) - ISNUM,SNAME
  210 CONTINUE
      RETURN
*
  220 CONTINUE
C     PROCEDURE (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
      IF (ISNUM.EQ.3) THEN
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = ALS .EQ. ALPHA
          ISAME(4) = .TRUE.
          IF (N.GT.0) ISAME(4) = LSE(AS,A,N,N,LDA)
          ISAME(5) = LDAS .EQ. LDA
          ISAME(6) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(6) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(7) = INCXS .EQ. INCX
          ISAME(8) = BLS .EQ. BETA
          ISAME(9) = .TRUE.
          IF (N.GT.0 .AND. INCY.NE.0) ISAME(9) = LSE(YS,Y,1,N,ABS(INCY))
          ISAME(10) = INCYS .EQ. INCY
*
      ELSE IF (ISNUM.EQ.4) THEN
C     COMPARE THE MATRIX IN THE SSBMV AND SSPMV DATA STRUCTURES WITH
C     THE SAVED COPY.
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = KS .EQ. K
          ISAME(4) = ALS .EQ. ALPHA
          ISAME(5) = .TRUE.
C     TEST THE MATRIX IN THE DATA STRUCTURE USED WITH SSBMV.
          IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
              KOFF = K
*
          ELSE
              KOFF = 0
          END IF
*
          IF (N.GT.0) THEN
              DO 240 J = 1,N
                 DO 230 I = MAX(1,J-K),MIN(N,J+K) 
                    IF (AS(1+ (I-1)+ (J-1)*LDA).NE.
     .                  A(1+ (KOFF+I-J)+ (J-1)*LDA)) THEN
                        ISAME(5) = .FALSE.
                        GO TO 250
*
                    END IF
*
  230            CONTINUE
  240         CONTINUE
  250         CONTINUE
          END IF
*
          ISAME(6) = LDAS .EQ. LDA
          ISAME(7) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(7) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(8) = INCXS .EQ. INCX
          ISAME(9) = BLS .EQ. BETA
          ISAME(10) = .TRUE.
          IF (N.GT.0 .AND. INCY.NE.0) ISAME(10) = LSE(YS,Y,1,N,
     .        ABS(INCY))
          ISAME(11) = INCYS .EQ. INCY
*
      ELSE IF (ISNUM.EQ.5) THEN
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = ALS .EQ. ALPHA
          ISAME(4) = .TRUE.
C     TEST THE MATRIX USING THE DATA STRUCTURE USED WITH SSPMV.
          IOFF = 0
          DO 270 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 260 I = ISTRT,IEND
                IOFF = IOFF + 1
                IF (A(IOFF).NE.AS(1+ (I-1)+ (J-1)*LDA)) THEN
                    ISAME(4) = .FALSE.
                    GO TO 280 
*
                END IF
*
  260        CONTINUE
*
  270     CONTINUE
  280     CONTINUE
          ISAME(5) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(5) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(6) = INCXS .EQ. INCX
          ISAME(7) = BLS .EQ. BETA
          ISAME(8) = .TRUE.
          IF (N.GT.0 .AND. INCY.NE.0) ISAME(8) = LSE(YS,Y,1,N,ABS(INCY))
          ISAME(9) = INCYS .EQ. INCY
      END IF
*
      GO TO IGO2
*
  290 CONTINUE
C     PROCEDURE (CALL SUBROUTINE)
C     SAVE EVERY DATUM BEFORE THE CALL. 
      ICHS = ICI
      NS = N
      KS = K
      ALS = ALPHA
      DO 300 I = 1,N*N
         AS(I) = A(I)
  300 CONTINUE
      LDAS = LDA
C     SAVE COPY OF THE X AND Y VECTORS. 
      IBX = 1
      IF (INCX.LT.0) IBX = 1 + (1-N)*INCX
      DO 310 J = 1,N
         XS(IBX+ (J-1)*INCX) = X(IBX+ (J-1)*INCX) 
  310 CONTINUE
      INCXS = INCX
      BLS = BETA
      IBY = 1
      IF (INCY.LT.0) IBY = 1 + (1-N)*INCY
      DO 320 I = 1,N
         YS(IBY+ (I-1)*INCY) = Y(IBY+ (I-1)*INCY) 
  320 CONTINUE
      INCYS = INCY
      IF (ISNUM.EQ.3) THEN
          CALL SSYMV(ICI,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
      ELSE IF (ISNUM.EQ.4) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH SSBMV.
          IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
              KOFF = K
*
          ELSE
              KOFF = 0
          END IF
*
          DO 340 J = 1,N
             DO 330 I = MAX(1,J-K),MIN(N,J+K)
                A(1+ (KOFF+I-J)+ (J-1)*LDA) = AS(1+ (I-1)+ (J-1)*LDA) 
  330        CONTINUE
  340     CONTINUE
          CALL SSBMV(ICI,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
      ELSE IF (ISNUM.EQ.5) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH SSPMV.
          IOFF = 0
          DO 360 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 350 I = ISTRT,IEND
                IOFF = IOFF + 1
                A(IOFF) = AS(1+ (I-1)+ (J-1)*LDA) 
  350        CONTINUE
*
  360     CONTINUE
          CALL SSPMV(ICI,N,ALPHA,A,X,INCX,BETA,Y,INCY)
      END IF
*
      GO TO IGO1
*
  370 CONTINUE
C     PROCEDURE (DEFINE A SET OF PROBLEM DATA)
C     DO NOTHING IF DIMENSIONS ARE NOT POSITIVE.
      IF (N.LE.0) GO TO IGO3
      TRANSL = ZERO 
      CALL SMAKE(A,N,N,LDA,RESET,TRANSL)
C     MAKE THE DATA MATRIX SYMMETRIC.
      DO 390 I = 1,N
         DO 380 J = I,N
            T = (A(1+ (I-1)+ (J-1)*LDA)+A(1+ (J-1)+ (I-1)*LDA))*HALF
            A(1+ (I-1)+ (J-1)*LDA) = T
            A(1+ (J-1)+ (I-1)*LDA) = T
  380    CONTINUE
  390 CONTINUE
C     TRIM AWAY ELEMENTS OUTSIDE THE BANDWIDTH FOR SSBMV.
      IF (ISNUM.EQ.4) THEN
          DO 410 J = 1,N
             DO 400 I = 1,N
                T = A(1+ (I-1)+ (J-1)*LDA)
                IF (J.GT.I .AND. J-I.GT.K) T = ZERO
                IF (I.GT.J .AND. I-J.GT.K) T = ZERO
                A(1+ (I-1)+ (J-1)*LDA) = T
  400        CONTINUE
  410     CONTINUE
      END IF
*
      TRANSL = 500.E0
      RESET = .FALSE.
      CALL SMAKE(X,1,N,MAX(1,ABS(INCX)),RESET,TRANSL)
      IF (N.GT.1 .AND. INCX.EQ.1) X(N/2) = ZERO
      TRANSL = ZERO 
      CALL SMAKE(Y,1,N,MAX(1,ABS(INCY)),RESET,TRANSL)
      GO TO IGO3
*
  420 CONTINUE
C     PROCEDURE (COMPUTE A CORRECT RESULT)
C     COMPUTE THE CONDITION NUMBER TO USE AS GAUGE FOR ACCURATE RESULTS.
C     THIS IS RETURNED IN XT(*).
C     COMPUTE THE APPROXIMATE CORRECT RESULT.
C     THIS IS RETURNED IN YT(*).
      IF (INCY.LT.0) THEN
          IBY = (1-N)*INCY + 1
*
      ELSE
          IBY = 1
      END IF
*
      DO 440 I = 1,N
         YT(I) = BETA*YS(IBY+ (I-1)*INCY)
         XT(I) = YS(IBY+ (I-1)*INCY)**2 
         IF (INCX.LT.0) THEN
             IBX = (1-N)*INCX + 1
*
         ELSE
             IBX = 1
         END IF
*
         DO 430 J = 1,N
            YT(I) = YT(I) + AS(1+ (I-1)+ (J-1)*LDA)*ALPHA*
     .              XS(IBX+ (J-1)*INCX) 
            XT(I) = XT(I) + AS(1+ (I-1)+ (J-1)*LDA)**2
  430    CONTINUE
         XT(I) = SQRT(XT(I))
  440 CONTINUE
      XN = BETA**2
      DO 450 J = 1,N
         XN = XN + XS(IBX+ (J-1)*INCX)**2
  450 CONTINUE
      XN = SQRT(XN) 
C     COMPUTE THE GAUGES FOR THE RESULTS.
      DO 460 I = 1,N
         XT(I) = XT(I)*XN
  460 CONTINUE
C     COMPUTE THE DIFFERENCES. THEY SHOULD BE SMALL FOR CORRECT RESULTS.
      DO 470 I = 1,N
         YT(I) = YT(I) - Y(IBY+ (I-1)*INCY)
  470 CONTINUE
C     COMPUTE THE GRADE OF THIS RESULT. 
      IGR = 0
      T = ONE
  480 CONTINUE
C     THIS TEST ALLOWS UP TO A LOSS OF FULL PRECISION BEFORE QUITTING.
      IF (IGR.GT.IG) GO TO 510
      DO 490 I = 1,N
         IF (SDIFF(T*ABS(YT(I))+XT(I),XT(I)).EQ.ZERO) GO TO 490
         T = T*HALF 
         IGR = IGR + 1
         GO TO 480
*
  490 CONTINUE
C     IF THE LOOP COMPLETES, ALL VALUES ARE 'SMALL.'  THE VALUE IGR/IG
C     IS THE GRADE ASSIGNED.  THE VALUE OF IGR IS MAXIMIZED OVER ALL THE
C     PROBLEMS.
  500 CONTINUE
      AVIGR = MAX(AVIGR,REAL(IGR))
      GO TO IGO4
*
  510 CONTINUE
      FATAL = .TRUE.
      GO TO 500
*
*     LAST EXECUTABLE LINE OF SCHCK2
 9001 FORMAT (' IN SUBPROGRAM ',A,/,' TESTS ACTIVE WITH OPTION = ',A,/,
     .  ' N = ',I4,/,' INCX = ',I2,', INCY = ',I2,/,' K =',I4)
 9011 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WITH INVALID INPUT.',/,' OPTION = ',A,/,' N = ',
     .  I4,/,' INCX = ',I2,', INCY = ',I2,/,' K = ',I4)
 9021 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WHILE COMPUTING',/,' OPTION = ',A,/,' N = ',I4,/,
     .  ' INCX = ',I2,', INCY = ',I2,/,' K = ',I4)
 9031 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'RECEIVED A LOSS GRADE OF ', 
     .  F5.2,' OUT OF ',I3)
 9041 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'FAILED.')
 9051 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'NOT TESTED.')
      END 
      SUBROUTINE SCHCK3(ISNUM,SNAME,IG,DOPE,NUNIT,AVIGR,FATAL)
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     TEST STRMV, 06, STBMV, 07, STPMV, 08,
C     STRSV, 09, STBSV, 10, AND STPSV, 11.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R.  J.  HANSON, SANDIA NATIONAL LABS.
C     THIS PROGRAM HAS TWO PARTS.  THE FIRST PART CHECKS TO SEE
C     IF ANY DATA GETS CHANGED WHEN NONE SHOULD.  FOR EXAMPLE WHEN
C     USING AN INVALID OPTION OR NONPOSITVE PROBLEM DIMENSIONS.
C     THE SECOND PART CHECKS THAT THE RESULTS ARE REASONABLY ACCURATE.
C     DIMENSION AND PROBLEM SIZE DATA.. 
      INTEGER INC(04),IDIM(06),NUNIT(2) 
      LOGICAL ISAME(13),LSE,FATAL,SAME,NCHNG,RESET
      CHARACTER *128 DOPE(2)
      CHARACTER *6 SNAME
      CHARACTER *3 ICHI,ICHJ,ICHK
      CHARACTER *1 ICIU,ICIT,ICID
      CHARACTER *1 ICIUS,ICITS,ICIDS
      INTEGER LA,LV 
      PARAMETER (LA=4096,LV=4096,LMN=2048)
      REAL A(LA),AS(LA),X(LV),XS(LV)
      REAL Y(LV),YS(LV),YT(LMN),XT(LMN) 
      PARAMETER (ZERO=0.E0,HALF=.5E0,ONE=1.E0)
      COMMON /ARRAYS/AR,AS,X,XS,Y,YS,XT,YT
      EXTERNAL SDIFF
      INTEGER LIJU_ZERO(3)
      DATA LIJU_ZERO/3*0/
*
      DATA INC/-2,-1,1,2/
      DATA IDIM/1,2,4,8,64,0/ 
      DATA ICHI/'LU/'/,ICHJ/'NT/'/,ICHK/'NU/'/
      FATAL = .FALSE.
C     CHECK TRIANGULAR MATRIX-VECTOR PRODUCT, X = A*X, 6-8, 
C     AND TRIANGULAR SOLVERS, 9-11.
      IF (ISNUM.LT.0) GO TO 180
      NC = 0
      RESET = .TRUE.
      AVIGR = ZERO
      IX = 0
   10 IX = IX + 1
      IF (IX.GT.4) GO TO 160
      INCX = INC(IX)
      NN = 0
   20 NN = NN + 1
      IF (NN.GT.6) GO TO 150
      N = IDIM(NN)
      IC = 0
   30 IC = IC + 1
      IF (IC.GT.3) GO TO 140
      IF (FATAL) GO TO 170
      ICIU = ICHI(IC:IC)
      ICIT = ICHJ(IC:IC)
      ICID = ICHK(IC:IC)
C     DEFINE DEFAULT VALUE OF K SO PRINTING IS OK.
      K = MAX(0,N-1)
C     DEFINE THE NUMBER OF ARGUMENTS AND THE X ARGUMENT NUMBER.
      LDA = MAX(N,1)
      IF (ICHAR(ICIT).EQ.ICHAR('T')) THEN
          INCRA = LDA
          INCCA = 1 
*
      ELSE
          INCRA = 1 
          INCCA = LDA
      END IF
*
      IF (ISNUM.EQ.6 .OR. ISNUM.EQ.9) THEN
          NARGS = 08
          IXARG = 07
*
      ELSE IF (ISNUM.EQ.7 .OR. ISNUM.EQ.10) THEN
          NARGS = 09
          IXARG = 08
C     DEFINE BANDWIDTH OF MATRIX FOR TEST OF STBMV.
          K = INT(SQRT(REAL(N))+HALF) - 1
*
      ELSE IF (ISNUM.EQ.8 .OR. ISNUM.EQ.11) THEN
          NARGS = 07
          IXARG = 06
      END IF
C     DO (PREPARE NOTES FOR THIS TEST)
C
C     PRINT A MESSAGE THAT GIVES DEBUGGING INFORMATION.  THIS
C     MESSAGE SAYS..
C     IN SUBPROGRAM XXXXX TESTS WERE ACTIVE WITH
C     OPTIONS = 'A' 'A' 'A'
C     N = IIII,
C     INCX = IIII C K = IIII. 
C     THE MAIN IDEA HERE IS THAT A SERIOUS BUG THAT OCCURS DURING THE 
C     TESTING OF THESE SUBPROGRAMS MAY LOSE PROGRAM CONTROL.  THIS
C     AUXILLIARY FILE CONTAINS THE DIMENSIONS THAT RESULTED IN THE LOSS
C     OF CONTROL. HENCE IT PROVIDES THE IMPLEMENTOR WITH MORE COMPLETE
C     INFORMATION ABOUT WHERE TO START TRACKING DOWN THE BUG.
      IF (NUNIT(1).GT.0) THEN 
C     IF UNIT IS NOT AVAILABLE WITH 'NEW' STATUS, OPEN WITH 
C    'OLD' AND THEN DELETE IT.
          ISTAT = 1 
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.1) GO TO 40
C    GET RID OF ANY OLD FILE. 
          CLOSE (UNIT=NUNIT(1),STATUS='DELETE',ERR=40)
   40     CONTINUE
          ISTAT = 2 
C    CREATE A NEW FILE FOR THE NEXT TEST.
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.0) GO TO 60
   50     CONTINUE
          NMESS = 7 
C     DO (PRINT A MESSAGE)
C     PRINT AN ERROR MESSAGE ABOUT OPENING THE NAME FILE.
          CALL SMESSG(LIJU_ZERO,1,NMESS)
          FATAL = .TRUE.
          GO TO 170 
*
   60     CONTINUE
          WRITE (NUNIT(1),9001) SNAME,ICIU,ICIT,ICID,N,INCX,K
C     CLOSE THE FILE SO USEFUL STATUS INFORMATION IS SEALED.
          CLOSE (UNIT=NUNIT(1))
      END IF
C     DO (DEFINE A SET OF PROBLEM DATA) 
      ASSIGN 70 TO IGO3
      GO TO 330
*
   70 CONTINUE
C     DO (CALL SUBROUTINE)
      ASSIGN 80 TO IGO1
      GO TO 260
*
   80 CONTINUE
      IF (N.LE.0 .OR. ICHAR(ICIU).EQ.ICHAR('/') .OR. ICHAR(ICIT).EQ.
     .    ICHAR('/') .OR. ICHAR(ICID).EQ.ICHAR('/')) THEN
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 90 TO IGO2
          GO TO 200 
*
   90     CONTINUE
C      IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 100 I = 1,NARGS
             SAME = SAME .AND. ISAME(I) 
             IF ( .NOT. ISAME(I)) THEN
                 WRITE (NUNIT(2),9011) SNAME,I,ICIU,ICIT,ICID,N,INCX,K
             END IF 
*
  100     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 170
*
          END IF
*
      ELSE
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 110 TO IGO2
          GO TO 200 
*
  110     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 120 I = 1,NARGS
             NCHNG = (I.EQ.IXARG .OR. ISAME(I))
             SAME = SAME .AND. NCHNG
             IF ( .NOT. NCHNG) THEN
                 WRITE (NUNIT(2),9021) SNAME,I,ICIU,ICIT,ICID,N,INCX,K
             END IF 
*
  120     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 170
*
          END IF
*
          NC = NC + 1
C     DO (COMPUTE A CORRECT RESULT)
          ASSIGN 130 TO IGO4
          GO TO 380 
*
  130     CONTINUE
C     IF GOT REALLY BAD ANSWER, PRINT NOTE AND GO BACK.
          IF (FATAL) GO TO 160
*
      END IF
*
      GO TO 30
*
  140 CONTINUE
      GO TO 20
*
  150 CONTINUE
      GO TO 10
*
  160 CONTINUE
C     REPORT ON ACCURACY OF DATA.
      WRITE (NUNIT(2),9031) ISNUM,SNAME,AVIGR,IG
      GO TO 190
*
  170 CONTINUE
      WRITE (NUNIT(2),9041) ISNUM,SNAME 
      GO TO 190
*
  180 CONTINUE
      WRITE (NUNIT(2),9051) - ISNUM,SNAME
  190 CONTINUE
      RETURN
*
  200 CONTINUE
C     PROCEDURE (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
      ISAME(1) = ICHAR(ICIU) .EQ. ICHAR(ICIUS)
      ISAME(2) = ICHAR(ICIT) .EQ. ICHAR(ICITS)
      ISAME(3) = ICHAR(ICID) .EQ. ICHAR(ICIDS)
      ISAME(4) = NS .EQ. N
      IF (ISNUM.EQ.6 .OR. ISNUM.EQ.9) THEN
          ISAME(5) = .TRUE.
          IF (N.GT.0) ISAME(5) = LSE(AS,A,N,N,LDA)
          ISAME(6) = LDAS .EQ. LDA
          ISAME(7) = .TRUE.
          IF (N.GT.0) ISAME(7) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(8) = INCXS .EQ. INCX
*
      ELSE IF (ISNUM.EQ.7 .OR. ISNUM.EQ.10) THEN
C     COMPARE THE MATRIX IN THE STBMV AND STPMV DATA STRUCTURES WITH
C     THE SAVED COPY.
          ISAME(5) = KS .EQ. K
          ISAME(6) = .TRUE.
          IF (N.GT.0) THEN
              DO 220 J = 1,N
                 IF (ICHAR(ICIU).EQ.ICHAR('U')) THEN
                     ISTRT = MAX(1,J-K) 
                     IEND = J 
*
                 ELSE
                     ISTRT = J
                     IEND = MIN(N,J+K)
                 END IF
*
                 DO 210 I = ISTRT,IEND
                    IF (AS(1+ (I-1)+ (J-1)*LDA).NE.
     .                  A(1+ (KOFF+I-J)+ (J-1)*LDA)) THEN
                        ISAME(6) = .FALSE.
                        GO TO 230
*
                    END IF
*
  210            CONTINUE
  220         CONTINUE
  230         CONTINUE
          END IF
*
          ISAME(7) = LDAS .EQ. LDA
          ISAME(8) = .TRUE.
          IF (N.GT.0) ISAME(8) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(9) = INCXS .EQ. INCX
*
      ELSE IF (ISNUM.EQ.8 .OR. ISNUM.EQ.11) THEN
          ISAME(5) = .TRUE.
          IOFF = 0
          DO 250 J = 1,N
             IF (ICHAR(ICIU).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 240 I = ISTRT,IEND
                IOFF = IOFF + 1
                IF (A(IOFF).NE.AS(1+ (I-1)+ (J-1)*
     .              LDA)) ISAME(5) = .FALSE.
  240        CONTINUE
*
  250     CONTINUE
          ISAME(6) = .TRUE.
          IF (N.GT.0) ISAME(6) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(7) = INCXS .EQ. INCX
      END IF
*
      GO TO IGO2
*
  260 CONTINUE
C     PROCEDURE (CALL SUBROUTINE)
C     SAVE EVERY DATUM BEFORE THE CALL. 
      ICIUS = ICIU
      ICITS = ICIT
      ICIDS = ICID
      NS = N
      KS = K
      DO 270 I = 1,N*N
         AS(I) = A(I)
  270 CONTINUE
      LDAS = LDA
C     SAVE COPY OF THE X VECTOR.
      IBX = 1
      IF (INCX.LT.0) IBX = 1 + (1-N)*INCX
      DO 280 J = 1,N
         XS(IBX+ (J-1)*INCX) = X(IBX+ (J-1)*INCX) 
  280 CONTINUE
      INCXS = INCX
      IF (ISNUM.EQ.6) THEN
          CALL STRMV(ICIU,ICIT,ICID,N,A,LDA,X,INCX)
*
      ELSE IF (ISNUM.EQ.9) THEN
          CALL STRSV(ICIU,ICIT,ICID,N,A,LDA,X,INCX)
*
      ELSE IF (ISNUM.EQ.7 .OR. ISNUM.EQ.10) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH STBMV.
          IF (ICHAR(ICIU).EQ.ICHAR('U')) THEN
              KOFF = K
*
          ELSE
              KOFF = 0
          END IF
*
          DO 300 J = 1,N
             IF (ICHAR(ICIU).EQ.ICHAR('U')) THEN
                 ISTRT = MAX(1,J-K)
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = MIN(N,J+K)
             END IF 
*
             DO 290 I = ISTRT,IEND
                A(1+ (KOFF+I-J)+ (J-1)*LDA) = AS(1+ (I-1)+ (J-1)*LDA) 
  290        CONTINUE
  300     CONTINUE
          IF (ISNUM.EQ.7) CALL STBMV(ICIU,ICIT,ICID,N,K,A,LDA,X,INCX) 
          IF (ISNUM.EQ.10) CALL STBSV(ICIU,ICIT,ICID,N,K,A,LDA,X,INCX)
*
      ELSE IF (ISNUM.EQ.8 .OR. ISNUM.EQ.11) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH STPMV.
          IOFF = 0
          DO 320 J = 1,N
             IF (ICHAR(ICIU).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 310 I = ISTRT,IEND
                IOFF = IOFF + 1
                A(IOFF) = AS(1+ (I-1)+ (J-1)*LDA) 
  310        CONTINUE
*
  320     CONTINUE
          IF (ISNUM.EQ.8) CALL STPMV(ICIU,ICIT,ICID,N,A,X,INCX)
          IF (ISNUM.EQ.11) CALL STPSV(ICIU,ICIT,ICID,N,A,X,INCX)
      END IF
*
      GO TO IGO1
*
  330 CONTINUE
C     PROCEDURE (DEFINE A SET OF PROBLEM DATA)
C     DO NOTHING IF DIMENSIONS ARE NOT POSITIVE.
      IF (N.LE.0) GO TO IGO3
      TRANSL = ZERO 
      CALL SMAKE(A,N,N,LDA,RESET,TRANSL)
C     MAKE THE DATA MATRIX TRIANGULAR.
      DO 350 I = 1,N
         DO 340 J = 1,N
            T = A(1+INCRA* (I-1)+ (J-1)*INCCA)
            S = A(1+INCRA* (J-1)+ (I-1)*INCCA)
C     SCALE TERMS SO THAT UNIT MATRICES ARE WELL-CONDITIONED.
            S = S/1000.E0
            T = T/1000.E0
            IF (ICHAR(ICIU).EQ.ICHAR('L') .AND. I.LT.J) T = ZERO
            IF (ICHAR(ICIU).EQ.ICHAR('U') .AND. I.GT.J) S = ZERO
            IF (ICHAR(ICID).EQ.ICHAR('U') .AND. I.EQ.J) THEN
                S = ONE
                T = ONE
            END IF
*
            A(1+INCRA* (I-1)+ (J-1)*INCCA) = T
            A(1+INCRA* (J-1)+ (I-1)*INCCA) = S
  340    CONTINUE
  350 CONTINUE
C     TRIM AWAY ELEMENTS OUTSIDE THE BANDWIDTH FOR STBMV.
      IF (ISNUM.EQ.7 .OR. ISNUM.EQ.10) THEN
          DO 370 I = 1,N
             DO 360 J = 1,N
                T = A(1+INCRA* (I-1)+ (J-1)*INCCA)
                IF (J.GT.I .AND. J-I.GT.K) T = ZERO
                IF (I.GT.J .AND. I-J.GT.K) T = ZERO
                A(1+INCRA* (I-1)+ (J-1)*INCCA) = T
  360        CONTINUE
  370     CONTINUE
      END IF
*
      TRANSL = 500.E0
      RESET = .FALSE.
      CALL SMAKE(X,1,N,MAX(1,ABS(INCX)),RESET,TRANSL)
      IF (N.GT.1 .AND. INCX.EQ.1) X(N/2) = ZERO
      GO TO IGO3
*
  380 CONTINUE
C     PROCEDURE (COMPUTE A CORRECT RESULT)
C     COMPUTE THE CONDITION NUMBER TO USE AS GAUGE FOR ACCURATE RESULTS.
C     THIS IS RETURNED IN XT(*).
C     COMPUTE THE APPROXIMATE CORRECT RESULT.
C     THIS IS RETURNED IN YT(*).
      DO 400 I = 1,N
         YT(I) = ZERO
         XT(I) = ZERO
         IF (INCX.LT.0) THEN
             IBX = (1-N)*INCX + 1
*
         ELSE
             IBX = 1
         END IF
*
         DO 390 J = 1,N
            T = XS(IBX+ (J-1)*INCX)
            IF (ISNUM.GE.9) T = X(IBX+ (J-1)*INCX)
            YT(I) = YT(I) + AS(1+ (I-1)*INCRA+ (J-1)*INCCA)*T
            XT(I) = XT(I) + AS(1+ (I-1)*INCRA+ (J-1)*INCCA)**2
  390    CONTINUE
         XT(I) = SQRT(XT(I))
  400 CONTINUE
      XN = ZERO
      DO 410 J = 1,N
         T = XS(IBX+ (J-1)*INCX)
         IF (ISNUM.GE.9) T = X(IBX+ (J-1)*INCX)
         XN = XN + T**2
  410 CONTINUE
      XN = SQRT(XN) 
C     COMPUTE THE GAUGES FOR THE RESULTS.
      DO 420 I = 1,N
         XT(I) = XT(I)*XN
  420 CONTINUE
C     COMPUTE THE DIFFERENCES. THEY SHOULD BE SMALL FOR CORRECT RESULTS.
      DO 430 I = 1,N
         T = X(IBX+ (I-1)*INCX)
         IF (ISNUM.GE.9) T = XS(IBX+ (I-1)*INCX)
         YT(I) = YT(I) - T
  430 CONTINUE
C     COMPUTE THE GRADE OF THIS RESULT. 
      IGR = 0
      T = ONE
  440 CONTINUE
C     THIS TEST ALLOWS UP TO A LOSS OF FULL PRECISION BEFORE QUITTING.
      IF (IGR.GE.IG) GO TO 470
      DO 450 I = 1,N
         IF (SDIFF(T*ABS(YT(I))+XT(I),XT(I)).EQ.ZERO) GO TO 450
         T = T*HALF 
         IGR = IGR + 1
         GO TO 440
*
  450 CONTINUE
C     IF THE LOOP COMPLETES, ALL VALUES ARE 'SMALL.'  THE VALUE IGR/IG
C     IS THE GRADE ASSIGNED.  THE VALUE OF IGR IS MAXIMIZED OVER ALL THE
C     PROBLEMS.
  460 CONTINUE
      AVIGR = MAX(AVIGR,REAL(IGR))
      GO TO IGO4
*
  470 CONTINUE
      FATAL = .TRUE.
      GO TO 460
*
*     LAST EXECUTABLE LINE OF SCHCK3
 9001 FORMAT (' IN SUBPROGRAM ',A,/,' TESTS ACTIVE WITH OPTIONS = ',
     .  3 (A,2X),/,' N = ',I4,/,' INCX = ',I2,/,' K =',I4)
 9011 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WITH INVALID INPUT.',/,' OPTIONS = ',3 (A,2X),/,
     .  ' N = ',I4,/,' INCX = ',I2,/,' K = ',I4)
 9021 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WHILE COMPUTING',/,' OPTIONS = ',3 (A,2X),/,
     .  ' N = ',I4,/,' INCX = ',I2,/,' K = ',I4)
 9031 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'RECEIVED A LOSS GRADE OF ', 
     .  F5.2,' OUT OF ',I3)
 9041 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'FAILED.')
 9051 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'NOT TESTED.')
      END 
      SUBROUTINE SCHCK4(ISNUM,SNAME,IG,DOPE,NUNIT,AVIGR,FATAL)
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     TEST SGER, 12.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C     THIS PROGRAM HAS TWO PARTS.  THE FIRST PART CHECKS TO SEE
C     IF ANY DATA GETS CHANGED WHEN NONE SHOULD.  FOR EXAMPLE WHEN
C     USING AN INVALID OPTION OR NONPOSITVE PROBLEM DIMENSIONS.
C     THE SECOND PART CHECKS THAT THE RESULTS ARE REASONABLY ACCURATE.
C     DIMENSION AND PROBLEM SIZE DATA.. 
      INTEGER INC(04),IDIM(08),NUNIT(2) 
      REAL ALF(04),SDIFF
      LOGICAL ISAME(13),LSE,FATAL,SAME,NCHNG,RESET
      CHARACTER *128 DOPE(2)
      CHARACTER *6 SNAME
      INTEGER LA,LV 
      PARAMETER (LA=4096,LV=4096,LMN=2048)
      REAL A(LA),AS(LA),X(LV),XS(LV)
      REAL Y(LV),YS(LV),YT(LMN),XT(LMN) 
      PARAMETER (ZERO=0.E0,HALF=.5E0,ONE=1.E0)
      COMMON /ARRAYS/AR,AS,X,XS,Y,YS,YT,XT
      EXTERNAL SDIFF
      INTEGER LIJU_ZERO(3)
      DATA LIJU_ZERO/3*0/
*
      DATA ALF/-1.E0,2.E0,.3E0,1.E0/
      DATA INC/-2,-1,1,2/
      DATA IDIM/1,2,4,8,64,128,2048,0/
      FATAL = .FALSE.
C     CHECK GENERAL RANK 1 UPDATE, 12.
      IF (ISNUM.LT.0) GO TO 200
      NC = 0
      RESET = .TRUE.
      AVIGR = ZERO
      IX = 0
   10 IX = IX + 1
      IF (IX.GT.4) GO TO 180
      INCX = INC(IX)
      ALPHA = ALF(IX)
      IY = 0
   20 IY = IY + 1
      IF (IY.GT.4) GO TO 170
      INCY = INC(IY)
      MM = 0
   30 MM = MM + 1
      IF (MM.GT.8) GO TO 160
      M = IDIM(MM)
      NN = 0
   40 NN = NN + 1
      IF (NN.GT.8) GO TO 150
      N = IDIM(NN)
      IF (FATAL) GO TO 190
      ML = N
      NL = M
      INCCA = M
      INCRA = 1
C     DEFINE THE NUMBER OF ARGUMENTS AND THE A ARGUMENT NUMBER.
      LDA = MAX(M,1)
      NARGS = 09
      IAARG = 08
C     IF NOT ENOUGH STORAGE, SKIP THIS CASE. (AVOID EXPLICT M*N).
      IF (SQRT(REAL(N))*SQRT(REAL(M)).GT.SQRT(REAL(LA))) GO TO 40
C     DO (PREPARE NOTES FOR THIS TEST)
C
C     PRINT A MESSAGE THAT GIVES DEBUGGING INFORMATION.  THIS
C     MESSAGE SAYS..
C     IN SUBPROGRAM XXXXX TESTS WERE ACTIVE WITH
C      M = IIII,     N = IIII,
C     INCX = IIII,  INCY = IIII,
C     THE MAIN IDEA HERE IS THAT A SERIOUS BUG THAT OCCURS DURING THE 
C     TESTING OF THESE SUBPROGRAMS MAY LOSE PROGRAM CONTROL.  THIS
C     AUXILLIARY FILE CONTAINS THE DIMENSIONS THAT RESULTED IN THE LOSS
C     OF CONTROL.  HENCE IT PROVIDES THE IMPLEMENTOR WITH MORE COMPLETE
C     INFORMATION ABOUT WHERE TO START TRACKING DOWN THE BUG.
      IF (NUNIT(1).GT.0) THEN 
C     IF UNIT IS NOT AVAILABLE WITH 'NEW' STATUS, OPEN WITH 
C     'OLD' AND THEN DELETE IT.
          ISTAT = 1 
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.1) GO TO 50
C     GET RID OF ANY OLD FILE.
          CLOSE (UNIT=NUNIT(1),STATUS='DELETE',ERR=50)
   50     CONTINUE
          ISTAT = 2 
C     CREATE A NEW FILE FOR THE NEXT TEST.
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.0) GO TO 70
   60     CONTINUE
          NMESS = 7 
C     DO (PRINT A MESSAGE)
C     PRINT AN ERROR MESSAGE ABOUT OPENING THE NAME FILE.
          CALL SMESSG(LIJU_ZERO,1,NMESS)
          FATAL = .TRUE.
          GO TO 190 
*
   70     CONTINUE
          WRITE (NUNIT(1),9001) SNAME,M,N,INCX,INCY
C     CLOSE THE FILE SO USEFUL STATUS INFORMATION IS SEALED.
          CLOSE (UNIT=NUNIT(1))
      END IF
C     DO (DEFINE A SET OF PROBLEM DATA) 
      ASSIGN 80 TO IGO3
      GO TO 270
*
   80 CONTINUE
C     DO (CALL SUBROUTINE)
      ASSIGN 90 TO IGO1
      GO TO 230
*
   90 CONTINUE
      IF (M.LE.0 .OR. N.LE.0) THEN
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 100 TO IGO2
          GO TO 220 
*
  100     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 110 I = 1,NARGS
             SAME = SAME .AND. ISAME(I) 
             IF ( .NOT. ISAME(I)) THEN
                 WRITE (NUNIT(2),9011) SNAME,I,M,N,INCX,INCY
             END IF 
*
  110     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 190
*
          END IF
*
      ELSE
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 120 TO IGO2
          GO TO 220 
*
  120     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 130 I = 1,NARGS
             NCHNG = (I.EQ.IAARG .OR. ISAME(I))
             SAME = SAME .AND. NCHNG
             IF ( .NOT. NCHNG) THEN
                 WRITE (NUNIT(2),9021) SNAME,I,M,N,INCX,INCY
             END IF 
*
  130     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 190
*
          END IF
*
          NC = NC + 1
C     DO (COMPUTE A CORRECT RESULT)
          ASSIGN 140 TO IGO4
          GO TO 280 
*
  140     CONTINUE
C     IF GOT REALLY BAD ANSWER, PRINT NOTE AND GO BACK.
          IF (FATAL) GO TO 180
*
      END IF
*
      GO TO 40
*
  150 CONTINUE
      GO TO 30
*
  160 CONTINUE
      GO TO 20
*
  170 CONTINUE
      GO TO 10
*
  180 CONTINUE
C     REPORT ON ACCURACY OF DATA.
      WRITE (NUNIT(2),9031) ISNUM,SNAME,AVIGR,IG
      GO TO 210
*
  190 CONTINUE
      WRITE (NUNIT(2),9041) ISNUM,SNAME 
      GO TO 210
*
  200 CONTINUE
      WRITE (NUNIT(2),9051) - ISNUM,SNAME
  210 CONTINUE
      RETURN
*
  220 CONTINUE
C     PROCEDURE (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
      ISAME(1) = MS .EQ. M
      ISAME(2) = NS .EQ. N
      ISAME(3) = ALS .EQ. ALPHA
      ISAME(4) = .TRUE.
      IF (NL.GT.0 .AND. INCX.NE.0) ISAME(4) = LSE(XS,X,1,NL,ABS(INCX))
      ISAME(5) = INCXS .EQ. INCX
      ISAME(6) = .TRUE.
      IF (ML.GT.0 .AND. INCY.NE.0) ISAME(6) = LSE(YS,Y,1,ML,ABS(INCY))
      ISAME(7) = INCYS .EQ. INCY
      ISAME(8) = .TRUE.
      IF (M.GT.0 .AND. N.GT.0) ISAME(8) = LSE(AS,A,M,N,LDA) 
      ISAME(9) = LDAS .EQ. LDA
*
      GO TO IGO2
*
  230 CONTINUE
C     PROCEDURE (CALL SUBROUTINE)
C     SAVE EVERY DATUM BEFORE THE CALL. 
      MS = M
      NS = N
      ALS = ALPHA
      DO 240 I = 1,M*N
         AS(I) = A(I)
  240 CONTINUE
      LDAS = LDA
C     SAVE COPY OF THE X AND Y VECTORS. 
      IBX = 1
      IF (INCX.LT.0) IBX = 1 + (1-NL)*INCX
      DO 250 J = 1,NL
         XS(IBX+ (J-1)*INCX) = X(IBX+ (J-1)*INCX) 
  250 CONTINUE
      INCXS = INCX
      IBY = 1
      IF (INCY.LT.0) IBY = 1 + (1-ML)*INCY
      DO 260 I = 1,ML
         YS(IBY+ (I-1)*INCY) = Y(IBY+ (I-1)*INCY) 
  260 CONTINUE
      INCYS = INCY
      CALL SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
      GO TO IGO1
*
  270 CONTINUE
C     PROCEDURE (DEFINE A SET OF PROBLEM DATA)
C     DO NOTHING IF BOTH DIMENSIONS ARE NOT POSITIVE.
      IF (M.LE.0 .OR. N.LE.0) GO TO IGO3
      TRANSL = ZERO 
      CALL SMAKE(A,M,N,LDA,RESET,TRANSL)
*
      TRANSL = 500.E0
      RESET = .FALSE.
      CALL SMAKE(X,1,NL,MAX(1,ABS(INCX)),RESET,TRANSL)
      IF (NL.GT.1 .AND. INCX.EQ.1) X(NL/2) = ZERO 
      TRANSL = ZERO 
      CALL SMAKE(Y,1,ML,MAX(1,ABS(INCY)),RESET,TRANSL)
      GO TO IGO3
*
  280 CONTINUE
C     PROCEDURE (COMPUTE A CORRECT RESULT)
C     COMPUTE THE CONDITION NUMBER TO USE AS GAUGE FOR ACCURATE RESULTS.
C     THIS IS RETURNED IN XT(*).
C     COMPUTE THE APPROXIMATE CORRECT RESULT.
C     THIS IS RETURNED IN YT(*), COLUMN BY COLUMN.
      IF (INCY.LT.0) THEN
          IBY = (1-ML)*INCY + 1
*
      ELSE
          IBY = 1
      END IF
*
      DO 340 J = 1,N
         DO 290 I = 1,M
            IF (INCX.LT.0) THEN
                IBX = (1-NL)*INCX + 1
*
            ELSE
                IBX = 1
            END IF
*
            YT(I) = AS(1+ (I-1)*INCRA+ (J-1)*INCCA) +
     .              ALPHA*XS(IBX+ (I-1)*INCX)*YS(IBY+ (J-1)*INCY)
            XT(I) = AS(1+ (I-1)*INCRA+ (J-1)*INCCA)**2 +
     .              ALPHA**2*XS(IBX+ (I-1)*INCX)**2*
     .              YS(IBY+ (J-1)*INCY)**2
C     COMPUTE THE GAUGES FOR THE RESULTS.
            XT(I) = SQRT(XT(I))
  290    CONTINUE
C     COMPUTE THE DIFFERENCES. THEY SHOULD BE SMALL FOR CORRECT RESULTS.
         DO 300 I = 1,M
            YT(I) = YT(I) - A(1+ (I-1)*INCRA+ (J-1)*INCCA)
  300    CONTINUE
C     COMPUTE THE GRADE OF THIS RESULT. 
         IGR = 0
         T = ONE
  310    CONTINUE
C     THIS TEST ALLOWS UP TO A LOSS OF FULL PRECISION BEFORE QUITTING.
         IF (IGR.GE.IG) GO TO 360
         DO 320 I = 1,M
            IF (SDIFF(T*ABS(YT(I))+XT(I),XT(I)).EQ.ZERO) GO TO 320
            T = T*HALF
            IGR = IGR + 1
            GO TO 310
*
  320    CONTINUE
C     IF THE LOOP COMPLETES, ALL VALUES ARE 'SMALL.'  THE VALUE IGR/IG
C     IS THE GRADE ASSIGNED.  THE VALUE OF IGR IS MAXIMIZED OVER ALL THE
C     PROBLEMS.
  330    CONTINUE
  340 CONTINUE
  350 AVIGR = MAX(AVIGR,REAL(IGR))
      GO TO IGO4
*
  360 CONTINUE
      FATAL = .TRUE.
      GO TO 350
*
*     LAST EXECUTABLE LINE OF SCHCK4
 9001 FORMAT (' IN SUBPROGRAM ',A,/,'  M =',I4,', N = ',I4,/,' INCX = ',
     .  I2,', INCY = ',I2)
 9011 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WITH INVALID INPUT.',/,' M =',I4,', N = ',I4,/, 
     .  ' INCX = ',I2,', INCY = ',I2)
 9021 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WHILE COMPUTING',/,' M =',I4,', N = ',I4,/,
     .  ' INCX = ',I2,', INCY = ',I2)
 9031 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'RECEIVED A LOSS GRADE OF ', 
     .  F5.2,' OUT OF ',I3)
 9041 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'FAILED.')
 9051 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'NOT TESTED.')
      END 
      SUBROUTINE SCHCK5(ISNUM,SNAME,IG,DOPE,NUNIT,AVIGR,FATAL)
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     TEST SSYR, 13, SSPR, 14, SSYR2, 15, AND SSPR2,16.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C     THIS PROGRAM HAS TWO PARTS.  THE FIRST PART CHECKS TO SEE
C     IF ANY DATA GETS CHANGED WHEN NONE SHOULD.  FOR EXAMPLE WHEN
C     USING AN INVALID OPTION OR NONPOSITVE PROBLEM DIMENSIONS.
C     THE SECOND PART CHECKS THAT THE RESULTS ARE REASONABLY ACCURATE.
C     DIMENSION AND PROBLEM SIZE DATA.. 
      INTEGER INC(04),IDIM(06),NUNIT(2) 
      REAL ALF(04)
      LOGICAL ISAME(13),LSE,FATAL,SAME,NCHNG,RESET
      CHARACTER *128 DOPE(2)
      CHARACTER *6 SNAME
      CHARACTER *3 ICH
      CHARACTER *1 ICHS,ICI
      INTEGER LA,LV 
      PARAMETER (LA=4096,LV=4096,LMN=2048)
      REAL A(LA),AS(LA),X(LV),XS(LV)
      REAL Y(LV),YS(LV),YT(LMN),XT(LMN) 
      PARAMETER (ZERO=0.E0,HALF=.5E0,ONE=1.E0)
      COMMON /ARRAYS/AR,AS,X,XS,Y,YS,YT,XT
      EXTERNAL SDIFF
      INTEGER LIJU_ZERO(3)
      DATA LIJU_ZERO/3*0/
*
      DATA ALF/-1.E0,2.E0,.3E0,1.E0/
      DATA INC/-2,-1,1,2/
      DATA IDIM/1,2,4,8,64,0/ 
      DATA ICH/'LU/'/
      FATAL = .FALSE.
C     CHECK SYMMETRIC MATRIX RANK 1 AND RANK 2 UPDATES.
      IF (ISNUM.LT.0) GO TO 200
      NC = 0
      RESET = .TRUE.
      AVIGR = ZERO
      IX = 0
   10 IX = IX + 1
      IF (IX.GT.4) GO TO 180
      INCX = INC(IX)
      ALPHA = ALF(IX)
      IY = 0
   20 IY = IY + 1
      IF (IY.GT.4) GO TO 170
      INCY = INC(IY)
      NN = 0
   30 NN = NN + 1
      IF (NN.GT.6) GO TO 160
      N = IDIM(NN)
      IC = 0
   40 IC = IC + 1
      IF (IC.GT.3) GO TO 150
      IF (FATAL) GO TO 190
      ICI = ICH(IC:IC)
C     DEFINE THE NUMBER OF ARGUMENTS AND THE Y ARGUMENT NUMBER.
      LDA = MAX(N,1)
      IF (ISNUM.EQ.13) THEN
          NARGS = 07
          IAARG = 06
*
      ELSE IF (ISNUM.EQ.14) THEN
          NARGS = 06
          IAARG = 06
*
      ELSE IF (ISNUM.EQ.15) THEN
          NARGS = 9 
          IAARG = 8 
*
      ELSE IF (ISNUM.EQ.16) THEN
          NARGS = 8 
          IAARG = 8 
      END IF
C     DO (PREPARE NOTES FOR THIS TEST)
C
C     PRINT A MESSAGE THAT GIVES DEBUGGING INFORMATION.  THIS
C     MESSAGE SAYS..
C     IN SUBPROGRAM XXXXX TESTS WERE ACTIVE WITH
C     OPTION = 'A'
C     N = IIII,
C     INCX = IIII,  INCY = IIII,
C     THE MAIN IDEA HERE IS THAT A SERIOUS BUG THAT OCCURS DURING THE 
C     TESTING OF THESE SUBPROGRAMS MAY LOSE PROGRAM CONTROL.  THIS
C     AUXILLIARY FILE CONTAINS THE DIMENSIONS THAT RESULTED IN THE LOSS
C     OF CONTROL.  HENCE IT PROVIDES THE IMPLEMENTOR WITH MORE COMPLETE
C     INFORMATION ABOUT WHERE TO START TRACKING DOWN THE BUG.
      IF (NUNIT(1).GT.0) THEN 
C     IF UNIT IS NOT AVAILABLE WITH 'NEW' STATUS, OPEN WITH 
C     'OLD' AND THEN DELETE IT.
          ISTAT = 1 
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.1) GO TO 50
C     GET RID OF ANY OLD FILE.
          CLOSE (UNIT=NUNIT(1),STATUS='DELETE',ERR=50)
   50     CONTINUE
          ISTAT = 2 
C     CREATE A NEW FILE FOR THE NEXT TEST.
          CALL SOPEN(NUNIT(1),DOPE(1),ISTAT,IERROR)
          IF (IERROR.EQ.0) GO TO 70
   60     CONTINUE
          NMESS = 7 
C     DO (PRINT A MESSAGE)
C     PRINT AN ERROR MESSAGE ABOUT OPENING THE NAME FILE.
          CALL SMESSG(LIJU_ZERO,1,NMESS)
          FATAL = .TRUE.
          GO TO 190 
*
   70     CONTINUE
          WRITE (NUNIT(1),9001) SNAME,ICI,N,INCX,INCY
C     CLOSE THE FILE SO USEFUL STATUS INFORMATION IS SEALED.
          CLOSE (UNIT=NUNIT(1))
      END IF
C     DO (DEFINE A SET OF PROBLEM DATA) 
      ASSIGN 80 TO IGO3
      GO TO 370
*
   80 CONTINUE
C     DO (CALL SUBROUTINE)
      ASSIGN 90 TO IGO1
      GO TO 290
*
   90 CONTINUE
      IF (N.LE.0 .OR. ICHAR(ICI).EQ.ICHAR('/')) THEN
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 100 TO IGO2
          GO TO 220 
*
  100     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 110 I = 1,NARGS
             SAME = SAME .AND. ISAME(I) 
             IF ( .NOT. ISAME(I)) THEN
                 WRITE (NUNIT(2),9011) SNAME,I,ICI,N,INCX,INCY
             END IF 
*
  110     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 190
*
          END IF
*
      ELSE
C     DO (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
          ASSIGN 120 TO IGO2
          GO TO 220 
*
  120     CONTINUE
C     IF DATA WAS INCORRECTLY CHANGED, MAKE NOTES AND RETURN.
          SAME = .TRUE.
          DO 130 I = 1,NARGS
             NCHNG = (I.EQ.IAARG .OR. ISAME(I))
             SAME = SAME .AND. NCHNG
             IF ( .NOT. NCHNG) THEN
                 WRITE (NUNIT(2),9021) SNAME,I,ICI,N,INCX,INCY
             END IF 
*
  130     CONTINUE
          IF ( .NOT. SAME) THEN
              FATAL = .TRUE.
              GO TO 190
*
          END IF
*
          NC = NC + 1
C     DO (COMPUTE A CORRECT RESULT)
          ASSIGN 140 TO IGO4
          GO TO 400 
*
  140     CONTINUE
C     IF GOT REALLY BAD ANSWER, PRINT NOTE AND GO BACK.
          IF (FATAL) GO TO 180
*
      END IF
*
      GO TO 40
*
  150 CONTINUE
      GO TO 30
*
  160 CONTINUE
      IF (ISNUM.GE.15) GO TO 20
      GO TO 10
*
  170 CONTINUE
      GO TO 10
*
  180 CONTINUE
C     REPORT ON ACCURACY OF DATA.
      WRITE (NUNIT(2),9031) ISNUM,SNAME,AVIGR,IG
      GO TO 210
*
  190 CONTINUE
      WRITE (NUNIT(2),9041) ISNUM,SNAME 
      GO TO 210
*
  200 CONTINUE
      WRITE (NUNIT(2),9051) - ISNUM,SNAME
  210 CONTINUE
      RETURN
*
  220 CONTINUE
C     PROCEDURE (SEE WHAT DATA CHANGED INSIDE SUBROUTINES)
      IF (ISNUM.EQ.13) THEN
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = ALS .EQ. ALPHA
          ISAME(4) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(4) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(5) = INCXS .EQ. INCX
          ISAME(6) = .TRUE.
          IF (N.GT.0) ISAME(6) = LSE(AS,A,N,N,LDA)
          ISAME(7) = LDAS .EQ. LDA
*
      ELSE IF (ISNUM.EQ.14) THEN
C     COMPARE THE MATRIX IN THE DATA STRUCTURES WITH THE SAVED COPY.
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = ALS .EQ. ALPHA
          ISAME(4) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(4) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(5) = INCXS .EQ. INCX
          ISAME(6) = .TRUE.
          IOFF = 0
          DO 240 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 230 I = ISTRT,IEND
                IOFF = IOFF + 1
                IF (A(IOFF).NE.AS(1+ (I-1)+ (J-1)*LDA)) THEN
                    ISAME(6) = .FALSE.
                    GO TO 250 
*
                END IF
*
  230        CONTINUE
  240     CONTINUE
  250     CONTINUE
*
      ELSE IF (ISNUM.EQ.15) THEN
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = ALS .EQ. ALPHA
          ISAME(4) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(4) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(5) = INCXS .EQ. INCX
          ISAME(6) = .TRUE.
          IF (N.GT.0 .AND. INCY.NE.0) ISAME(6) = LSE(YS,Y,1,N,ABS(INCY))
          ISAME(7) = INCYS .EQ. INCY
          ISAME(8) = .TRUE.
          IF (N.GT.0) ISAME(8) = LSE(AS,A,N,N,LDA)
          ISAME(9) = LDAS .EQ. LDA
*
      ELSE IF (ISNUM.EQ.16) THEN
          ISAME(1) = ICHAR(ICI) .EQ. ICHAR(ICHS)
          ISAME(2) = NS .EQ. N
          ISAME(3) = ALS .EQ. ALPHA
          ISAME(4) = .TRUE.
          IF (N.GT.0 .AND. INCX.NE.0) ISAME(4) = LSE(XS,X,1,N,ABS(INCX))
          ISAME(5) = INCXS .EQ. INCX
          ISAME(6) = .TRUE.
          IF (N.GT.0 .AND. INCY.NE.0) ISAME(6) = LSE(YS,Y,1,N,ABS(INCY))
          ISAME(7) = INCYS .EQ. INCY
          ISAME(8) = .TRUE.
          IOFF = 0
          DO 270 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 260 I = ISTRT,IEND
                IOFF = IOFF + 1
                IF (A(IOFF).NE.AS(1+ (I-1)+ (J-1)*LDA)) THEN
                    ISAME(8) = .FALSE.
                    GO TO 280 
*
                END IF
*
  260        CONTINUE
  270     CONTINUE
  280     CONTINUE
      END IF
*
      GO TO IGO2
*
  290 CONTINUE
C     PROCEDURE (CALL SUBROUTINE)
C     SAVE EVERY DATUM BEFORE THE CALL. 
      ICHS = ICI
      NS = N
      ALS = ALPHA
      DO 300 I = 1,N*N
         AS(I) = A(I)
  300 CONTINUE
      LDAS = LDA
C     SAVE COPY OF THE X AND Y VECTORS. 
      IBX = 1
      IF (INCX.LT.0) IBX = 1 + (1-N)*INCX
      DO 310 J = 1,N
         XS(IBX+ (J-1)*INCX) = X(IBX+ (J-1)*INCX) 
  310 CONTINUE
      INCXS = INCX
      IBY = 1
      IF (INCY.LT.0) IBY = 1 + (1-N)*INCY
      DO 320 I = 1,N
         YS(IBY+ (I-1)*INCY) = Y(IBY+ (I-1)*INCY) 
  320 CONTINUE
      INCYS = INCY
      IF (ISNUM.EQ.13) THEN
          CALL SSYR(ICI,N,ALPHA,X,INCX,A,LDA)
*
      ELSE IF (ISNUM.EQ.14) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH SSPR.
          IOFF = 0
          DO 340 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 330 I = ISTRT,IEND
                IOFF = IOFF + 1
                A(IOFF) = AS(1+ (I-1)+ (J-1)*LDA) 
  330        CONTINUE
*
  340     CONTINUE
          CALL SSPR(ICI,N,ALPHA,X,INCX,A)
*
      ELSE IF (ISNUM.EQ.15) THEN
*
          CALL SSYR2(ICI,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
      ELSE IF (ISNUM.EQ.16) THEN
C     TRANSFER THE MATRIX TO THE DATA STRUCTURE USED WITH SSPR2.
          IOFF = 0
          DO 360 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 350 I = ISTRT,IEND
                IOFF = IOFF + 1
                A(IOFF) = AS(1+ (I-1)+ (J-1)*LDA) 
  350        CONTINUE
*
  360     CONTINUE
          CALL SSPR2(ICI,N,ALPHA,X,INCX,Y,INCY,A) 
      END IF
*
      GO TO IGO1
*
  370 CONTINUE
C     PROCEDURE (DEFINE A SET OF PROBLEM DATA)
C     DO NOTHING IF DIMENSIONS ARE NOT POSITIVE.
      IF (N.LE.0) GO TO IGO3
      TRANSL = ZERO 
      CALL SMAKE(A,N,N,LDA,RESET,TRANSL)
C     MAKE THE DATA MATRIX SYMMETRIC.
      DO 390 I = 1,N
         DO 380 J = I,N
            T = (A(1+ (I-1)+ (J-1)*LDA)+A(1+ (J-1)+ (I-1)*LDA))*HALF
            A(1+ (I-1)+ (J-1)*LDA) = T
            A(1+ (J-1)+ (I-1)*LDA) = T
  380    CONTINUE
  390 CONTINUE
*
      TRANSL = 500.E0
      RESET = .FALSE.
      CALL SMAKE(X,1,N,MAX(1,ABS(INCX)),RESET,TRANSL)
      IF (N.GT.1 .AND. INCX.EQ.1) X(N/2) = ZERO
      TRANSL = ZERO 
      CALL SMAKE(Y,1,N,MAX(1,ABS(INCY)),RESET,TRANSL)
      GO TO IGO3
*
  400 CONTINUE
C     PROCEDURE (COMPUTE A CORRECT RESULT)
C     COMPUTE THE CONDITION NUMBER TO USE AS GAUGE FOR ACCURATE RESULTS.
C     THIS IS RETURNED IN XT(*).
C     COMPUTE THE APPROXIMATE CORRECT RESULT.
      IF (ISNUM.EQ.13 .OR. ISNUM.EQ.14) THEN
          IF (INCX.LT.0) THEN 
              IBX = (1-N)*INCX + 1
*
          ELSE
              IBX = 1
          END IF
*
          IOFF = 0
          DO 450 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 410 I = ISTRT,IEND
                YT(I) = AS(1+ (I-1)+ (J-1)*LDA) + 
     .                  ALPHA*XS(IBX+ (J-1)*INCX)*XS(IBX+ (I-1)*INCX) 
                XT(I) = AS(1+ (I-1)+ (J-1)*LDA)**2 +
     .                  ALPHA**2*XS(IBX+ (I-1)*INCX)**2*
     .                  XS(IBX+ (J-1)*INCX)**2
  410        CONTINUE
C     COMPUTE THE DIFFERENCES. THEY SHOULD BE SMALL FOR CORRECT RESULTS.
             DO 420 I = ISTRT,IEND
                XT(I) = SQRT(XT(I))
                IF (ISNUM.EQ.13) THEN
                    YT(I) = YT(I) - A(1+ (I-1)+ (J-1)*LDA)
*
                ELSE IF (ISNUM.EQ.14) THEN
                    IOFF = IOFF + 1
                    YT(I) = YT(I) - A(IOFF)
                END IF
*
  420        CONTINUE
C     COMPUTE THE GRADE OF THIS RESULT. 
             IGR = 0
             T = ONE
             DO 440 I = ISTRT,IEND
  430           CONTINUE
C     THIS TEST ALLOWS UP TO A LOSS OF FULL PRECISION BEFORE QUITTING.
                IF (IGR.GE.IG) GO TO 520
                IF (SDIFF(T*ABS(YT(I))+XT(I),XT(I)).EQ.ZERO) GO TO 440
                T = T*HALF
                IGR = IGR + 1 
                GO TO 430
*
C     IF THE LOOP COMPLETES, ALL VALUES ARE 'SMALL.'  THE VALUE IGR/IG
C     IS THE GRADE ASSIGNED.  THE VALUE OF IGR IS MAXIMIZED OVER ALL THE
C     PROBLEMS.
  440        CONTINUE
  450     CONTINUE
*
      ELSE IF (ISNUM.EQ.15 .OR. ISNUM.EQ.16) THEN 
          IF (INCX.LT.0) THEN 
              IBX = (1-N)*INCX + 1
*
          ELSE
              IBX = 1
          END IF
*
          IF (INCY.LT.0) THEN 
              IBY = (1-N)*INCY + 1
*
          ELSE
              IBY = 1
          END IF
*
          IOFF = 0
          DO 500 J = 1,N
             IF (ICHAR(ICI).EQ.ICHAR('U')) THEN
                 ISTRT = 1
                 IEND = J
*
             ELSE
                 ISTRT = J
                 IEND = N
             END IF 
*
             DO 460 I = ISTRT,IEND
                YT(I) = AS(1+ (I-1)+ (J-1)*LDA) + 
     .                  ALPHA*XS(IBX+ (J-1)*INCX)*YS(IBY+ (I-1)*INCY) +
     .                  ALPHA*XS(IBX+ (I-1)*INCX)*YS(IBY+ (J-1)*INCY) 
                XT(I) = AS(1+ (I-1)+ (J-1)*LDA)**2 +
     .                  ALPHA**2*XS(IBX+ (I-1)*INCX)**2*
     .                  YS(IBY+ (J-1)*INCY)**2 +
     .                  ALPHA**2*XS(IBX+ (J-1)*INCX)**2*
     .                  YS(IBY+ (I-1)*INCY)**2
  460        CONTINUE
C     COMPUTE THE DIFFERENCES. THEY SHOULD BE SMALL FOR CORRECT RESULTS.
             DO 470 I = ISTRT,IEND
                XT(I) = SQRT(XT(I))
                IF (ISNUM.EQ.15) THEN
                    YT(I) = YT(I) - A(1+ (I-1)+ (J-1)*LDA)
*
                ELSE IF (ISNUM.EQ.16) THEN
                    IOFF = IOFF + 1
                    YT(I) = YT(I) - A(IOFF)
                END IF
*
  470        CONTINUE
C     COMPUTE THE GRADE OF THIS RESULT. 
             IGR = 0
             T = ONE
             DO 490 I = ISTRT,IEND
  480           CONTINUE
C     THIS TEST ALLOWS UP TO A LOSS OF FULL PRECISION BEFORE QUITTING.
                IF (IGR.GE.IG) GO TO 520
                IF (SDIFF(T*ABS(YT(I))+XT(I),XT(I)).EQ.ZERO) GO TO 490
                T = T*HALF
                IGR = IGR + 1 
                GO TO 480
*
C     IF THE LOOP COMPLETES, ALL VALUES ARE 'SMALL.'  THE VALUE IGR/IG
C     IS THE GRADE ASSIGNED.  THE VALUE OF IGR IS MAXIMIZED OVER ALL THE
C     PROBLEMS.
  490        CONTINUE
  500     CONTINUE
      END IF
*
  510 CONTINUE
      AVIGR = MAX(AVIGR,REAL(IGR))
      GO TO IGO4
*
  520 CONTINUE
      FATAL = .TRUE.
      GO TO 510
*
*     LAST EXECUTABLE LINE OF SCHCK5
 9001 FORMAT (' IN SUBPROGRAM ',A,/,' TESTS ACTIVE WITH OPTION = ',A,/,
     .  ' N = ',I4,/,' INCX = ',I2,', INCY = ',I2)
 9011 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WITH INVALID INPUT.',/,' OPTION = ',A,/,' N = ',
     .  I4,/,' INCX = ',I2,', INCY = ',I2)
 9021 FORMAT (' IN SUBPROGRAM ',A,/,' ARGUMENT ',I2,
     .  ' WAS CHANGED WHILE COMPUTING',/,' OPTION = ',A,/,' N = ',I4,/,
     .  ' INCX = ',I2,', INCY = ',I2)
 9031 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'RECEIVED A LOSS GRADE OF ', 
     .  F5.2,' OUT OF ',I3)
 9041 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'FAILED.')
 9051 FORMAT (1X,I2,' SUBPROGRAM ',A,T24,'NOT TESTED.')
      END 
      SUBROUTINE SMAKE(A,M,N,LDA,RESET,TRANS)
C     GENERATE VALUES FOR AN M BY N MATRIX A.
C     RESET THE GENERATOR IF FLAG RESET = .TRUE.
C     TRANSLATE THE VALUES WITH TRANS.
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      REAL A(LDA,*),TRANS,ANOISE
      REAL ZERO,HALF,ONE
      PARAMETER (ZERO=0.E0,HALF=.5E0,ONE=1.E0,THREE=3.E0)
      LOGICAL RESET 
      IF (RESET) THEN
          ANOISE = -ONE
          ANOISE = SBEG(ANOISE)
          ANOISE = ZERO
      END IF
*
      IC = 0
      DO 20 I = 1,M 
         DO 10 J = 1,N
            IC = IC + 1
C     BREAK UP PERIODICITIES THAT ARE MULTIPLES OF 5.
            IF (MOD(IC,5).EQ.0) A(I,J) = SBEG(ANOISE)
            A(I,J) = SBEG(ANOISE) - TRANS
C     HERE THE PERTURBATION IN THE LAST BIT POSITION IS MADE.
            A(I,J) = A(I,J) + ONE/THREE 
            ANOISE = 0.E0
   10    CONTINUE
   20 CONTINUE
      RETURN
*     LAST EXECUTABLE LINE OF SMAKE
      END 
      SUBROUTINE SOPEN(IUNIT,NAME,ISTAT,IERROR)
C     OPEN UNIT IUNIT WITH FILE NAMED NAME.
C     ISTAT=1 FOR 'OLD', =2 FOR 'NEW', =3 FOR 'UNKNOWN'.
C     THE RETURN FLAG IERROR=0 FOR SUCCESS, =1 FOR FAILURE. 
C     A BAD VALUE OF ISTAT CAN ALSO INDICATE FAILURE.
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      CHARACTER * (*)  NAME
      IF (ISTAT.EQ.1) OPEN (UNIT=IUNIT,FILE=NAME,STATUS='OLD',ERR=10) 
      IF (ISTAT.EQ.2) OPEN (UNIT=IUNIT,FILE=NAME,STATUS='NEW',ERR=10) 
      IF (ISTAT.EQ.3) OPEN (UNIT=IUNIT,FILE=NAME,STATUS='UNKNOWN',
     .                     ERR=10)
      GO TO (20,20,20),ISTAT
*
   10 CONTINUE
      IERROR = 1
      GO TO 30
*
   20 CONTINUE
      IERROR = 0
   30 CONTINUE
      RETURN
*     LAST EXECUTABLE LINE OF SOPEN
      END 
      FUNCTION SDIFF(X,Y)
C     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUNE 7
C     APPEARED IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C     THIS IS USED AS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      SDIFF = X - Y 
      RETURN
*     LAST EXECUTABLE LINE OF SDIFF
      END 
*
      FUNCTION SBEG(ANOISE)
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      SAVE
C          GENERATE NUMBERS FOR CONSTRUCTION OF TEST CASES. 
      IF (ANOISE) 10,30,20
   10 MI = 891
      MJ = 457
      I = 7
      J = 7
      AJ = 0.
      SBEG = 0.
      RETURN
*
   20 J = J*MJ
      J = J - 997* (J/997)
      AJ = J - 498
C     THE SEQUENCE OF VALUES OF I  IS BOUNDED BETWEEN 1 AND 999
C     IF INITIAL I = 1,2,3,6,7, OR 9,  THE PERIOD WILL BE 50
C     IF INITIAL I = 4 OR 8   THE PERIOD WILL BE 25
C     IF INITIAL I = 5        THE PERIOD WILL BE 10
   30 I = I*MI
      I = I - 1000* (I/1000)
      AI = I - 500
      SBEG = AI + AJ*ANOISE
      RETURN
*     LAST EXECUTABLE LINE OF SBEG
      END 
*
      LOGICAL FUNCTION LSE(RI,RJ,M,N,LDI)
C     TEST IF TWO REAL ARRAYS ARE IDENTICAL.
C     THE ARRAYS ARE M BY N.
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      REAL RI(LDI,*),RJ(LDI,*)
      DO 20 I = 1,M 
         DO 10 J = 1,N
            IF (RI(I,J).NE.RJ(I,J)) THEN
                LSE = .FALSE. 
                GO TO 30
*
            END IF
*
   10    CONTINUE
   20 CONTINUE
      LSE = .TRUE.
   30 CONTINUE
      RETURN
*     LAST EXECUTABLE LINE OF LSE
      END 
*
      LOGICAL FUNCTION LDE(DI,DJ,M,N,LDI)
C     TEST IF TWO DOUBLE PRECISION ARRAYS ARE IDENTICAL.
C     THE ARRAYS ARE M BY N.
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      DOUBLE PRECISION DI(LDI,*),DJ(LDI,*)
      DO 20 I = 1,M 
         DO 10 J = 1,N
            IF (DI(I,J).NE.DJ(I,J)) THEN
                LDE = .FALSE. 
                GO TO 30
*
            END IF
*
   10    CONTINUE
   20 CONTINUE
      LDE = .TRUE.
   30 CONTINUE
      RETURN
*     LAST EXECUTABLE LINE OF LDE
      END 
*
      LOGICAL FUNCTION LCE(CI,CJ,M,N,LDI)
C     TEST IF TWO COMPLEX ARRAYS ARE IDENTICAL.
C     THE ARRAYS ARE M BY N.
C     THIS IS A TEST SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
      COMPLEX CI(LDI,*),CJ(LDI,*)
      DO 20 I = 1,M 
         DO 10 J = 1,N
            IF (REAL(CI(I,J)).NE.REAL(CJ(I,J)) .OR. AIMAG(CI(I,J)).NE.
     .          AIMAG(CJ(I,J))) THEN
                LCE = .FALSE. 
                GO TO 30
*
            END IF
*
   10    CONTINUE
   20 CONTINUE
      LCE = .TRUE.
   30 CONTINUE
      RETURN
*     LAST EXECUTABLE LINE OF LCE
      END 
C
C***********************************************************************
C
C     File of the REAL              Level 2 BLAS routines:  
C
C      SGEMV, SGBMV, SSYMV, SSBMV, SSPMV, STRMV, STBMV, STPMV,
C      SGER , SSYR , SSPR ,
C      SSYR2, SSPR2,
C      STRSV, STBSV, STPSV.
C
C     See: 
C
C        Dongarra J. J., Du Croz J. J., Hammarling S. and Hanson R. J.. 
C        A proposal for an extended set of Fortran Basic Linear Algebra
C        Subprograms. Technical Memorandum No.41 (revision 1),
C        Mathematics and Computer Science Division, Argone National
C        Laboratory, 9700 South Cass Avenue, Argonne, Illinois 60439,
C        USA, or NAG Technical Report TR4/85, Numerical Algorithms Group
C        Inc., 1101 31st Street, Suite 100, Downers Grove, Illinois
C        60606-1263, USA.
C
C***********************************************************************
C
      SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 TRANS
      INTEGER M,N,LDA,INCX,INCY
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, 
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N'  y := alpha*A*x + beta*y.
*
*              TRANS = 'T'  y := alpha*A'*x + beta*y.
*
*              TRANS = 'C'  y := alpha*A'*x + beta*y
*.
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           max(m,1).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the 
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Note that TRANS, M, N and LDA must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KX,KY,LENX,LENY 
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. ((M.GT.0) .AND. (N.GT.0) .AND.
     .     (LDA.GE.M))
*
*     Quick return if possible.
*
      IF (((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE)) .OR. .NOT. OK) RETURN
*
*     Set LENX and LENY, the lengths of the vectors x and y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
*
      ELSE
          LENX = M
          LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,LENY
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,LENY
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (LENX-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (LENY-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,LENY
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,LENY
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 50,I = 1,M
                        Y(I) = Y(I) + TEMP*A(I,J) 
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IY = KY
                     DO 70,I = 1,M
                        Y(IY) = Y(IY) + TEMP*A(I,J)
                        IY = IY + INCY
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP = ZERO
                 DO 90,I = 1,M
                    TEMP = TEMP + A(I,J)*X(I)
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP
  100         CONTINUE
*
          ELSE
              JY = KY
              DO 120,J = 1,N
                 TEMP = ZERO
                 IX = KX
                 DO 110,I = 1,M
                    TEMP = TEMP + A(I,J)*X(IX)
                    IX = IX + INCX
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SGEMV .
*
      END 
      SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 TRANS
      INTEGER M,N,KL,KU,LDA,INCX,INCY
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SGBMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, 
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n band matrix, with kl sub-diagonals and ku super-diagonals. 
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N'  y := alpha*A*x + beta*y.
*
*              TRANS = 'T'  y := alpha*A'*x + beta*y.
*
*              TRANS = 'C'  y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  KL     - INTEGER.
*           On entry, KL specifies the number of sub-diagonals of the
*           matrix A. KL must satisfy  0 .le. KL.
*           Unchanged on exit.
*
*  KU     - INTEGER.
*           On entry, KU specifies the number of super-diagonals of the
*           matrix A. KU must satisfy   0 .le. KU. 
*           Unchanged on exit.
*
*  Users may find that efficiency of their application is enhanced by
*  adjusting the values of m and n so that KL .ge. max(0,m-n) and
*  KU .ge. max(0,n-m) or KL and KU so that KL .lt. m and KU .lt. n.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading ( kl + ku + 1 ) by n part of the
*           array A must contain the matrix of coefficients, supplied
*           column by column, with the leading diagonal of the matrix in
*           row ( ku + 1 ) of the array, the first super-diagonal
*           starting at position 2 in row ku, the first sub-diagonal
*           starting at position 1 in row ( ku + 2 ), and so on.
*           This placement of the data can be realized with the
*           following loops: 
*               DO 20 J =1,N
*                    K=KU+1-J 
*                    DO 10 I =MAX(1,J-KU),MIN(M,J+KL)
*                         A(K+I,J)=matrix entry of row I, column J. 
*     10             CONTINUE 
*     20        CONTINUE
*           Elements in the array A that do not correspond to elements
*           in the band matrix (such as the top left ku by ku triangle)
*           are not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( kl + ku + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the 
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry, the incremented array Y must contain the 
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-Sept-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER K,KUP1,KX,KY,LENX,LENY
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (M.GT.0) .AND. (N.GT.0) .AND.
     .     (KL.GE.0) .AND. (KU.GE.0) .AND.
     .     (LDA.GE. (KL+KU+1))
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
*
      ELSE
          LENX = M
          LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the band part of A.
*
*     First form  y := beta*y  and set up the start points in  X  and  Y
*     if the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,LENY
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,LENY
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (LENX-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (LENY-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,LENY
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,LENY
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      KUP1 = KU + 1 
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     K = KUP1 - J
                     DO 50,I = MAX(1,J-KU),MIN(M,J+KL)
                        Y(I) = Y(I) + TEMP*A(K+I,J)
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IY = KY
                     K = KUP1 - J
                     DO 70,I = MAX(1,J-KU),MIN(M,J+KL)
                        Y(IY) = Y(IY) + TEMP*A(K+I,J)
                        IY = IY + INCY
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
                 IF (J.GT.KU) KY = KY + INCY
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP = ZERO
                 K = KUP1 - J 
                 DO 90,I = MAX(1,J-KU),MIN(M,J+KL)
                    TEMP = TEMP + A(K+I,J)*X(I)
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP
  100         CONTINUE
*
          ELSE
              JY = KY
              DO 120,J = 1,N
                 TEMP = ZERO
                 IX = KX
                 K = KUP1 - J 
                 DO 110,I = MAX(1,J-KU),MIN(M,J+KL)
                    TEMP = TEMP + A(K+I,J)*X(IX)
                    IX = IX + INCX
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP
                 JY = JY + INCY
                 IF (J.GT.KU) KX = KX + INCX
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SGBMV .
*
      END 
      SUBROUTINE SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 UPLO
      INTEGER N,LDA,INCX,INCY 
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SSYMV  performs the matrix-vector  operation 
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows: 
*
*              UPLO = 'U'          Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L'          Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(n,1).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-Sept-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KX,KY 
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (LDA.GE.N)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,N
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,N
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,N
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,N
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when A is stored in upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 DO 50,I = 1,J - 1
                    Y(I) = Y(I) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(I)
   50            CONTINUE
                 Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
*
          ELSE
              IX = KX - INCX
              DO 80,J = 1,N
                 TEMP1 = ALPHA*X(IX+INCX)
                 TEMP2 = ZERO 
                 IX = KX
                 IY = KY
                 DO 70,I = 1,J - 1
                    Y(IY) = Y(IY) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
   70            CONTINUE
                 Y(IY) = Y(IY) + TEMP1*A(J,J) + ALPHA*TEMP2 
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 Y(J) = Y(J) + TEMP1*A(J,J)
                 DO 90,I = J + 1,N
                    Y(I) = Y(I) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(I)
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 120,J = 1,N
                 TEMP1 = ALPHA*X(JX)
                 TEMP2 = ZERO 
                 Y(JY) = Y(JY) + TEMP1*A(J,J)
                 IX = JX
                 IY = JY
                 DO 110,I = J + 1,N
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(IX)
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP2
                 JX = JX + INCX
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSYMV .
*
      END 
      SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 UPLO
      INTEGER N,K,LDA,INCX,INCY
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SSBMV  performs the matrix-vector  operation 
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric band matrix, with k super-diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the band matrix A is being supplied as
*           follows: 
*
*              UPLO = 'U'          The upper triangular part of A is
*                                  being supplied.
*
*              UPLO = 'L'          The lower triangular part of A is
*                                  being supplied.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry, K specifies the number of super-diagonals of the
*           matrix A. K must satisfy  0 .le. K .lt. n.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the symmetric matrix, supplied column by 
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           Before entry with UPLO = 'L', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the symmetric matrix, supplied column by 
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the 
*           array A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the 
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the 
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KPLUS1,KX,KY,L
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (K.GE.0) .AND. (K.LT.N) .AND. (LDA.GE. (K+1))
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Start the operations. In this version the elements of the array A
*     are accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,N
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,N
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,N
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,N
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when upper triangle of A is stored.
*
          KPLUS1 = K + 1
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 I = MAX(1,J-K)
                 DO 50,L = KPLUS1 + I - J,K
                    Y(I) = Y(I) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(I)
                    I = I + 1 
   50            CONTINUE
                 Y(J) = Y(J) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
   60         CONTINUE
*
          ELSE
              IX = KX - INCX
              DO 80,J = 1,N
                 TEMP1 = ALPHA*X(IX+INCX)
                 TEMP2 = ZERO 
                 IX = KX
                 IY = KY
                 DO 70,L = 1 + MAX(KPLUS1-J,0),K
                    Y(IY) = Y(IY) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
   70            CONTINUE
                 Y(IY) = Y(IY) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
                 IF (J.GT.K) THEN
                     KX = KX + INCX
                     KY = KY + INCY
                 END IF
*
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y  when lower triangle of A is stored.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 Y(J) = Y(J) + TEMP1*A(1,J)
                 I = J + 1
                 DO 90,L = 2,1 + MIN(K,N-J)
                    Y(I) = Y(I) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(I)
                    I = I + 1 
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 120,J = 1,N
                 TEMP1 = ALPHA*X(JX)
                 TEMP2 = ZERO 
                 Y(JY) = Y(JY) + TEMP1*A(1,J)
                 IX = JX
                 IY = JY
                 DO 110,L = 2,1 + MIN(K,N-J)
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(IX)
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP2
                 JX = JX + INCX
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSBMV .
*
      END 
      SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      CHARACTER *1 UPLO
      INTEGER N,INCX,INCY
      REAL ALPHA,AP(*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SSPMV  performs the matrix-vector operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows: 
*
*              UPLO = 'U'          The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L'          The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with UPLO = 'U', the array AP must
*           contain the upper triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) 
*           and a( 2, 2 ) respectively, and so on.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) 
*           and a( 3, 1 ) respectively, and so on.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-Sept-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER K,KK,KX,KY
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,N
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,N
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,N
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,N
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      K = 1
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when AP contains the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 DO 50,I = 1,J - 1
                    Y(I) = Y(I) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(I)
                    K = K + 1 
   50            CONTINUE
                 Y(J) = Y(J) + TEMP1*AP(K) + ALPHA*TEMP2
                 K = K + 1
   60         CONTINUE
*
          ELSE
              IX = KX - INCX
              DO 80,J = 1,N
                 TEMP1 = ALPHA*X(IX+INCX)
                 TEMP2 = ZERO 
                 IX = KX
                 IY = KY
                 KK = K
                 DO 70,K = KK,KK + J - 2
                    Y(IY) = Y(IY) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
   70            CONTINUE
                 Y(IY) = Y(IY) + TEMP1*AP(K) + ALPHA*TEMP2
                 K = K + 1
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y  when AP contains the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 Y(J) = Y(J) + TEMP1*AP(K)
                 K = K + 1
                 DO 90,I = J + 1,N
                    Y(I) = Y(I) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(I)
                    K = K + 1 
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 120,J = 1,N
                 TEMP1 = ALPHA*X(JX)
                 TEMP2 = ZERO 
                 Y(JY) = Y(JY) + TEMP1*AP(K)
                 IX = JX
                 IY = JY
                 KK = K + 1
                 DO 110,K = KK,KK + N - (J+1)
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(IX)
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP2
                 JX = JX + INCX
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSPMV .
*
      END 
      SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x, 
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N' x := A*x.
*
*              TRANS = 'T' x := A'*x. 
*
*              TRANS = 'C' x := A'*x. 
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U'          A is assumed to be unit triangular.
*
*              DIAG = 'N'          A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(n,1).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. 
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (LDA.GE.N) 
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := A*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         DO 10,I = 1,J - 1
                            X(I) = X(I) + X(J)*A(I,J)
   10                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(J,J)
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX
                  DO 40,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 30,I = 1,J - 1
                            X(IX) = X(IX) + X(JX)*A(I,J)
                            IX = IX + INCX
   30                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     END IF
*
                     JX = JX + INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         DO 50,I = N,J + 1,-1
                            X(I) = X(I) + X(J)*A(I,J)
   50                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(J,J)
                     END IF
*
   60             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 70,I = N,J + 1,-1
                            X(IX) = X(IX) + X(JX)*A(I,J)
                            IX = IX - INCX
   70                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     END IF
*
                     JX = JX - INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := A'*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100,J = N,1,-1
                     IF (NOUNIT) X(J) = X(J)*A(J,J)
                     DO 90,I = J - 1,1,-1
                        X(J) = X(J) + A(I,J)*X(I) 
   90                CONTINUE 
  100             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120,J = N,1,-1
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     DO 110,I = J - 1,1,-1
                        IX = IX - INCX
                        X(JX) = X(JX) + A(I,J)*X(IX)
  110                CONTINUE 
                     JX = JX - INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = 1,N
                     IF (NOUNIT) X(J) = X(J)*A(J,J)
                     DO 130,I = J + 1,N 
                        X(J) = X(J) + A(I,J)*X(I) 
  130                CONTINUE 
  140             CONTINUE
*
              ELSE
                  JX = KX
                  DO 160,J = 1,N
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     DO 150,I = J + 1,N 
                        IX = IX + INCX
                        X(JX) = X(JX) + A(I,J)*X(IX)
  150                CONTINUE 
                     JX = JX + INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STRMV .
*
      END 
      SUBROUTINE STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,K,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STBMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x, 
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U' A is an upper triangular matrix.
*
*              UPLO = 'L' A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N' x := A*x.
*
*              TRANS = 'T' x := A'*x. 
*
*              TRANS = 'C' x := A'*x. 
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U' A is assumed to be unit triangular.
*
*              DIAG = 'N' A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U', K specifies the number of
*           super-diagonals of the matrix A. 
*           On entry with UPLO = 'L', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           Before entry with UPLO = 'L', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the 
*           array A is not referenced.
*           Note that when DIAG = 'U' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. 
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 5-November-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KPLUS1,KX
      INTEGER L
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (K.GE.0) .AND.
     .     (LDA.GE. (K+1))
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX   too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*         Form  x := A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         I = MAX(1,J-K) 
                         DO 10,L = KPLUS1 + I - J,K
                            X(I) = X(I) + X(J)*A(L,J)
                            I = I + 1
   10                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(KPLUS1,J)
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX
                  DO 40,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 30,L = 1 + MAX(KPLUS1-J,0),K
                            X(IX) = X(IX) + X(JX)*A(L,J)
                            IX = IX + INCX
   30                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(KPLUS1,J)
                     END IF
*
                     JX = JX + INCX
                     IF (J.GT.K) KX = KX + INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         I = MIN(N,J+K) 
                         DO 50,L = 1 + I - J,2,-1 
                            X(I) = X(I) + X(J)*A(L,J)
                            I = I - 1
   50                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(1,J)
                     END IF
*
   60             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 70,L = 1 + MIN(K,N-J),2,-1
                            X(IX) = X(IX) + X(JX)*A(L,J)
                            IX = IX - INCX
   70                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(1,J)
                     END IF
*
                     JX = JX - INCX
                     IF ((N-J).GE.K) KX = KX - INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := A'*x.
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100,J = N,1,-1
                     I = J
                     IF (NOUNIT) X(J) = X(J)*A(KPLUS1,J)
                     DO 90,L = K,1 + MAX(KPLUS1-J,0),-1
                        I = I - 1
                        X(J) = X(J) + A(L,J)*X(I) 
   90                CONTINUE 
  100             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 120,J = N,1,-1
                     KX = KX - INCX
                     IX = KX
                     IF (NOUNIT) X(JX) = X(JX)*A(KPLUS1,J)
                     DO 110,L = K,1 + MAX(KPLUS1-J,0),-1
                        X(JX) = X(JX) + A(L,J)*X(IX)
                        IX = IX - INCX
  110                CONTINUE 
                     JX = JX - INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = 1,N
                     I = J
                     IF (NOUNIT) X(J) = X(J)*A(1,J)
                     DO 130,L = 2,1 + MIN(K,N-J)
                        I = I + 1
                        X(J) = X(J) + A(L,J)*X(I) 
  130                CONTINUE 
  140             CONTINUE
*
              ELSE
                  JX = KX
                  DO 160,J = 1,N
                     KX = KX + INCX
                     IX = KX
                     IF (NOUNIT) X(JX) = X(JX)*A(1,J)
                     DO 150,L = 2,1 + MIN(K,N-J)
                        X(JX) = X(JX) + A(L,J)*X(IX)
                        IX = IX + INCX
  150                CONTINUE 
                     JX = JX + INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STBMV .
*
      END 
      SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,INCX
      REAL AP(*),X(*)
*
*  Purpose
*  =======
*
*  STPMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x, 
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N' x := A*x.
*
*              TRANS = 'T' x := A'*x. 
*
*              TRANS = 'C' x := A'*x. 
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U'          A is assumed to be unit triangular.
*
*              DIAG = 'N'          A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
*           respectively, and so on.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
*           respectively, and so on.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. 
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*  Note that UPLO, TRANS, DIAG and N must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 2-October-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,K,KK
      INTEGER KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0)
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of AP are
*     accessed sequentially with one pass through AP.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x:= A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 20,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         DO 10,I = 1,J - 1
                            X(I) = X(I) + X(J)*AP(K)
                            K = K + 1
   10                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*AP(K)
                         K = K + 1
*
                     ELSE
                         K = K + J
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX
                  DO 40,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         KK = K
                         DO 30,K = KK,KK + J - 2
                            X(IX) = X(IX) + X(JX)*AP(K)
                            IX = IX + INCX
   30                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*AP(K)
                         K = K + 1
*
                     ELSE
                         K = K + J
                     END IF
*
                     JX = JX + INCX
   40             CONTINUE
              END IF
*
          ELSE
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 60,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         DO 50,I = N,J + 1,-1
                            X(I) = X(I) + X(J)*AP(K)
                            K = K - 1
   50                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*AP(K)
                         K = K - 1
*
                     ELSE
                         K = K - (N-J+1)
                     END IF
*
   60             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         KK = K
                         DO 70,K = KK,KK - (N- (J+1)),-1
                            X(IX) = X(IX) + X(JX)*AP(K)
                            IX = IX - INCX
   70                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*AP(K)
                         K = K - 1
*
                     ELSE
                         K = K - (N-J+1)
                     END IF
*
                     JX = JX - INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := A'*x.
*
          IF (LSAME(UPLO,'U')) THEN
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 100,J = N,1,-1
                     IF (NOUNIT) X(J) = X(J)*AP(K)
                     K = K - 1
                     DO 90,I = J - 1,1,-1
                        X(J) = X(J) + AP(K)*X(I)
                        K = K - 1
   90                CONTINUE 
  100             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120,J = N,1,-1
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*AP(K)
                     KK = K - 1
                     DO 110,K = KK,KK - J + 2,-1
                        IX = IX - INCX
                        X(JX) = X(JX) + AP(K)*X(IX)
  110                CONTINUE 
                     JX = JX - INCX
  120             CONTINUE
              END IF
*
          ELSE
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 140,J = 1,N
                     IF (NOUNIT) X(J) = X(J)*AP(K)
                     K = K + 1
                     DO 130,I = J + 1,N 
                        X(J) = X(J) + AP(K)*X(I)
                        K = K + 1
  130                CONTINUE 
  140             CONTINUE
*
              ELSE
                  JX = KX
                  DO 160,J = 1,N
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*AP(K)
                     KK = K + 1
                     DO 150,K = KK,KK + N - (J+1) 
                        IX = IX + INCX
                        X(JX) = X(JX) + AP(K)*X(IX)
  150                CONTINUE 
                     JX = JX + INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STPMV .
*
      END 
      SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STRSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows: 
*
*              TRANS = 'N' A*x = b.
*
*              TRANS = 'T' A'*x = b.
*
*              TRANS = 'C' A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U'          A is assumed to be unit triangular.
*
*              DIAG = 'N'          A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(n,1).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (LDA.GE.N) 
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := inv( A )*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(J,J)
                         DO 10,I = J - 1,1,-1
                            X(I) = X(I) - X(J)*A(I,J)
   10                    CONTINUE
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                         IX = JX
                         DO 30,I = J - 1,1,-1
                            IX = IX - INCX
                            X(IX) = X(IX) - X(JX)*A(I,J)
   30                    CONTINUE
                     END IF
*
                     JX = JX - INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(J,J)
                         DO 50,I = J + 1,N
                            X(I) = X(I) - X(J)*A(I,J)
   50                    CONTINUE
                     END IF
*
   60             CONTINUE
*
              ELSE
                  JX = KX
                  DO 80,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                         IX = JX
                         DO 70,I = J + 1,N
                            IX = IX + INCX
                            X(IX) = X(IX) - X(JX)*A(I,J)
   70                    CONTINUE
                     END IF
*
                     JX = JX + INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := inv( A' )*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100,J = 1,N
                     DO 90,I = 1,J - 1
                        X(J) = X(J) - A(I,J)*X(I) 
   90                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(J,J)
  100             CONTINUE
*
              ELSE
                  JX = KX
                  DO 120,J = 1,N
                     IX = KX
                     DO 110,I = 1,J - 1 
                        X(JX) = X(JX) - A(I,J)*X(IX)
                        IX = IX + INCX
  110                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                     JX = JX + INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = N,1,-1
                     DO 130,I = N,J + 1,-1
                        X(J) = X(J) - A(I,J)*X(I) 
  130                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(J,J)
  140             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160,J = N,1,-1
                     IX = KX
                     DO 150,I = N,J + 1,-1
                        X(JX) = X(JX) - A(I,J)*X(IX)
                        IX = IX - INCX
  150                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                     JX = JX - INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STRSV .
*
      END 
      SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,K,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STBSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular band matrix, with ( k + 1 )
*  diagonals.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows: 
*
*              TRANS = 'N' A*x = b.
*
*              TRANS = 'T' A'*x = b.
*
*              TRANS = 'C' A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U' A is assumed to be unit triangular.
*
*              DIAG = 'N' A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U', K specifies the number of
*           super-diagonals of the matrix A. 
*           On entry with UPLO = 'L', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           Before entry with UPLO = 'L', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the 
*           array A is not referenced.
*           Note that when DIAG = 'U' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 7-November-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KPLUS1,KX
      INTEGER L
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (K.GE.0) .AND.
     .     (LDA.GE. (K+1))
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed by sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := inv( A )*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                         I = J
                         DO 10,L = K,1 + MAX(KPLUS1-J,0),-1 
                            I = I - 1
                            X(I) = X(I) - X(J)*A(L,J)
   10                    CONTINUE
                     END IF
*
   20             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40,J = N,1,-1
                     KX = KX - INCX
                     IX = KX
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                         DO 30 L = K,1 + MAX(KPLUS1-J,0),-1 
                            X(IX) = X(IX) - X(JX)*A(L,J)
                            IX = IX - INCX
   30                    CONTINUE
                     END IF
*
                     JX = JX - INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(1,J)
                         I = J
                         DO 50,L = 2,1 + MIN(K,N-J)
                            I = I + 1
                            X(I) = X(I) - X(J)*A(L,J)
   50                    CONTINUE
                     END IF
*
   60             CONTINUE
*
              ELSE
                  JX = KX
                  DO 80,J = 1,N
                     KX = KX + INCX
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                         IX = KX
                         DO 70,L = 2,1 + MIN(K,N-J)
                            X(IX) = X(IX) - X(JX)*A(L,J)
                            IX = IX + INCX
   70                    CONTINUE
                     END IF
*
                     JX = JX + INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := inv( A')*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100,J = 1,N
                     I = MAX(1,J-K)
                     DO 90,L = KPLUS1 + I - J,K
                        X(J) = X(J) - A(L,J)*X(I) 
                        I = I + 1
   90                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
  100             CONTINUE
*
              ELSE
                  JX = KX
                  DO 120,J = 1,N
                     IX = KX
                     DO 110,L = 1 + MAX(KPLUS1-J,0),K
                        X(JX) = X(JX) - A(L,J)*X(IX)
                        IX = IX + INCX
  110                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                     JX = JX + INCX
                     IF (J.GT.K) KX = KX + INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = N,1,-1
                     I = MIN(N,J+K)
                     DO 130,L = 1 + I - J,2,-1
                        X(J) = X(J) - A(L,J)*X(I) 
                        I = I - 1
  130                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(1,J)
  140             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160,J = N,1,-1
                     IX = KX
                     DO 150,L = 1 + MIN(K,N-J),2,-1
                        X(JX) = X(JX) - A(L,J)*X(IX)
                        IX = IX - INCX
  150                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                     JX = JX - INCX
                     IF ((N-J).GE.K) KX = KX - INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STBSV .
*
      END 
      SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,INCX
      REAL AP(*),X(*)
*
*  Purpose
*  =======
*
*  STPSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U' A is an upper triangular matrix.
*
*              UPLO = 'L' A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows: 
*
*              TRANS = 'N' A*x = b.
*
*              TRANS = 'T' A'*x = b.
*
*              TRANS = 'C' A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U' A is assumed to be unit triangular.
*
*              DIAG = 'N' A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
*           respectively, and so on.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
*           respectively, and so on.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 11-November-1985. 
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,K,KK
      INTEGER KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of AP are
*     accessed sequentially with one pass through AP.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := inv( A )*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 20,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/AP(K)
                         K = K - 1
                         DO 10,I = J - 1,1,-1
                            X(I) = X(I) - X(J)*AP(K)
                            K = K - 1
   10                    CONTINUE
*
                     ELSE
                         K = K - J
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/AP(K)
                         IX = JX
                         KK = K - 1
                         DO 30,K = KK,KK - J + 2,-1
                            IX = IX - INCX
                            X(IX) = X(IX) - X(JX)*AP(K)
   30                    CONTINUE
*
                     ELSE
                         K = K - J
                     END IF
*
                     JX = JX - INCX
   40             CONTINUE
              END IF
*
          ELSE
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 60,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/AP(K)
                         K = K + 1
                         DO 50,I = J + 1,N
                            X(I) = X(I) - X(J)*AP(K)
                            K = K + 1
   50                    CONTINUE
*
                     ELSE
                         K = K + N - J + 1
                     END IF
*
   60             CONTINUE
*
              ELSE
                  JX = KX
                  DO 80,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/AP(K)
                         IX = JX
                         KK = K + 1
                         DO 70,K = KK,KK + N - (J+1)
                            IX = IX + INCX
                            X(IX) = X(IX) - X(JX)*AP(K)
   70                    CONTINUE
*
                     ELSE
                         K = K + N - J + 1
                     END IF
*
                     JX = JX + INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := inv( A' )*x.
*
          IF (LSAME(UPLO,'U')) THEN
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 100,J = 1,N
                     DO 90,I = 1,J - 1
                        X(J) = X(J) - AP(K)*X(I)
                        K = K + 1
   90                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/AP(K)
                     K = K + 1
  100             CONTINUE
*
              ELSE
                  JX = KX
                  DO 120,J = 1,N
                     IX = KX
                     KK = K
                     DO 110,K = KK,KK + J - 2
                        X(JX) = X(JX) - AP(K)*X(IX)
                        IX = IX + INCX
  110                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/AP(K)
                     K = K + 1
                     JX = JX + INCX
  120             CONTINUE
              END IF
*
          ELSE
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 140,J = N,1,-1
                     DO 130,I = N,J + 1,-1
                        X(J) = X(J) - AP(K)*X(I)
                        K = K - 1
  130                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/AP(K)
                     K = K - 1
  140             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160,J = N,1,-1
                     IX = KX
                     KK = K
                     DO 150,K = KK,KK - (N- (J+1)),-1
                        X(JX) = X(JX) - AP(K)*X(IX)
                        IX = IX - INCX
  150                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/AP(K)
                     K = K - 1
                     JX = JX - INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STPSV .
*
      END 
      SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      INTEGER M,N,INCX,INCY,LDA
      REAL ALPHA,X(*),Y(*),A(LDA,*)
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix. 
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,m).
*           Unchanged on exit.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,J,JY,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP
      LOGICAL OK
      OK = (M.GT.0) .AND. (N.GT.0) .AND. (LDA.GE.M)
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          DO 20,J = 1,N
             IF (Y(J).NE.ZERO) THEN
                 TEMP = ALPHA*Y(J)
                 DO 10,I = 1,M
                    A(I,J) = A(I,J) + X(I)*TEMP
   10            CONTINUE
             END IF 
*
   20     CONTINUE
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              JY = 1
*
          ELSE
              JY = 1 - (N-1)*INCY
          END IF
*
          DO 40,J = 1,N
             IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                 IX = KX
                 DO 30,I = 1,M
                    A(I,J) = A(I,J) + X(IX)*TEMP
                    IX = IX + INCX
   30            CONTINUE
             END IF 
*
             JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of SGER  .
*
      END 
      SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
      CHARACTER *1 UPLO
      INTEGER N,INCX,LDA
      REAL ALPHA,X(*),A(LDA,*)
*
*  Purpose
*  =======
*
*  SSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows: 
*
*              UPLO = 'U' Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,n).
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,J,JX,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (LDA.GE.N)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in upper triangle.
*
          IF (INCX.EQ.1) THEN 
              DO 20,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 10,I = 1,J
                        A(I,J) = A(I,J) + X(I)*TEMP
   10                CONTINUE 
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              DO 40,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = KX
                     DO 30,I = 1,J
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
   30                CONTINUE 
                 END IF
*
                 JX = JX + INCX
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
          IF (INCX.EQ.1) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 50,I = J,N
                        A(I,J) = A(I,J) + X(I)*TEMP
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = JX
                     DO 70,I = J,N
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSYR  .
*
      END 
      SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
      CHARACTER *1 UPLO
      INTEGER N,INCX
      REAL ALPHA,X(*),AP(*)
*
*  Purpose
*  =======
*
*  SSPR    performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows: 
*
*              UPLO = 'U' The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) 
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) 
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,J,JX,K,KK
      INTEGER KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      K = 1
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when upper triangle is stored in AP.
*
          IF (INCX.EQ.1) THEN 
              DO 20,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 10,I = 1,J
                        AP(K) = AP(K) + X(I)*TEMP 
                        K = K + 1
   10                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              DO 40,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = KX
                     KK = K
                     DO 30,K = KK,KK + J - 1
                        AP(K) = AP(K) + X(IX)*TEMP
                        IX = IX + INCX
   30                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
                 JX = JX + INCX
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
          IF (INCX.EQ.1) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 50,I = J,N
                        AP(K) = AP(K) + X(I)*TEMP 
                        K = K + 1
   50                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = JX
                     KK = K
                     DO 70,K = KK,KK + N - J
                        AP(K) = AP(K) + X(IX)*TEMP
                        IX = IX + INCX
   70                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
                 JX = JX + INCX
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSPR  .
*
      END 
      SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      CHARACTER *1 UPLO
      INTEGER N,INCX,INCY,LDA 
      REAL ALPHA,X(*),Y(*),A(LDA,*)
*
*  Purpose
*  =======
*
*  SSYR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an n
*  by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows: 
*
*              UPLO = 'U' Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,n).
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KX,KY 
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (LDA.GE.N)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set up the start points in X and Y if the increments are not both 
*     unity.
*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 20,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 10,I = 1,J
                        A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                CONTINUE 
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 40,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = KX
                     IY = KY
                     DO 30,I = 1,J
                        A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   30                CONTINUE 
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when A is stored in the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 50,I = J,N
                        A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 80,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = JX
                     IY = JY
                     DO 70,I = J,N
                        A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSYR2 .
*
      END 
      SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      CHARACTER *1 UPLO
      INTEGER N,INCX,INCY
      REAL ALPHA,X(*),Y(*),AP(*)
*
*  Purpose
*  =======
*
*  SSPR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows: 
*
*              UPLO = 'U' The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) 
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) 
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER K,KK,KX,KY
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set up the start points in X and Y if the increments are not both 
*     unity.
*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      K = 1
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when upper triangle is stored in AP.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 20,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 10,I = 1,J
                        AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                        K = K + 1
   10                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 40,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = KX
                     IY = KY
                     KK = K
                     DO 30,K = KK,KK + J - 1
                        AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   30                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 50,I = J,N
                        AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                        K = K + 1
   50                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 80,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = JX
                     IY = JY
                     KK = K
                     DO 70,K = KK,KK + N - J
                        AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   70                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSPR2 .
*
      END 
      
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module ZDSCAL from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 14:48:38 1998.
* ======================================================================
      subroutine  zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*)
      double precision da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end

* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module ZDOTC from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 15:01:53 1998.
* ======================================================================
      double complex function zdotc(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end

* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module ZAXPY from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 15:02:56 1998.
* ======================================================================
      subroutine zaxpy(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),za
      integer i,incx,incy,ix,iy,n
      double precision dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end

* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module ZHPR2 from package BLAS2.
* Retrieved from NETLIB on Mon Jun 15 15:03:49 1998.
* ======================================================================
      SUBROUTINE ZHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX*16         AP( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHPR2  performs the hermitian rank 2 operation
*
*     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an
*  n by n hermitian matrix, supplied in packed form.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  AP     - COMPLEX*16       array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set, they are assumed to be zero, and on exit they
*           are set to zero.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DBLE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZHPR2 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when upper triangle is stored in AP.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( J ) )
                  TEMP2 = DCONJG( ALPHA*X( J ) )
                  K     = KK
                  DO 10, I = 1, J - 1
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   10             CONTINUE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) ) +
     $                               DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) )
               END IF
               KK = KK + J
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( JY ) )
                  TEMP2 = DCONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  DO 30, K = KK, KK + J - 2
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   30             CONTINUE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) ) +
     $                               DBLE( X( JX )*TEMP1 +
     $                                     Y( JY )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1   = ALPHA*DCONJG( Y( J ) )
                  TEMP2   = DCONJG( ALPHA*X( J ) )
                  AP( KK ) = DBLE( AP( KK ) ) +
     $                       DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
                  K        = KK               + 1
                  DO 50, I = J + 1, N
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   50             CONTINUE
               ELSE
                  AP( KK ) = DBLE( AP( KK ) )
               END IF
               KK = KK + N - J + 1
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1    = ALPHA*DCONJG( Y( JY ) )
                  TEMP2    = DCONJG( ALPHA*X( JX ) )
                  AP( KK ) = DBLE( AP( KK ) ) +
     $                       DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX       = JX
                  IY       = JY
                  DO 70, K = KK + 1, KK + N - J
                     IX      = IX      + INCX
                     IY      = IY      + INCY
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  AP( KK ) = DBLE( AP( KK ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + N - J + 1
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZHPR2 .
*
      END
      
* ======================================================================
* NIST Guide to Available Math Software.
* Source for module ZGEMV from package BLAS2.
* Retrieved from NETLIB on Mon Jun 15 15:22:37 1998.
* ======================================================================
      SUBROUTINE ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
*
*     y := alpha*conjg( A' )*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
      NOCONJ = LSAME( TRANS, 'T' )
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZGEMV .
*
      END
      
* ======================================================================
* NIST Guide to Available Math Software.
* Source for module ZGERC from package BLAS2.
* Retrieved from NETLIB on Mon Jun 15 15:24:20 1998.
* ======================================================================
      SUBROUTINE ZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGERC  performs the rank 1 operation
*
*     A := alpha*x*conjg( y' ) + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGERC ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*DCONJG( Y( JY ) )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*DCONJG( Y( JY ) )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERC .
*
      END

* ======================================================================
* NIST Guide to Available Math Software.
* Source for module DZNRM2 from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 15:24:54 1998.
* ======================================================================
      DOUBLE PRECISION FUNCTION DZNRM2( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      COMPLEX*16                        X( * )
*     ..
*
*  DZNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DZNRM2 := sqrt( conjg( x' )*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to ZLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      NORM, SCALE, SSQ, TEMP
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, DIMAG, DBLE, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DZNRM2 = NORM
      RETURN
*
*     End of DZNRM2.
*
      END
      
* ======================================================================
* NIST Guide to Available Math Software.
* Source for module ZSCAL from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 15:25:27 1998.
* ======================================================================
      subroutine  zscal(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex za,zx(*)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end

* ======================================================================
* NIST Guide to Available Math Software.
* Source for module ZSWAP from package BLAS1.
* Retrieved from NETLIB on Mon Jun 15 15:25:56 1998.
* ======================================================================
      subroutine  zswap (n,zx,incx,zy,incy)
c
c     interchanges two vectors.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
