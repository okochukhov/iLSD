! =============================================================================================
!
!   iLSD: improved least-squares deconvolution of stellar intensity and polarization spectra
!
!   Developed by Oleg Kochukhov (Uppsala University), oleg.kochukhov@physics.uu.se
!
!   Users of this code should cite Kochukhov et al. (2010, A&A, 524, A5)
!
! =============================================================================================
!
! ilsd.f90: Fortran 90 version of the least-squares deconvolution code.  The observed
! spectrum is described as a superposition of identical profiles, scaled by the weights of
! a line mask and shifted to the position of every line, so that the model spectrum is the
! product of a line-pattern matrix M and a mean profile Z.  The code fills M for a given
! line mask and velocity grid and solves the corresponding linear least-squares problem,
! Z = (M^T S^2 M + Lambda R)^-1 M^T S^2 Y, where S is the diagonal matrix of inverse error
! bars and R is the first-order Tikhonov regularization matrix.  Several mean profiles can
! be reconstructed simultaneously (multiprofile LSD) by concatenating the line-pattern
! matrices built from different columns of line weights.
!
! Required input data:
!
!  1. Configuration file with one item per line
!       'file.obs'           = name of the file with the observed spectrum
!       'file.lin'           = name of the file with the line mask
!       Vfirst, Vlast, Vstep = velocity limits and step of the LSD profile grid (km / s)
!       Pol                  = 0 for intensity, 1 for polarization, 2 for log(intensity)
!     and two optional lines
!       Reg                  = Tikhonov regularization parameter, negative for the default
!                              smoothing of the reconstructed profiles, 0 or absent for none
!       Iout                 = amount of output; each level adds files to the previous one:
!                              0 or absent for .lsd alone, 1 to add .wln and .mod, 2 to add .cov
!
!  2. Observed spectrum: 2 or 3 column table of wavelength (A), continuum normalized
!     intensity or polarization, and optional error bar.  A single default error bar is
!     assumed when the third column is absent.
!
!  3. Line mask: table of the line central wavelength (A) followed by one or more columns
!     of line weights.  Each column of weights yields one LSD profile.
!
! Output data (formatted text files, <prefix> is the configuration file name without its
! extension; which of them are written is controlled by Iout, except for .cor, which is
! always written when the line weights are iterated):
!   <prefix>.lsd  = velocity, LSD profile and its error bar, one block per profile
!   <prefix>.wln  = line wavelength and total observational weight of every mask line
!   <prefix>.mod  = wavelength, LSD model, observed spectrum and its error bar
!   <prefix>.cov  = covariance matrix of the LSD profiles
!   <prefix>.cor  = accumulated line weight corrections, written when Niter > 0
!
! Usage:
!   1. Reconstruct the LSD profiles:
!   ilsd <configuration_file>
!
!   2. Reconstruct the LSD profiles, adjusting the line weights in Niter iterations:
!   ilsd <configuration_file> <Niter>
!
! =============================================================================================

! Observed spectrum, line mask, LSD velocity grid and run configuration
MODULE LsdData

  IMPLICIT NONE

  REAL(8), PARAMETER :: Vc = 299792.458d0          ! speed of light (km / s)

  ! observed spectrum
  INTEGER :: Nwl = 0                               ! number of spectral points
  REAL(8), ALLOCATABLE :: Wl(:)                    ! wavelengths (A)
  REAL(8), ALLOCATABLE :: Sp(:)                    ! observed intensity or polarization
  REAL(8), ALLOCATABLE :: SigSp(:)                 ! inverse error bars, 1 / sigma

  ! line mask
  INTEGER :: Nlin = 0                              ! number of lines
  REAL(8), ALLOCATABLE :: WlcLin(:)                ! line central wavelengths (A)
  REAL(8), ALLOCATABLE :: WgtLin(:,:)              ! line weights (Nlsd, Nlin)
  REAL(8), ALLOCATABLE :: WgtIni(:,:)              ! line weights of the input mask

  ! velocity grid of the LSD profiles
  INTEGER :: Nlsd = 0                              ! number of LSD profiles
  INTEGER :: Nvlsd = 0                             ! number of velocity bins per profile
  REAL(8) :: Vfirst = 0d0                          ! first velocity of the grid (km / s)
  REAL(8) :: Vlast = 0d0                           ! last velocity of the grid (km / s)
  REAL(8) :: Vstep = 0d0                           ! velocity step of the grid (km / s)
  REAL(8), ALLOCATABLE :: Vlsd(:)                  ! velocity grid (km / s)

  ! run configuration
  INTEGER :: Pol = 0                               ! 0 = intensity, 1 = polarization, 2 = log
  INTEGER :: Niter = 0                             ! number of line weight iterations
  INTEGER :: Iout = 0                              ! amount of output, see ReportOutput
  REAL(8) :: Reg = 0d0                             ! Tikhonov regularization parameter

  ! output levels selected by Iout; every level adds files to those of the previous one
  INTEGER, PARAMETER :: Iolsd = 0                  ! .lsd alone
  INTEGER, PARAMETER :: Iomod = 1                  ! adds .wln and .mod
  INTEGER, PARAMETER :: Iocov = 2                  ! adds .cov

  ! calculation options fixed at compile time
  LOGICAL, PARAMETER :: Lscl = .TRUE.              ! scale LSD error bars by sqrt(chi2)
  LOGICAL, PARAMETER :: Verify = .FALSE.           ! check the accuracy of the inversion

  ! line weight adjustement parameters
  REAL(8), PARAMETER :: Wmu = 0.1d0                ! prior on the adjusted line weights
  REAL(8), PARAMETER :: Chitol = 1.05d0            ! tolerated rise of the reduced chi-square

END MODULE LsdData

! Error handling, string manipulation and input of free-format numerical tables
MODULE LsdUtils

  IMPLICIT NONE

  INTEGER, PARAMETER :: Ucfg = 1                   ! unit of the configuration file
  INTEGER, PARAMETER :: Udat = 2                   ! unit of the input data files
  INTEGER, PARAMETER :: Uout = 11                  ! unit of the output files
  INTEGER, PARAMETER :: Ucov = 12                  ! unit of the covariance matrix file

  CONTAINS

  ! Report an error message and terminate the program
  SUBROUTINE ErrStop(Msg)

    CHARACTER(*), INTENT(IN) :: Msg

    WRITE(*, '(A)') ' *** ERROR: ' // TRIM(Msg)
    STOP 1

  END SUBROUTINE ErrStop

  ! Remove all blanks from a character string, padding the result with blanks on the right
  FUNCTION Compress(Name) RESULT(Out)

    CHARACTER(*), INTENT(IN) :: Name
    CHARACTER(LEN(Name)) :: Out
    INTEGER :: i, j

    Out = ' '
    j = 0
    DO i = 1, LEN(Name)
      IF (Name(i:i) == ' ') CYCLE
      j = j + 1
      Out(j:j) = Name(i:i)
    END DO

  END FUNCTION Compress

  ! Count the blank- or tab-separated tokens of a character string
  INTEGER FUNCTION Ntokens(Line)

    CHARACTER(*), INTENT(IN) :: Line
    INTEGER :: i
    LOGICAL :: Inside

    Ntokens = 0
    Inside = .FALSE.
    DO i = 1, LEN(Line)
      IF (Line(i:i) == ' ' .OR. IACHAR(Line(i:i)) == 9) THEN
        Inside = .FALSE.
      ELSE IF (.NOT. Inside) THEN
        Inside = .TRUE.
        Ntokens = Ntokens + 1
      END IF
    END DO

  END FUNCTION Ntokens

  ! Keep only the elements listed in Idx of an allocatable array, shrinking it to their number
  SUBROUTINE Shrink(A, Idx)

    REAL(8), ALLOCATABLE, INTENT(INOUT) :: A(:)
    INTEGER, INTENT(IN) :: Idx(:)
    REAL(8), ALLOCATABLE :: Tmp(:)

    ALLOCATE(Tmp(SIZE(Idx)))
    Tmp = A(Idx)
    DEALLOCATE(A)
    ALLOCATE(A(SIZE(Tmp)))
    A = Tmp
    DEALLOCATE(Tmp)

  END SUBROUTINE Shrink

  ! Read a free-format numerical table, allocating the output array to the data size
  SUBROUTINE ReadColumns(Dfile, Dat, Nx, Ny)

    CHARACTER(*), INTENT(IN) :: Dfile
    REAL(8), ALLOCATABLE, INTENT(OUT) :: Dat(:,:)
    INTEGER, INTENT(OUT) :: Nx, Ny
    INTEGER :: j, Ios
    REAL(8), ALLOCATABLE :: Row(:)
    CHARACTER(4096) :: Sline

    ! open the file
    OPEN(UNIT=Udat, FILE=Dfile, STATUS='OLD', ACTION='READ', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadColumns: cannot open input data file: ' // TRIM(Dfile))

    ! the number of columns is given by the first non-empty record
    DO
      READ(Udat, '(A)', IOSTAT=Ios) Sline
      IF (Ios /= 0) CALL ErrStop('ReadColumns: no data found in file: ' // TRIM(Dfile))
      IF (LEN_TRIM(Sline) > 0) EXIT
    END DO
    Nx = Ntokens(Sline)
    IF (Nx < 1) CALL ErrStop('ReadColumns: no columns found in file: ' // TRIM(Dfile))

    ! count the rows by reading the whole file the same way as it is read below
    ALLOCATE(Row(Nx))
    REWIND(Udat)
    Ny = 0
    DO
      READ(Udat, *, IOSTAT=Ios) (Row(j), j = 1, Nx)
      IF (Ios > 0) CALL ErrStop('ReadColumns: error reading file: ' // TRIM(Dfile))
      IF (Ios /= 0) EXIT
      Ny = Ny + 1
    END DO
    DEALLOCATE(Row)
    IF (Ny < 1) CALL ErrStop('ReadColumns: no data rows found in file: ' // TRIM(Dfile))

    ! read the data
    ALLOCATE(Dat(Nx,Ny), STAT=Ios)
    IF (Ios /= 0) CALL ErrStop('ReadColumns: cannot allocate data array for ' // TRIM(Dfile))
    REWIND(Udat)
    DO j = 1, Ny
      READ(Udat, *, IOSTAT=Ios) Dat(:,j)
      IF (Ios /= 0) CALL ErrStop('ReadColumns: error reading file: ' // TRIM(Dfile))
    END DO
    CLOSE(Udat)

  END SUBROUTINE ReadColumns

END MODULE LsdUtils

! Line-pattern matrix in compressed sparse row form and the three products it enters
MODULE Sparse

  USE LsdUtils

  IMPLICIT NONE

  ! Sparse matrix in compressed sparse row form: row r holds the elements Val(k) in columns
  ! Col(k) for k = Ptr(r) ... Ptr(r+1)-1.  The line-pattern matrix is very sparse, because
  ! every line contributes to only two velocity bins of every LSD profile at any spectral
  ! point, and only a handful of lines reach a given point.  A column may appear more than
  ! once in a row when two lines land in the same velocity bin; the row is then the sum of
  ! those elements, which is exactly what the three products below compute, so duplicates
  ! need not be merged.
  TYPE SparseType
    INTEGER :: Nrow = 0                            ! rows, i.e. spectral points
    INTEGER :: Ncol = 0                            ! columns, i.e. LSD profile points
    INTEGER :: Nnz = 0                             ! stored elements
    INTEGER, ALLOCATABLE :: Ptr(:)                 ! first element of every row, size Nrow+1
    INTEGER, ALLOCATABLE :: Col(:)                 ! column index of every element
    REAL(8), ALLOCATABLE :: Val(:)                 ! value of every element
  END TYPE SparseType

  CONTAINS

  ! Allocate a sparse matrix with Nrow rows, Ncol columns and room for Nnz elements
  SUBROUTINE SparseInit(M, Nrow, Ncol, Nnz)

    TYPE(SparseType), INTENT(OUT) :: M
    INTEGER, INTENT(IN) :: Nrow, Ncol, Nnz
    INTEGER :: Ios

    M%Nrow = Nrow
    M%Ncol = Ncol
    M%Nnz = Nnz
    ALLOCATE(M%Ptr(Nrow+1), M%Col(Nnz), M%Val(Nnz), STAT=Ios)
    IF (Ios /= 0) CALL ErrStop('SparseInit: cannot allocate the sparse matrix')
    M%Ptr = 1

  END SUBROUTINE SparseInit

  ! Release the storage of a sparse matrix
  SUBROUTINE SparseFree(M)

    TYPE(SparseType), INTENT(INOUT) :: M

    IF (ALLOCATED(M%Ptr)) DEALLOCATE(M%Ptr)
    IF (ALLOCATED(M%Col)) DEALLOCATE(M%Col)
    IF (ALLOCATED(M%Val)) DEALLOCATE(M%Val)
    M%Nrow = 0
    M%Ncol = 0
    M%Nnz = 0

  END SUBROUTINE SparseFree

  ! Autocorrelation matrix Am = M^T # S^2 # M, filled in both triangles.  Every spectral
  ! point contributes the outer product of the non-zero elements of its own row of M,
  ! scaled by the inverse variance of that point.  The dense equivalent is
  !   DO p = 1, Ncol
  !     Am(p,p:Ncol) = MATMUL(M(:,p) * S2, M(:,p:Ncol))
  !   END DO
  ! followed by filling the lower triangle by symmetry.
  SUBROUTINE SparseAutoCorr(M, S2, Am)

    TYPE(SparseType), INTENT(IN) :: M
    REAL(8), INTENT(IN) :: S2(M%Nrow)
    REAL(8), INTENT(OUT) :: Am(M%Ncol,M%Ncol)
    INTEGER :: r, k, l, p
    REAL(8) :: Vk

    Am = 0d0
    DO r = 1, M%Nrow
      DO k = M%Ptr(r), M%Ptr(r+1) - 1
        p = M%Col(k)
        Vk = M%Val(k) * S2(r)
        DO l = M%Ptr(r), M%Ptr(r+1) - 1
          Am(M%Col(l),p) = Am(M%Col(l),p) + Vk * M%Val(l)
        END DO
      END DO
    END DO

  END SUBROUTINE SparseAutoCorr

  ! Weighted cross-correlation vector Rhs = M^T # S^2 # Y.  The dense equivalent is
  !   Rhs = MATMUL(Y * S2, M)
  SUBROUTINE SparseCrossCorr(M, S2, Y, Rhs)

    TYPE(SparseType), INTENT(IN) :: M
    REAL(8), INTENT(IN) :: S2(M%Nrow), Y(M%Nrow)
    REAL(8), INTENT(OUT) :: Rhs(M%Ncol)
    INTEGER :: r, k
    REAL(8) :: Wr

    Rhs = 0d0
    DO r = 1, M%Nrow
      Wr = Y(r) * S2(r)
      DO k = M%Ptr(r), M%Ptr(r+1) - 1
        Rhs(M%Col(k)) = Rhs(M%Col(k)) + Wr * M%Val(k)
      END DO
    END DO

  END SUBROUTINE SparseCrossCorr

  ! Model spectrum Y = M # Z.  The dense equivalent is
  !   Y = MATMUL(M, Z)
  SUBROUTINE SparseModel(M, Z, Y)

    TYPE(SparseType), INTENT(IN) :: M
    REAL(8), INTENT(IN) :: Z(M%Ncol)
    REAL(8), INTENT(OUT) :: Y(M%Nrow)
    INTEGER :: r, k
    REAL(8) :: Sum1

    DO r = 1, M%Nrow
      Sum1 = 0d0
      DO k = M%Ptr(r), M%Ptr(r+1) - 1
        Sum1 = Sum1 + M%Val(k) * Z(M%Col(k))
      END DO
      Y(r) = Sum1
    END DO

  END SUBROUTINE SparseModel

END MODULE Sparse

! Least-squares deconvolution: input, line-pattern matrix, inversion, weight adjustment
MODULE LsdCalc

  USE LsdData
  USE LsdUtils
  USE Sparse

  IMPLICIT NONE

  ! inversion of a symmetric positive definite matrix, provided by lineq4.f
  INTERFACE
    SUBROUTINE LINEQ(Am, Ai, N, Info)
      INTEGER :: N, Info
      REAL(8) :: Am(N,N), Ai(N,N)
    END SUBROUTINE LINEQ
  END INTERFACE

  CONTAINS

  ! Read the configuration file, the observed spectrum and the line mask
  SUBROUTINE ReadInput(Prefix)

    CHARACTER(*), INTENT(OUT) :: Prefix
    INTEGER :: iv, Nx, Ny, Ios, Idot, Islash
    REAL(8), ALLOCATABLE :: Dat(:,:)
    CHARACTER(1024) :: Stream, Spfile, Linfile, Arg

    ! configuration file name from the command line or from the terminal
    Stream = ' '
    IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, Stream)
    DO
      IF (Stream == ' ') THEN
        WRITE(*, '(A)', ADVANCE='NO') ' Enter input filename ? ... '
        READ(*, '(A)', IOSTAT=Ios) Stream
        IF (Ios /= 0) CALL ErrStop('ReadInput: cannot read the input file name')
      END IF
      Stream = Compress(Stream)
      OPEN(UNIT=Ucfg, FILE=Stream, STATUS='OLD', ACTION='READ', IOSTAT=Ios)
      IF (Ios == 0) EXIT
      WRITE(*, '(A)') ' Cannot open input data file: ' // TRIM(Stream)
      Stream = ' '
    END DO

    ! prefix of the output files: configuration file name without its extension
    Islash = INDEX(Stream, '/', BACK=.TRUE.)
    Idot = INDEX(Stream, '.', BACK=.TRUE.)
    IF (Idot > Islash + 1) THEN
      Prefix = Stream(1:Idot-1)
    ELSE
      Prefix = TRIM(Stream)
    END IF

    ! optional number of line weight adjustment iterations
    Niter = 0
    IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
      CALL GET_COMMAND_ARGUMENT(2, Arg)
      IF (Arg /= ' ') THEN
        READ(Arg, *, IOSTAT=Ios) Niter
        IF (Ios /= 0) CALL ErrStop('ReadInput: error reading iteration number: ' // TRIM(Arg))
        Niter = MAX(Niter, 0)
      END IF
    END IF
    IF (Niter > 0) WRITE(*, '(A,I0,A,T75,A)') ' Will make ', Niter, ' iterations', &
                                              '*** ReadInput'

    ! file name of the observed spectrum
    READ(Ucfg, *, IOSTAT=Ios) Spfile
    IF (Ios /= 0) CALL ErrStop('ReadInput: error reading filename with observed spectrum')
    Spfile = Compress(Spfile)

    ! observed spectrum, with the error bar inverted
    CALL ReadColumns(Spfile, Dat, Nx, Ny)
    IF (Nx < 2) CALL ErrStop('ReadInput: less than two columns in ' // TRIM(Spfile))
    Nwl = Ny
    ALLOCATE(Wl(Nwl), Sp(Nwl), SigSp(Nwl))
    Wl = Dat(1,:)
    Sp = Dat(2,:)
    IF (Nx == 3) THEN
      IF (ANY(Dat(3,:) <= 0d0)) &
        CALL ErrStop('ReadInput: non-positive error bar in ' // TRIM(Spfile))
      SigSp = 1d0 / Dat(3,:)
    ELSE
      SigSp = 1d2
    END IF
    DEALLOCATE(Dat)
    IF (ANY(Wl <= 0d0)) CALL ErrStop('ReadInput: non-positive wavelength in ' // TRIM(Spfile))
    WRITE(*, '(A,I0,A,T75,A)') ' ' // TRIM(Spfile) // ': ', Nwl, ' wavelength points', &
                               '*** ReadInput'

    ! file name of the line mask
    READ(Ucfg, *, IOSTAT=Ios) Linfile
    IF (Ios /= 0) CALL ErrStop('ReadInput: error reading filename with line list')
    Linfile = Compress(Linfile)

    ! line mask: one column of central wavelengths and one column of weights per profile
    CALL ReadColumns(Linfile, Dat, Nx, Ny)
    IF (Nx < 2) CALL ErrStop('ReadInput: less than two columns in ' // TRIM(Linfile))
    Nlin = Ny
    Nlsd = Nx - 1
    ALLOCATE(WlcLin(Nlin), WgtLin(Nlsd,Nlin))
    WlcLin = Dat(1,:)
    WgtLin = Dat(2:Nx,:)
    DEALLOCATE(Dat)
    IF (ANY(WlcLin <= 0d0)) &
      CALL ErrStop('ReadInput: non-positive line wavelength in ' // TRIM(Linfile))
    WRITE(*, '(A,I0,A,I0,A,T75,A)') ' ' // TRIM(Linfile) // ': ', Nlin, ' x ', Nlsd, &
                                    ' line weights', '*** ReadInput'

    ! velocity limits and step of the LSD profile grid
    READ(Ucfg, *, IOSTAT=Ios) Vfirst, Vlast, Vstep
    IF (Ios /= 0) CALL ErrStop('ReadInput: error in formatting velocity limits and/or step')
    IF (Vstep <= 0d0) CALL ErrStop('ReadInput: velocity step must be positive')
    IF (Vlast - Vfirst < Vstep) &
      CALL ErrStop('ReadInput: velocity range is smaller than the velocity step')

    ! polarization flag
    READ(Ucfg, *, IOSTAT=Ios) Pol
    IF (Ios /= 0) CALL ErrStop('ReadInput: error in formatting polarization flag')

    ! bring the observed spectrum and the line weights to the scale implied by the flag
    SELECT CASE (Pol)

    CASE (0)
      WRITE(*, '(A,T75,A)') ' Input interpreted as intensity', '*** ReadInput'
      Sp = 1d0 - Sp

    CASE (1)
      WRITE(*, '(A,T75,A)') ' Input interpreted as polarization', '*** ReadInput'

    CASE (2)
      WRITE(*, '(A,T75,A)') ' Input converted to log(intensity)', '*** ReadInput'
      IF (ANY(Sp <= 0d0)) &
        CALL ErrStop('ReadInput: non-positive intensity cannot be converted to log scale')
      IF (ANY(WgtLin >= 1d0)) &
        CALL ErrStop('ReadInput: line weight >= 1 cannot be converted to log scale')
      SigSp = SigSp * Sp
      Sp = LOG(Sp)
      WgtLin = -LOG(1d0 - WgtLin)

    CASE DEFAULT
      CALL ErrStop('ReadInput: polarization flag must be 0, 1 or 2')

    END SELECT

    ! keep the weights of the input mask, on the scale just established, as the reference
    ! that the prior of AdjustLsd pulls the adjusted weights towards
    ALLOCATE(WgtIni(Nlsd,Nlin))
    WgtIni = WgtLin

    ! equidistant velocity grid of the LSD profiles
    Nvlsd = NINT((Vlast - Vfirst) / Vstep) + 1
    IF (Nvlsd < 2) CALL ErrStop('ReadInput: less than two velocity points in LSD profile')
    ALLOCATE(Vlsd(Nvlsd))
    Vlsd = (/ (Vfirst + Vstep * (iv - 1), iv = 1, Nvlsd) /)
    Vlast = Vlsd(Nvlsd)
    WRITE(*, '(A,I0,A,F0.1,A,F0.1,A,T75,A)') ' ', Nvlsd, ' velocity points in LSD profile ' &
      // 'from ', Vlsd(1), ' to ', Vlsd(Nvlsd), ' km/s', '*** ReadInput'

    ! optional regularization parameter and, after it, an optional flag selecting how much
    ! output is written; every level adds files to the ones written by the previous level
    Reg = 0d0
    Iout = 0
    READ(Ucfg, *, IOSTAT=Ios) Reg
    IF (Ios /= 0) THEN
      Reg = 0d0
    ELSE
      IF (Reg < 0d0) THEN
        WRITE(*, '(A,T75,A)') ' Default regularization', '*** ReadInput'
      ELSE
        WRITE(*, '(A,ES8.2,T75,A)') ' Regularization parameter = ', Reg, '*** ReadInput'
      END IF
      READ(Ucfg, *, IOSTAT=Ios) Iout
      IF (Ios /= 0) Iout = 0
      Iout = MAX(Iout, 0)
    END IF

    ! close the configuration file
    CLOSE(Ucfg)

  END SUBROUTINE ReadInput

  ! Remove the spectral pixels that no line of the mask contributes to
  SUBROUTINE PreInit

    INTEGER :: iwl, iline, Nwluse
    INTEGER, ALLOCATABLE :: Iuse(:)
    REAL(8) :: Dv

    ! collect the indices of the pixels covered by at least one line
    ALLOCATE(Iuse(Nwl))
    Nwluse = 0
    DO iwl = 1, Nwl
      DO iline = 1, Nlin
        Dv = (Wl(iwl) - WlcLin(iline)) / WlcLin(iline) * Vc
        IF (Dv > Vfirst .AND. Dv < Vlast) THEN
          Nwluse = Nwluse + 1
          Iuse(Nwluse) = iwl
          EXIT
        END IF
      END DO
    END DO
    IF (Nwluse < 1) CALL ErrStop('PreInit: no spectral point is covered by the line mask')

    ! keep only the used wavelength points, shrinking the spectrum arrays to their number
    CALL Shrink(Wl, Iuse(1:Nwluse))
    CALL Shrink(Sp, Iuse(1:Nwluse))
    CALL Shrink(SigSp, Iuse(1:Nwluse))
    DEALLOCATE(Iuse)

    ! report the number of used and removed wavelength points
    WRITE(*, '(A,I0,A,I0,A,T75,A)') ' ', Nwluse, ' points are used, ', Nwl - Nwluse, &
                                    ' are removed', '*** PreInit'
    Nwl = Nwluse

  END SUBROUTINE PreInit

  ! Compute the line-pattern matrix and the total observational weight of every line
  SUBROUTINE LineMatrix(Mm, Luse, Msig)

    TYPE(SparseType), INTENT(OUT) :: Mm
    INTEGER, INTENT(OUT) :: Luse(3,Nlin)
    REAL(8), INTENT(OUT) :: Msig(Nlin)
    INTEGER :: iwl, iline, ilsd, iv, ii, k, Nnz, Nlinuse
    REAL(8) :: Dv, Fv, Ssig, Wgt, Wtot, Rnnz
    REAL(8) :: Mwgt(Nlsd)

    ! line contribution counters and mean line 1/sigma^2 accumulator
    Luse = 0
    Msig = 0d0

    ! first pass: loop through the wavelength points
    DO iwl = 1, Nwl
      Ssig = SigSp(iwl) * SigSp(iwl)

      ! loop through the line list
      DO iline = 1, Nlin

        ! consider the line if it falls within the LSD velocity grid from the current point
        Dv = (Wl(iwl) - WlcLin(iline)) / WlcLin(iline) * Vc
        IF (Dv <= Vfirst .OR. Dv >= Vlast) CYCLE

        ! first and last spectral pixel this line contributes to, and 1/sigma^2 accumulator
        IF (Luse(1,iline) == 0) Luse(2,iline) = iwl
        Luse(3,iline) = iwl
        Luse(1,iline) = Luse(1,iline) + 1
        Msig(iline) = Msig(iline) + Ssig
      END DO
    END DO

    ! every line stores two elements per LSD profile at every point it contributes to
    Rnnz = 2d0 * Nlsd * SUM(REAL(Luse(1,:), 8))
    IF (Rnnz > REAL(HUGE(Nnz), 8)) &
      CALL ErrStop('LineMatrix: too many elements in the line-pattern matrix')
    Nnz = NINT(Rnnz)
    CALL SparseInit(Mm, Nwl, Nlsd * Nvlsd, Nnz)

    ! second pass: store the line-pattern matrix row by row, in order of increasing row
    k = 0
    DO iwl = 1, Nwl
      Mm%Ptr(iwl) = k + 1
      DO iline = 1, Nlin
        Dv = (Wl(iwl) - WlcLin(iline)) / WlcLin(iline) * Vc
        IF (Dv <= Vfirst .OR. Dv >= Vlast) CYCLE

        ! velocity bin below the current point and its linear interpolation weight; the
        ! index is clipped so that rounding cannot push it outside the profile block
        iv = MIN(MAX(INT((Dv - Vfirst) / Vstep) + 1, 1), Nvlsd - 1)
        Fv = (Dv - Vlsd(iv)) / Vstep

        ! store the bidiagonal contribution of this line to every LSD profile block
        DO ilsd = 1, Nlsd
          ii = iv + (ilsd - 1) * Nvlsd
          Mm%Col(k+1) = ii
          Mm%Val(k+1) = WgtLin(ilsd,iline) * (1d0 - Fv)
          Mm%Col(k+2) = ii + 1
          Mm%Val(k+2) = WgtLin(ilsd,iline) * Fv
          k = k + 2
        END DO
      END DO
    END DO
    Mm%Ptr(Nwl+1) = k + 1

    ! count the lines used in the line-pattern matrix and accumulate the weighted mean
    ! line weight of every LSD profile
    Nlinuse = COUNT(Luse(1,:) > 0)
    Wtot = 0d0
    Mwgt = 0d0
    DO iline = 1, Nlin
      IF (Luse(1,iline) <= 0) CYCLE
      Wgt = Msig(iline) / Luse(1,iline)
      Wtot = Wtot + Wgt
      Mwgt = Mwgt + WgtLin(:,iline) * Wgt
    END DO
    IF (Nlinuse < 1) CALL ErrStop('LineMatrix: no line contributes to the observed spectrum')

    ! report the number of lines used and the weighted mean line weights
    WRITE(*, '(A,I0,A,T75,A)') ' ', Nlinuse, ' lines are used in line-pattern matrix', &
                               '*** LineMatrix'
    DO ilsd = 1, Nlsd
      WRITE(*, '(A,I0,A,ES11.4,T75,A)') ' Mean weight for LSD profile #', ilsd, ': ', &
                                        Mwgt(ilsd) / Wtot, '*** LineMatrix'
    END DO

  END SUBROUTINE LineMatrix

  ! Reconstruct the LSD profiles by solving the regularized least-squares problem
  SUBROUTINE InvertLsd(Mm, Plsd, SigLsd, Msp, Chi, Mde, Covfile)

    TYPE(SparseType), INTENT(IN) :: Mm
    CHARACTER(*), INTENT(IN) :: Covfile
    REAL(8), INTENT(OUT) :: Plsd(Nlsd*Nvlsd), SigLsd(Nlsd*Nvlsd), Msp(Nwl)
    REAL(8), INTENT(OUT) :: Chi, Mde
    INTEGER :: i, j, ilsd, i1, i2, N, Info, Ios
    CHARACTER(256) :: Msg
    REAL(8) :: Reg1, Da1, Da2, MaxDev
    REAL(8) :: Tr1(Nlsd), Tr2(Nlsd), Psm(Nvlsd)
    REAL(8), ALLOCATABLE :: Am(:,:), Ai(:,:), Am0(:,:), Cov(:,:), Prod(:,:)
    REAL(8), ALLOCATABLE :: Rhs(:), Sig2(:), Wcol(:)

    ! total number of elements in the LSD profiles
    N = Nlsd * Nvlsd
    IF (Nwl <= N) CALL ErrStop('InvertLsd: fewer spectral points than LSD profile points')

    ALLOCATE(Am(N,N), Ai(N,N), Rhs(N), Sig2(Nwl), Wcol(Nwl), STAT=Ios)
    IF (Ios /= 0) CALL ErrStop('InvertLsd: cannot allocate the autocorrelation matrix')

    ! inverse variance of every spectral pixel
    Sig2 = SigSp * SigSp

    ! square symmetric autocorrelation matrix Am = M^T S^2 M, both triangles filled
    CALL SparseAutoCorr(Mm, Sig2, Am)

    ! keep the unregularised normal matrix, needed for the covariance
    IF (Iout >= Iocov .AND. Reg > 0d0) THEN
      ALLOCATE(Am0(N,N), STAT=Ios)
      IF (Ios /= 0) CALL ErrStop('InvertLsd: cannot allocate the normal matrix copy')
      Am0 = Am
    END IF

    ! normalized traces of the main diagonal and sub-diagonal of every profile block
    DO ilsd = 1, Nlsd
      i1 = (ilsd - 1) * Nvlsd + 1
      i2 = i1 + Nvlsd - 1
      Tr1(ilsd) = SUM((/ (Am(j,j), j = i1, i2) /)) / Nvlsd
      Tr2(ilsd) = SUM((/ (Am(j,j+1), j = i1, i2 - 1) /)) / (Nvlsd - 1)
    END DO

    ! apply optional first-order Tikhonov regularization by adding the tri-diagonal matrix
    ! diag(1, 2 ... 2, 1) with -1 on both sub-diagonals, scaled by the mean diagonal element
    IF (Reg > 0d0) THEN
      DO ilsd = 1, Nlsd
        Reg1 = Reg * Tr1(ilsd)
        i1 = (ilsd - 1) * Nvlsd + 1
        i2 = i1 + Nvlsd - 1
        DO j = i1, i2
          Am(j,j) = Am(j,j) + 2d0 * Reg1
        END DO
        Am(i1,i1) = Am(i1,i1) - Reg1
        Am(i2,i2) = Am(i2,i2) - Reg1
        DO j = i1, i2 - 1
          Am(j,j+1) = Am(j,j+1) - Reg1
          Am(j+1,j) = Am(j+1,j) - Reg1
        END DO
      END DO
    END IF

    ! Invert the autocorrelation matrix.  A failure here means the matrix is not positive
    ! definite, so the LSD profiles are not determined by the data and anything computed
    ! from them would be meaningless.
    CALL LINEQ(Am, Ai, N, Info)
    IF (Info /= 0) THEN
      WRITE(Msg, '(A,I0,A)') 'InvertLsd: autocorrelation matrix is not positive definite ' &
        // '(LINEQ returned info = ', Info, '), the LSD profiles are not determined by ' &
        // 'the data; check the line mask for repeated or linearly dependent weight columns'
      CALL ErrStop(Msg)
    END IF

    ! maximum deviations of Am * Ai from the unity matrix
    IF (Verify) THEN
      ALLOCATE(Prod(N,N), STAT=Ios)
      IF (Ios /= 0) CALL ErrStop('InvertLsd: cannot allocate the verification matrix')
      Prod = MATMUL(Am, Ai)
      Da2 = 0d0
      DO i = 1, N
        Da2 = MAX(Da2, ABS(Prod(i,i) - 1d0))
        Prod(i,i) = 0d0
      END DO
      Da1 = MAXVAL(ABS(Prod))
      DEALLOCATE(Prod)
      WRITE(*, '(A,2ES12.4,T75,A)') ' Maximum matrix inversion errors: ', Da2, Da1, &
                                    '*** InvertLsd'
    END IF

    ! Write the covariance matrix of the LSD profiles.  Without regularisation this is the
    ! inverse of the autocorrelation matrix; with it the covariance is Ai # (M^T S^2 M) # Ai,
    ! which is what is written.  Element (I,J) with I=(ILSD-1)*NVLSD+IV refers to velocity
    ! bin IV of profile ILSD.  The error bars in the .lsd file are SQRT(diagonal) scaled by
    ! SQRT(chi2) from the header of that file; the matrix here is NOT scaled.
    IF (Iout >= Iocov) THEN
      ALLOCATE(Cov(N,N), STAT=Ios)
      IF (Ios /= 0) CALL ErrStop('InvertLsd: cannot allocate the covariance matrix')
      IF (Reg > 0d0) THEN
        Cov = MATMUL(Ai, MATMUL(Am0, Ai))
        DEALLOCATE(Am0)
      ELSE
        Cov = Ai
      END IF
      OPEN(UNIT=Ucov, FILE=Covfile, STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=Ios)
      IF (Ios /= 0) CALL ErrStop('InvertLsd: cannot open output file: ' // TRIM(Covfile))
      WRITE(Ucov, '(2I6)') Nvlsd, Nlsd
      DO i = 1, N
        WRITE(Ucov, '(1P,5E16.7)', IOSTAT=Ios) Cov(i,:)
        IF (Ios /= 0) CALL ErrStop('InvertLsd: cannot write file: ' // TRIM(Covfile))
      END DO
      CLOSE(Ucov)
      DEALLOCATE(Cov)
    END IF

    ! weighted cross-correlation of the line mask with the observed spectrum, M^T S^2 Y
    CALL SparseCrossCorr(Mm, Sig2, Sp, Rhs)

    ! multiply the inverse autocorrelation matrix by the cross-correlation vector
    Plsd = MATMUL(Ai, Rhs)

    ! apply the default regularization directly to the LSD profile: a normalized three-point
    ! smoothing whose width follows from the sub-diagonal of the autocorrelation matrix
    IF (Reg < 0d0) THEN
      DO ilsd = 1, Nlsd
        IF (Tr1(ilsd) + Tr2(ilsd) <= 0d0) CYCLE
        i1 = (ilsd - 1) * Nvlsd + 1
        i2 = i1 + Nvlsd - 1
        Da1 = Tr1(ilsd) / (Tr1(ilsd) + Tr2(ilsd))
        Da2 = Tr2(ilsd) / (Tr1(ilsd) + Tr2(ilsd))
        Psm(2:Nvlsd-1) = 0.5d0 * Da2 * (Plsd(i1:i2-2) + Plsd(i1+2:i2)) &
                       + Da1 * Plsd(i1+1:i2-1)
        Plsd(i1+1:i2-1) = Psm(2:Nvlsd-1)
      END DO
    END IF

    ! LSD error bars from the diagonal of the inverse autocorrelation matrix
    SigLsd = SQRT(MAX((/ (Ai(i,i), i = 1, N) /), 0d0))

    ! model spectrum M # Z, reduced chi^2, mean standard deviation and maximum deviation
    CALL SparseModel(Mm, Plsd, Msp)
    Chi = SUM(((Msp - Sp) * SigSp)**2) / (Nwl - N)
    IF (Pol == 2) THEN
      ! in Pol=2 mode chi^2 is computed in the log(F) scale, MaxDev and Mde in the F scale
      Wcol = EXP(Msp) - EXP(Sp)
      MaxDev = MAXVAL(ABS(Wcol * SigSp / EXP(Sp)))
      Mde = SQRT(SUM(Wcol * Wcol) / Nwl)
    ELSE
      MaxDev = MAXVAL(ABS((Msp - Sp) * SigSp))
      Mde = SQRT(SUM((Msp - Sp)**2) / Nwl)
    END IF
    WRITE(*, '(A,ES10.4,A,ES10.4,A,ES10.4)') ' MSD = ', Mde, ', MaxDev = ', MaxDev, &
                                             ', reduced chi2 = ', Chi

    ! scale the LSD error bars to make the reduced chi^2 of the model equal to one
    IF (Lscl) THEN
      SigLsd = SigLsd * SQRT(Chi)
      WRITE(*, '(A,ES10.4)') ' LSD error bars are scaled by sqrt(chi2) = ', SQRT(Chi)
    END IF

    DEALLOCATE(Am, Ai, Rhs, Sig2, Wcol)

  END SUBROUTINE InvertLsd

  ! Adjust the individual line weights to improve the fit of the LSD model spectrum
  SUBROUTINE AdjustLsd(Plsd, Msp, Luse, Lcor)

    REAL(8), INTENT(IN) :: Plsd(Nlsd*Nvlsd)
    INTEGER, INTENT(IN) :: Luse(3,Nlin)
    REAL(8), INTENT(INOUT) :: Msp(Nwl)
    REAL(8), INTENT(OUT) :: Lcor(Nlin)
    INTEGER :: i, j, iw, ilsd, iline, i0, i1, i2, Nadj
    REAL(8) :: Dv, Fv, Aq, Bq, X0, Xmin, Xmax
    REAL(8), ALLOCATABLE :: Pvec(:)

    ! lines that do not contribute to any pixel keep their weight unchanged
    Lcor = 0d0
    Nadj = 0
    Xmin = 1d5
    Xmax = -1d5
    ALLOCATE(Pvec(Nwl))

    ! loop over the spectral lines
    DO iline = 1, Nlin
      IF (Luse(1,iline) <= 0) CYCLE

      ! find the first non-zero line component to adjust (different for different lines)
      iw = Nlsd
      DO ilsd = 1, Nlsd
        IF (WgtLin(ilsd,iline) /= 0d0) THEN
          iw = ilsd
          EXIT
        END IF
      END DO
      i0 = (iw - 1) * Nvlsd

      ! first and last spectral point this line reaches
      i1 = Luse(2,iline)
      i2 = Luse(3,iline)

      ! Contribution of a unit increment of this weight to every point the line reaches, and
      ! the two sums that make chi^2 the quadratic function Aq*X^2 + 2*Bq*X + const of the
      ! weight increment X.  Points between the first and the last are tested against the
      ! velocity limits rather than assumed to belong to the line, so that a spectrum whose
      ! wavelengths are not monotonic, as when echelle orders overlap, is treated correctly.
      Aq = 0d0
      Bq = 0d0
      DO i = i1, i2
        Dv = (Wl(i) - WlcLin(iline)) / WlcLin(iline) * Vc
        IF (Dv <= Vfirst .OR. Dv >= Vlast) THEN
          Pvec(i) = 0d0
          CYCLE
        END IF
        j = MIN(MAX(INT((Dv - Vfirst) / Vstep) + 1, 1), Nvlsd - 1)
        Fv = (Dv - Vlsd(j)) / Vstep
        Pvec(i) = Plsd(i0+j) * (1d0 - Fv) + Plsd(i0+j+1) * Fv
        Aq = Aq + (Pvec(i) * SigSp(i))**2
        Bq = Bq + (Msp(i) - Sp(i)) * Pvec(i) * SigSp(i)**2
      END DO
      IF (Aq <= 0d0) CYCLE

      ! minimum of chi^2 plus the prior Wmu*Aq*(weight - input weight)^2, which keeps the
      ! correction of blended lines from running away into their mutual degeneracy
      X0 = (-Bq / Aq - Wmu * (WgtLin(iw,iline) - WgtIni(iw,iline))) / (1d0 + Wmu)

      ! in the intensity modes a line weight is a central depth and cannot become negative;
      ! polarization weights are allowed either sign
      IF (Pol == 0 .OR. Pol == 2) X0 = MAX(X0, -WgtLin(iw,iline))

      ! fold the correction into the model spectrum before moving on to the next line, so
      ! that a line blended with this one sees the residual that is actually left
      Msp(i1:i2) = Msp(i1:i2) + X0 * Pvec(i1:i2)

      ! apply the correction to the initial weight
      WgtLin(iw,iline) = WgtLin(iw,iline) + X0
      Lcor(iline) = X0
      Nadj = Nadj + 1
      Xmax = MAX(Xmax, X0)
      Xmin = MIN(Xmin, X0)
    END DO
    DEALLOCATE(Pvec)

    ! report the number of adjusted lines and the range of the corrections
    IF (Nadj > 0) WRITE(*, '(A,I0,A,2F8.3,T75,A)') ' Adjusted ', Nadj, &
      ' lines, correction range ', Xmin, Xmax, '*** AdjustLsd'

  END SUBROUTINE AdjustLsd

  ! Write the LSD profiles and their error bars
  SUBROUTINE WriteProfiles(Outfile, Plsd, SigLsd, Chi)

    CHARACTER(*), INTENT(IN) :: Outfile
    REAL(8), INTENT(IN) :: Plsd(Nlsd*Nvlsd), SigLsd(Nlsd*Nvlsd), Chi
    INTEGER :: i, j, i0, Ios
    REAL(8) :: Msig

    OPEN(UNIT=Uout, FILE=Outfile, STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('WriteProfiles: cannot open output file: ' // TRIM(Outfile))

    ! the number of velocity bins and of LSD profiles open the first record
    WRITE(Uout, '(I6,I6)', ADVANCE='NO', IOSTAT=Ios) Nvlsd, Nlsd
    IF (Ios /= 0) CALL ErrStop('WriteProfiles: cannot write file: ' // TRIM(Outfile))

    DO i = 1, Nlsd
      i0 = (i - 1) * Nvlsd

      ! average S/N ratio of the LSD profile
      Msig = Nvlsd / SUM(SigLsd(i0+1:i0+Nvlsd))
      WRITE(*, '(A,I0,A,I0,T75,A)') ' Average S/N of LSD profile #', i, ': ', NINT(Msig), &
                                    '*** WriteProfiles'
      WRITE(Uout, '(I7,ES11.4)') NINT(Msig), Chi

      DO j = 1, Nvlsd
        SELECT CASE (Pol)

        CASE (0)
          WRITE(Uout, '(F8.2,2ES14.5)') Vlsd(j), 1d0 - Plsd(j+i0), SigLsd(j+i0)

        CASE (1)
          WRITE(Uout, '(F8.2,2ES14.5)') Vlsd(j), Plsd(j+i0), SigLsd(j+i0)

        CASE (2)
          WRITE(Uout, '(F8.2,2ES14.5)') Vlsd(j), EXP(Plsd(j+i0)), &
            SigLsd(j+i0) * EXP(Plsd(j+i0))

        END SELECT
      END DO
    END DO
    CLOSE(Uout)

  END SUBROUTINE WriteProfiles

  ! Write the LSD model of the retained part of the observed spectrum
  SUBROUTINE WriteModel(Outfile, Msp)

    CHARACTER(*), INTENT(IN) :: Outfile
    REAL(8), INTENT(IN) :: Msp(Nwl)
    INTEGER :: iwl, Ios

    OPEN(UNIT=Uout, FILE=Outfile, STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('WriteModel: cannot open output file: ' // TRIM(Outfile))

    DO iwl = 1, Nwl
      SELECT CASE (Pol)

      CASE (0)
        WRITE(Uout, '(F12.6,3E14.5)') Wl(iwl), 1d0 - Msp(iwl), 1d0 - Sp(iwl), &
          1d0 / SigSp(iwl)

      CASE (1)
        WRITE(Uout, '(F12.6,3E14.5)') Wl(iwl), Msp(iwl), Sp(iwl), 1d0 / SigSp(iwl)

      CASE (2)
        WRITE(Uout, '(F12.6,3E14.5)') Wl(iwl), EXP(Msp(iwl)), EXP(Sp(iwl)), &
          EXP(Sp(iwl)) / SigSp(iwl)

      END SELECT
    END DO
    CLOSE(Uout)

  END SUBROUTINE WriteModel

  ! Write a table of the line central wavelengths and an accompanying per-line quantity
  SUBROUTINE WriteLines(Outfile, Val)

    CHARACTER(*), INTENT(IN) :: Outfile
    REAL(8), INTENT(IN) :: Val(Nlin)
    INTEGER :: iline, Ios

    OPEN(UNIT=Uout, FILE=Outfile, STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=Ios)
    IF (Ios /= 0) CALL ErrStop('WriteLines: cannot open output file: ' // TRIM(Outfile))

    DO iline = 1, Nlin
      WRITE(Uout, '(F11.4,1X,ES11.4)') WlcLin(iline), Val(iline)
    END DO
    CLOSE(Uout)

  END SUBROUTINE WriteLines

  ! Report which output files have been written and what they contain
  SUBROUTINE ReportOutput(Prefix)

    CHARACTER(*), INTENT(IN) :: Prefix
    INTEGER :: N

    N = Nlsd * Nvlsd

    WRITE(*, '(/,A)') ' Output written to:'
    WRITE(*, '(A)') '   ' // TRIM(Prefix) // '.lsd = LSD profiles, error bars, average ' &
      // 'S/N and reduced chi2'
    IF (Iout >= Iomod) THEN
      WRITE(*, '(A)') '   ' // TRIM(Prefix) // '.wln = total observational weight of ' &
        // 'every mask line'
      WRITE(*, '(A)') '   ' // TRIM(Prefix) // '.mod = LSD model, observed spectrum ' &
        // 'and its error bar'
    END IF
    IF (Iout >= Iocov) WRITE(*, '(A,I0,A,I0,A)') '   ' // TRIM(Prefix) // '.cov = ' &
      // 'covariance matrix of the LSD profiles (', N, ' x ', N, ')'
    IF (Niter > 0) WRITE(*, '(A)') '   ' // TRIM(Prefix) // '.cor = accumulated line ' &
      // 'weight corrections'

  END SUBROUTINE ReportOutput

END MODULE LsdCalc

! Main code for least-squares deconvolution
PROGRAM Ilsd

  USE LsdData
  USE LsdUtils
  USE Sparse
  USE LsdCalc

  IMPLICIT NONE

  INTEGER :: iter, N
  REAL(8) :: Chi, Mde, Chiprev
  TYPE(SparseType) :: Mm
  INTEGER, ALLOCATABLE :: Luse(:,:)
  REAL(8), ALLOCATABLE :: Plsd(:), SigLsd(:), Msp(:)
  REAL(8), ALLOCATABLE :: Lwgt(:), Lcor(:), Lcor1(:)
  CHARACTER(1024) :: Prefix

  ! read input data and initialize arrays
  CALL ReadInput(Prefix)

  ! remove spectral pixels not covered by any line
  CALL PreInit

  ! allocate the profiles and the line bookkeeping arrays
  N = Nlsd * Nvlsd
  ALLOCATE(Luse(3,Nlin), Lwgt(Nlin), Lcor(Nlin), Lcor1(Nlin))
  ALLOCATE(Plsd(N), SigLsd(N), Msp(Nwl))
  Lcor1 = 0d0
  Chiprev = HUGE(Chiprev)

  ! reconstruct the LSD profiles, optionally iterating the line weight adjustment
  DO iter = 0, Niter

    ! iterative line weight adjustment
    IF (iter > 0) THEN
      WRITE(*, '(/,A,I0,T75,A)') ' Iteration ', iter, '*** Ilsd'
      CALL AdjustLsd(Plsd, Msp, Luse, Lcor)

      ! accumulate the correction in Lcor1
      Lcor1 = Lcor1 + Lcor
    END IF

    ! compute the line-pattern matrix and save the total weight of every line
    CALL LineMatrix(Mm, Luse, Lwgt)
    IF (Iout >= Iomod) CALL WriteLines(TRIM(Prefix) // '.wln', Lwgt)

    ! compute the LSD profiles; from the Iocov level their covariance matrix is written too
    CALL InvertLsd(Mm, Plsd, SigLsd, Msp, Chi, Mde, TRIM(Prefix) // '.cov')

    ! warn about iteration that makes the fit clearly worse
    IF (iter > 0 .AND. Chi > Chiprev * Chitol) WRITE(*, '(A,ES10.4,A,ES10.4,T75,A)') &
      ' *** WARNING: reduced chi2 rose from ', Chiprev, ' to ', Chi, '*** Ilsd'
    Chiprev = Chi

    ! write the output LSD profiles and the model spectrum
    CALL WriteProfiles(TRIM(Prefix) // '.lsd', Plsd, SigLsd, Chi)
    IF (Iout >= Iomod) CALL WriteModel(TRIM(Prefix) // '.mod', Msp)
  END DO

  ! save the resulting line weight corrections
  IF (Niter > 0) CALL WriteLines(TRIM(Prefix) // '.cor', Lcor1)

  ! report the output files
  CALL ReportOutput(Prefix)

  CALL SparseFree(Mm)
  DEALLOCATE(Luse, Lwgt, Lcor, Lcor1, Plsd, SigLsd, Msp)

END PROGRAM Ilsd
