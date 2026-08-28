! test_hob_su3.f90
!
! Validation driver for the SU(3)-basis HOB code (hob_su3.f90),
! arXiv:2501.19131.
!
! For every E = 0..Emax, L = 0..E and a set of mass-ratio parameters d:
!   * compares every bracket with the classical Talmi-Moshinsky formula
!     evaluated independently by the compact expression of Kamuntavicius et
!     al. (2001), implemented in tmb_kam.f90
!   * checks HOB x HOB^T = I (orthonormality, as in the paper, sec. 6)
!   * checks HOB(d)^2 = I: T(d) is a reflection, so its lift is an
!     involution.  This is independent of orthonormality.
!   * checks the S(3) relations carried by the brackets at d = 1/3, where
!     P23 = HOB(1/3) and P12 = diag((-1)^l1):
!         P23^2 = I        and       (P12 P23)^3 = I
!     Neither test uses any external bracket implementation, and unlike
!     orthonormality they are NOT blind to a spurious orthogonal rotation.
!   * the alpha0 multiplicity check against the Racah formula, eq. (38),
!     is done inside hob_prepare
!
! Usage:  test_hob_su3 [Emax]      (default Emax = 8)
!
program test_hob_su3
    use hob_su3
    use tmb_kam, only: TMB, tmb_kam_init
    implicit none

    integer :: Emax
    character(len=32) :: arg
    real(kind=8), parameter :: dvals(3) = [1.d0, 2.d0, 0.5d0]
    type(hob_solver_t) :: S
    real(kind=8), allocatable :: H(:, :), G(:, :), P(:, :), Q(:, :)
    real(kind=8) :: d, ref, orthmax, invmax, diff, s3max, braidmax, t0, t1, tsum
    integer :: idv, E, L, i, j, n, nbad

    Emax = 8
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) Emax
    end if

    call hob_init(Emax)
    call tmb_kam_init()

    write (*, '(a)') 'SU(3)-basis HOB code (arXiv:2501.19131) - validation'
    write (*, '(a,i3)') 'Emax =', Emax
    write (*, '(a)') ''
    write (*, '(a)') '   d      max|HOB.HOB^T - I|   max|HOB^2 - I|'// &
        '   max|SU3HOB - Kam01|'
    write (*, '(a)') '          (orthonormality)     (involution)'// &
        '     (classical formula)'

    nbad = 0
    do idv = 1, size(dvals)
        d = dvals(idv)
        orthmax = 0.d0
        invmax = 0.d0
        diff = 0.d0
        do E = 0, Emax
            do L = 0, E
                call hob_prepare(E, L, S)
                if (S%b%ntot == 0) cycle
                if (.not. S%ok) then
                    nbad = nbad + 1
                    write (*, '(a,2i4)') 'WARNING: alpha0 mismatch at E, L =', E, L
                end if
                call hob_eval(S, d, H)
                n = S%b%ntot
                G = matmul(H, transpose(H))
                do i = 1, n
                    G(i, i) = G(i, i) - 1.d0
                end do
                orthmax = max(orthmax, maxval(abs(G)))
                G = matmul(H, H)
                do i = 1, n
                    G(i, i) = G(i, i) - 1.d0
                end do
                invmax = max(invmax, maxval(abs(G)))
                do i = 1, n
                    do j = 1, n
                        ref = TMB(E, L, d, &
                                      S%b%e1(i), S%b%l1(i), E - S%b%e1(i), S%b%l2(i), &
                                      S%b%e1(j), S%b%l1(j), E - S%b%e1(j), S%b%l2(j))
                        diff = max(diff, abs(H(i, j) - ref))
                    end do
                end do
            end do
        end do
        write (*, '(f6.2,3(8x,es14.6))') d, orthmax, invmax, diff
    end do
    if (nbad == 0) then
        write (*, '(a)') ''
        write (*, '(a)') 'alpha0 multiplicities agree with the Racah formula,'// &
            ' eq. (38), for all (E, J, L).'
    end if

    ! S(3) relations carried by the d = 1/3 brackets:  P23 = HOB(1/3),
    ! P12 = diag((-1)^l1).  Both identities are exact for exact brackets and,
    ! unlike orthonormality, are sensitive to any error at all.
    write (*, '(a)') ''
    write (*, '(a)') 'S(3) relations at d = 1/3   (P23 = HOB(1/3), P12 = diag((-1)^l1)):'
    write (*, '(a)') '   E    max|P23^2 - I|   max|(P12.P23)^3 - I|'
    do E = 0, Emax
        s3max = 0.d0
        braidmax = 0.d0
        do L = 0, E
            call hob_prepare(E, L, S)
            n = S%b%ntot
            if (n == 0) cycle
            call hob_eval(S, 1.d0/3.d0, H)
            if (allocated(G)) deallocate (G)
            allocate (G(n, n))
            G = matmul(H, H)
            do i = 1, n
                G(i, i) = G(i, i) - 1.d0
            end do
            s3max = max(s3max, maxval(abs(G)))
            ! P = P12 . P23
            if (allocated(P)) deallocate (P)
            if (allocated(Q)) deallocate (Q)
            allocate (P(n, n), Q(n, n))
            do j = 1, n
                do i = 1, n
                    P(i, j) = H(i, j)
                    if (mod(S%b%l1(i), 2) /= 0) P(i, j) = -P(i, j)
                end do
            end do
            Q = matmul(P, matmul(P, P))
            do i = 1, n
                Q(i, i) = Q(i, i) - 1.d0
            end do
            braidmax = max(braidmax, maxval(abs(Q)))
        end do
        write (*, '(i5,2(6x,es14.6))') E, s3max, braidmax
    end do

    ! timing of the SU(3) scheme (cf. Table 1 of the paper): full HOB
    ! matrix for all L at given E
    write (*, '(a)') ''
    write (*, '(a)') 'Timing (all L blocks per E, d = 1):'
    write (*, '(a)') '   E    prepare+eval [s]'
    do E = 0, Emax
        call cpu_time(t0)
        tsum = 0.d0
        do L = 0, E
            call hob_prepare(E, L, S)
            if (S%b%ntot == 0) cycle
            call hob_eval(S, 1.d0, H)
            tsum = tsum + H(1, 1)   ! prevent optimizing the loop away
        end do
        call cpu_time(t1)
        write (*, '(i4,4x,es12.4)') E, t1 - t0
    end do

end program test_hob_su3
