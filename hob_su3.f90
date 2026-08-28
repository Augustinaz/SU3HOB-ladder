! hob_su3.f90
!
! Harmonic-oscillator (Talmi-Moshinsky) brackets computed in the SU(3)
! basis, following:
!
!   R.K. Kalinauskas, A. Stepsys, D. Germanas, S. Mickevicius,
!   "Calculation of harmonic oscillator brackets in SU(3) basis",
!   arXiv:2501.19131 [nucl-th].
!
! Central result implemented here is eq. (49) of the paper:
!
!   <e1 l1 e2 l2 L| T(d) |e1' l1' e2' l2' L> =
!       sum_{J}  Delta0^E2 * d^J_{M,M'}(d) *
!       sum_{alpha=1..alpha0} <l1 l2|J alpha>^{e1 e2 L}
!                             <l1' l2'|J alpha>^{e1' e2' L}
!
! with  E = e1+e2 = e1'+e2',  M = (e1-e2)/2,  M' = (e1'-e2')/2,
! E2 = E/2 - J and Delta0 = det(TM) = -1 (eq. (28)).  d^J_{M,M'} is the
! Wigner D-matrix, eq. (31), of the reflection-containing Talmi-Moshinsky
! matrix, i.e. with a1 = cos(theta), a2 = a3 = sin(theta),
! a4 = -cos(theta)  (eqs. (32)-(34)) and cot^2(theta) = d (eq. (24)).
!
! The SU(3) > SO(3) isofactors <l1 l2|J alpha>^{e1 e2 L} of the group
! chain U(6) > U(3) x U(2) (eq. (12)) are generated numerically:
! the U(2) raising operator J+ = a+(1).a(2), which is an SO(3) scalar,
! is built in the coupled basis |e1 l1, e2 l2 : L M_L>; for every (E,J,L)
! the highest weight states (kernel of J+ in the block M = J) are found
! and lowered with J- = (J+)^T, which fixes a consistent alpha labelling
! across all M.  The multiplicity alpha0 is verified against the Racah
! formula, eq. (38).
!
! The reduced matrix elements of a+ used for J+ are
!     <e+1 l+1||a+||e l> = sqrt( (l+1)(e+l+3) / (2l+3) )
!     <e+1 l-1||a+||e l> = sqrt(  l   (e-l+2) / (2l-1) )
! (Clebsch-Gordan convention  <n'l'm'|a+_mu|nlm> = <l m 1 mu|l' m'> <||.||>).
! Taking both branches positive selects the radial phase convention of the
! reference HOB codes (Buck-Merchant / Kamuntavicius et al. [5]); with the
! (-1)^n convention of eq. (2) the l-1 branch is negative, which changes
! every bracket by (-1)^(n1+n2+n1'+n2').
!
module hob_su3
    use WignerSymbol, only: CG, wigner_init
    implicit none
    private

    public :: hob_init, hob_build_basis, hob_prepare, hob_eval
    public :: hob_bracket, hob_state_index, alpha_mult
    public :: hob_basis_t, hob_solver_t

    type :: mat_t
        real(kind=8), allocatable :: a(:, :)
    end type mat_t

    ! coupled two-HO basis |e1 l1, e2 l2 : L> at fixed E = e1+e2,
    ! grouped in blocks of fixed e1 (i.e. fixed M = (e1-e2)/2)
    type :: hob_basis_t
        integer :: E = -1, L = -1
        integer :: ntot = 0
        integer, allocatable :: nblk(:)              ! (0:E) block sizes
        integer, allocatable :: off(:)               ! (0:E) block offsets
        integer, allocatable :: e1(:), l1(:), l2(:)  ! (ntot) state labels
    end type hob_basis_t

    ! one SU(3) irrep tower:  the isofactor vectors <l1 l2|J alpha M>
    ! for all M (blocks e1lo..e1hi), nblk(e1) x nalpha each
    type :: tower_t
        integer :: twoJ = 0, nalpha = 0
        integer :: e1lo = 0, e1hi = 0
        type(mat_t), allocatable :: V(:)
    end type tower_t

    type :: hob_solver_t
        type(hob_basis_t) :: b
        integer :: ntower = 0
        type(tower_t), allocatable :: tw(:)
        logical :: ok = .true.    ! alpha0 agreed with Racah formula eq. (38)
    end type hob_solver_t

contains

    subroutine hob_init(Emax)
        integer, intent(in) :: Emax
        call wigner_init(max(100, 8*Emax + 20), "nmax", 3)
    end subroutine hob_init

    ! Racah multiplicity formula, eqs. (37)-(38):
    ! number alpha0 of occurrences of L in the SU(3) irrep with
    ! lambda = 2J, mu = (E-2J)/2
    integer function alpha_mult(E, twoJ, L) result(a0)
        integer, intent(in) :: E, twoJ, L
        integer :: lm, mu
        lm = twoJ
        mu = (E - twoJ)/2
        a0 = max(0, (lm + mu + 2 - L)/2) - max(0, (lm + 1 - L)/2) &
             - max(0, (mu + 1 - L)/2)
    end function alpha_mult

    subroutine hob_build_basis(E, L, b)
        integer, intent(in) :: E, L
        type(hob_basis_t), intent(out) :: b
        integer :: ie, jl1, jl2, cnt, idx

        b%E = E
        b%L = L
        allocate (b%nblk(0:E), b%off(0:E))
        b%ntot = 0
        do ie = 0, E
            cnt = 0
            do jl1 = mod(ie, 2), ie, 2
                do jl2 = mod(E - ie, 2), E - ie, 2
                    if (jl2 >= abs(L - jl1) .and. jl2 <= L + jl1) cnt = cnt + 1
                end do
            end do
            b%nblk(ie) = cnt
            b%off(ie) = b%ntot
            b%ntot = b%ntot + cnt
        end do
        allocate (b%e1(b%ntot), b%l1(b%ntot), b%l2(b%ntot))
        idx = 0
        do ie = 0, E
            do jl1 = mod(ie, 2), ie, 2
                do jl2 = mod(E - ie, 2), E - ie, 2
                    if (jl2 >= abs(L - jl1) .and. jl2 <= L + jl1) then
                        idx = idx + 1
                        b%e1(idx) = ie
                        b%l1(idx) = jl1
                        b%l2(idx) = jl2
                    end if
                end do
            end do
        end do
    end subroutine hob_build_basis

    integer function hob_state_index(b, e1, l1, l2) result(idx)
        type(hob_basis_t), intent(in) :: b
        integer, intent(in) :: e1, l1, l2
        integer :: i
        idx = 0
        if (e1 < 0 .or. e1 > b%E) return
        do i = b%off(e1) + 1, b%off(e1) + b%nblk(e1)
            if (b%l1(i) == l1 .and. b%l2(i) == l2) then
                idx = i
                return
            end if
        end do
    end function hob_state_index

    ! <e+1 lp || a+ || e l>, CG convention (see header note on phases)
    real(kind=8) function spa(e, l, lp) result(a)
        integer, intent(in) :: e, l, lp
        a = 0.d0
        if (e < 0 .or. l < 0 .or. l > e) return
        if (lp == l + 1) then
            a = sqrt(dble((l + 1)*(e + l + 3))/dble(2*l + 3))
        else if (lp == l - 1 .and. l >= 1) then
            a = sqrt(dble(l*(e - l + 2))/dble(2*l - 1))
        end if
    end function spa

    ! <e1+1 jl1, e2-1 jl2 : L L| J+ |e1 il1, e2 il2 : L L>,  e2 = E-e1
    real(kind=8) function jp_elem(L, E, ie, il1, il2, jl1, jl2) result(s)
        integer, intent(in) :: L, E, ie, il1, il2, jl1, jl2
        real(kind=8) :: A1, A2
        integer :: m1, mu, m1p, m2, m2p
        s = 0.d0
        A1 = spa(ie, il1, jl1)
        if (A1 == 0.d0) return
        A2 = spa(E - ie - 1, jl2, il2)
        if (A2 == 0.d0) return
        do m1 = -il1, il1
            m2 = L - m1
            if (abs(m2) > il2) cycle
            do mu = -1, 1
                m1p = m1 + mu
                if (abs(m1p) > jl1) cycle
                m2p = L - m1p
                if (abs(m2p) > jl2) cycle
                s = s + CG(2*jl1, 2*jl2, 2*L, 2*m1p, 2*m2p, 2*L) &
                      *CG(2*il1, 2*il2, 2*L, 2*m1, 2*m2, 2*L) &
                      *CG(2*il1, 2, 2*jl1, 2*m1, 2*mu, 2*m1p) &
                      *CG(2*jl2, 2, 2*il2, 2*m2p, 2*mu, 2*m2)
            end do
        end do
        s = s*A1*A2
    end function jp_elem

    ! Wigner D-matrix element of the Talmi-Moshinsky matrix, eq. (31)
    ! with a1 = cos(theta), a2 = a3 = sin(theta), a4 = -cos(theta),
    ! cot^2(theta) = d (eqs. (22)-(24), (32)-(34))
    real(kind=8) function dtm(twoJ, twoM, twoMp, d) result(val)
        integer, intent(in) :: twoJ, twoM, twoMp
        real(kind=8), intent(in) :: d
        real(kind=8) :: ct, st, pref, term
        integer :: jpm, jmm, jpmp, jmmp, dmm, k, kmin, kmax
        val = 0.d0
        jpm = (twoJ + twoM)/2
        jmm = (twoJ - twoM)/2
        jpmp = (twoJ + twoMp)/2
        jmmp = (twoJ - twoMp)/2
        if (jpm < 0 .or. jmm < 0 .or. jpmp < 0 .or. jmmp < 0) return
        ct = sqrt(d/(1.d0 + d))
        st = sqrt(1.d0/(1.d0 + d))
        dmm = (twoM - twoMp)/2
        pref = 0.5d0*(lfact(jpm) + lfact(jmm) + lfact(jpmp) + lfact(jmmp))
        kmin = max(0, -dmm)
        kmax = min(jmm, jpmp)
        do k = kmin, kmax
            term = exp(pref - lfact(k) - lfact(jmm - k) - lfact(jpmp - k) &
                       - lfact(k + dmm))
            term = term*ct**(jpmp - k + jmm - k)*st**(2*k + dmm)
            if (mod(jmm - k, 2) /= 0) term = -term
            val = val + term
        end do
    end function dtm

    real(kind=8) function lfact(n)
        integer, intent(in) :: n
        lfact = log_gamma(dble(n) + 1.d0)
    end function lfact

    ! eigen decomposition of a small symmetric matrix (cyclic Jacobi)
    subroutine jacobi_eig(n, A, V, w)
        integer, intent(in) :: n
        real(kind=8), intent(inout) :: A(n, n)
        real(kind=8), intent(out) :: V(n, n), w(n)
        integer :: i, p, q, sweep
        real(kind=8) :: off, frob, theta, t, c, s
        real(kind=8) :: tp(n), tq(n)

        V = 0.d0
        do i = 1, n
            V(i, i) = 1.d0
        end do
        if (n <= 1) then
            if (n == 1) w(1) = A(1, 1)
            return
        end if
        frob = sqrt(sum(A*A))
        do sweep = 1, 200
            off = 0.d0
            do p = 1, n - 1
                do q = p + 1, n
                    off = off + A(p, q)**2
                end do
            end do
            if (sqrt(2.d0*off) <= 1.d-14*max(frob, 1.d-30)) exit
            do p = 1, n - 1
                do q = p + 1, n
                    if (abs(A(p, q)) <= 1.d-300) cycle
                    theta = (A(q, q) - A(p, p))/(2.d0*A(p, q))
                    t = sign(1.d0, theta)/(abs(theta) + sqrt(theta**2 + 1.d0))
                    c = 1.d0/sqrt(t*t + 1.d0)
                    s = t*c
                    tp = A(:, p)
                    tq = A(:, q)
                    A(:, p) = c*tp - s*tq
                    A(:, q) = s*tp + c*tq
                    tp = A(p, :)
                    tq = A(q, :)
                    A(p, :) = c*tp - s*tq
                    A(q, :) = s*tp + c*tq
                    tp = V(:, p)
                    tq = V(:, q)
                    V(:, p) = c*tp - s*tq
                    V(:, q) = s*tp + c*tq
                end do
            end do
        end do
        do i = 1, n
            w(i) = A(i, i)
        end do
    end subroutine jacobi_eig

    ! orthonormal basis of the null space of an m x n matrix A
    subroutine null_space(m, n, A, K, nk)
        integer, intent(in) :: m, n
        real(kind=8), intent(in) :: A(m, n)
        real(kind=8), allocatable, intent(out) :: K(:, :)
        integer, intent(out) :: nk
        real(kind=8) :: G(n, n), Vv(n, n), w(n), tol
        integer :: i, j

        G = matmul(transpose(A), A)
        call jacobi_eig(n, G, Vv, w)
        ! nonzero eigenvalues of (J+)^T(J+) are J'(J'+1)-M(M+1) >= 2,
        ! so the spectral gap is clean
        tol = 1.d-8*max(1.d0, maxval(w))
        nk = count(w < tol)
        allocate (K(n, nk))
        j = 0
        do i = 1, n
            if (w(i) < tol) then
                j = j + 1
                K(:, j) = Vv(:, i)
            end if
        end do
    end subroutine null_space

    ! build the basis and all isofactor towers for given (E, L);
    ! this is the d-independent part (SO(3) -> SU(3) transformation),
    ! reusable for any number of Talmi-Moshinsky parameters d
    subroutine hob_prepare(E, L, S)
        integer, intent(in) :: E, L
        type(hob_solver_t), intent(out) :: S
        type(mat_t), allocatable :: jp(:)
        real(kind=8), allocatable :: K(:, :)
        real(kind=8) :: xn
        integer :: ie, i, j, gi, gj, twoJ, e1h, e1lo, nh, nk, it, twoMup

        call hob_build_basis(E, L, S%b)
        allocate (S%tw(E/2 + 1))
        S%ntower = 0
        if (S%b%ntot == 0) return

        ! J+ block matrices: jp(ie) maps block ie -> block ie+1
        allocate (jp(0:max(0, E - 1)))
        do ie = 0, E - 1
            allocate (jp(ie)%a(S%b%nblk(ie + 1), S%b%nblk(ie)))
            do j = 1, S%b%nblk(ie)
                gj = S%b%off(ie) + j
                do i = 1, S%b%nblk(ie + 1)
                    gi = S%b%off(ie + 1) + i
                    jp(ie)%a(i, j) = jp_elem(L, E, ie, S%b%l1(gj), S%b%l2(gj), &
                                             S%b%l1(gi), S%b%l2(gi))
                end do
            end do
        end do

        it = 0
        do twoJ = E, mod(E, 2), -2
            e1h = (E + twoJ)/2
            e1lo = (E - twoJ)/2
            nh = S%b%nblk(e1h)
            if (nh == 0) then
                if (alpha_mult(E, twoJ, L) /= 0) S%ok = .false.
                cycle
            end if
            ! highest weight states: kernel of J+ in the block M = J
            if (e1h == E) then
                nk = nh
                if (allocated(K)) deallocate (K)
                allocate (K(nh, nh))
                K = 0.d0
                do i = 1, nh
                    K(i, i) = 1.d0
                end do
            else if (S%b%nblk(e1h + 1) == 0) then
                nk = nh
                if (allocated(K)) deallocate (K)
                allocate (K(nh, nh))
                K = 0.d0
                do i = 1, nh
                    K(i, i) = 1.d0
                end do
            else
                call null_space(S%b%nblk(e1h + 1), nh, jp(e1h)%a, K, nk)
            end if
            if (nk /= alpha_mult(E, twoJ, L)) S%ok = .false.
            if (nk == 0) cycle
            it = it + 1
            S%tw(it)%twoJ = twoJ
            S%tw(it)%nalpha = nk
            S%tw(it)%e1lo = e1lo
            S%tw(it)%e1hi = e1h
            allocate (S%tw(it)%V(e1lo:e1h))
            S%tw(it)%V(e1h)%a = K(:, 1:nk)
            ! lower with J- = (J+)^T down to M = -J
            do ie = e1h - 1, e1lo, -1
                twoMup = 2*(ie + 1) - E
                xn = sqrt(dble(twoJ*(twoJ + 2) - twoMup*(twoMup - 2))/4.d0)
                S%tw(it)%V(ie)%a = matmul(transpose(jp(ie)%a), &
                                          S%tw(it)%V(ie + 1)%a)/xn
            end do
        end do
        S%ntower = it
    end subroutine hob_prepare

    ! full HOB matrix for the prepared (E, L) and given d, eq. (49)
    subroutine hob_eval(S, d, H)
        type(hob_solver_t), intent(in) :: S
        real(kind=8), intent(in) :: d
        real(kind=8), allocatable, intent(out) :: H(:, :)
        integer :: it, ie, je, twoJ, ni, nj, oi, oj
        real(kind=8) :: w, ph

        allocate (H(S%b%ntot, S%b%ntot))
        H = 0.d0
        do it = 1, S%ntower
            twoJ = S%tw(it)%twoJ
            ! Delta0^E2 with Delta0 = det(TM) = -1, E2 = (E-2J)/2
            ph = 1.d0
            if (mod((S%b%E - twoJ)/2, 2) /= 0) ph = -1.d0
            do ie = S%tw(it)%e1lo, S%tw(it)%e1hi
                ni = S%b%nblk(ie)
                if (ni == 0) cycle
                oi = S%b%off(ie)
                do je = S%tw(it)%e1lo, S%tw(it)%e1hi
                    nj = S%b%nblk(je)
                    if (nj == 0) cycle
                    oj = S%b%off(je)
                    w = ph*dtm(twoJ, 2*ie - S%b%E, 2*je - S%b%E, d)
                    H(oi + 1:oi + ni, oj + 1:oj + nj) = &
                        H(oi + 1:oi + ni, oj + 1:oj + nj) &
                        + w*matmul(S%tw(it)%V(ie)%a, transpose(S%tw(it)%V(je)%a))
                end do
            end do
        end do
    end subroutine hob_eval

    ! single bracket <e1 l1, e2 l2 : L| T(d) |e1p l1p, e2p l2p : L>
    real(kind=8) function hob_bracket(S, e1, l1, e2, l2, e1p, l1p, e2p, l2p, d) &
        result(val)
        type(hob_solver_t), intent(in) :: S
        integer, intent(in) :: e1, l1, e2, l2, e1p, l1p, e2p, l2p
        real(kind=8), intent(in) :: d
        integer :: it, ib, ik, twoJ
        real(kind=8) :: ph

        val = 0.d0
        if (e1 + e2 /= S%b%E .or. e1p + e2p /= S%b%E) return
        ib = hob_state_index(S%b, e1, l1, l2)
        ik = hob_state_index(S%b, e1p, l1p, l2p)
        if (ib == 0 .or. ik == 0) return
        ib = ib - S%b%off(e1)
        ik = ik - S%b%off(e1p)
        do it = 1, S%ntower
            if (e1 < S%tw(it)%e1lo .or. e1 > S%tw(it)%e1hi) cycle
            if (e1p < S%tw(it)%e1lo .or. e1p > S%tw(it)%e1hi) cycle
            twoJ = S%tw(it)%twoJ
            ph = 1.d0
            if (mod((S%b%E - twoJ)/2, 2) /= 0) ph = -1.d0
            val = val + ph*dtm(twoJ, 2*e1 - S%b%E, 2*e1p - S%b%E, d) &
                  *dot_product(S%tw(it)%V(e1)%a(ib, :), S%tw(it)%V(e1p)%a(ik, :))
        end do
    end function hob_bracket

end module hob_su3
