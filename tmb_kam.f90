! tmb_kam.f90 -- the general harmonic-oscillator transformation bracket.
!
! A fresh, self-contained implementation of the compact expression of
!
!   G. P. Kamuntavicius, R. K. Kalinauskas, B. R. Barrett, S. Mickevicius,
!   D. Germanas, "The general harmonic oscillator brackets: compact
!   expression, symmetries, sums and Fortran code",
!   Nucl. Phys. A695 (2001) 191.
!
! The bracket for two oscillator quanta coupled to L, with mass ratio d, is
! that paper's Eq. (26),
!
!   <e1 l1, e2 l2 : L | e1' l1', e2' l2' : L>_d
!     = (-1)^((l1'+l2'-l1-l2)/2) sqrt( d^(e1-e2') / (1+d)^(e1+e2) )
!       * sum_{ed} (-d)^ed sum_{la lb lc ld}
!            { la lb l1 }
!            { lc ld l2 }   G(e1,l1;ea,la,eb,lb) G(e2,l2;ec,lc,ed,ld)
!            { l1' l2' L }  G(e1',l1';ea,la,ec,lc) G(e2',l2';eb,lb,ed,ld)
!
! with ea = e1-e2'+ed, eb = e2'-ed, ec = e2-ed, and the elementary factor
!
!   G(e,l;ea,la,eb,lb) = <la 0, lb 0 | l 0>
!                        sqrt( (2la+1)(2lb+1) T(e-l, ea-la, eb-lb)
!                                             T(e+l+1, ea+la+1, eb+lb+1) )
!
!   T(i,j,k) = i!! / ( j!! k!! ).
!
! The angular-momentum coefficients are taken from
!
!   D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii, "Quantum Theory
!   of Angular Momentum", World Scientific (1988),
!
! namely the 6j of Sec. 9.2.1 (Racah's single sum), the 9j of Sec. 10.2.4
! (a single sum over products of three 6j), the general Clebsch-Gordan
! coefficient of Sec. 8.2.1 (Racah's single sum), and
! its zero-projection closed form of Sec. 8.5.1.
!
! Everything with a factorial or a double factorial in it is evaluated
! through a table of logarithms.  This is the one substantive departure from
! the published Fortran, which tabulates the binomial and trinomial
! coefficients themselves: those are exact only while they stay below 2**53,
! which caps the classical route at E ~ 12, whereas the log form carries the
! same expression to E = 50 and beyond.  Only the alternating sums of the 6j
! and the ed-sum of the bracket are done in linear arithmetic, and those are
! summed with a running compensation.
!
! Arguments: e, l are ordinary (undoubled) oscillator and orbital quantum
! numbers.  c6j and c9j take DOUBLED angular momenta, as in the published
! code, so that they remain usable for half-integer spins.
!
module tmb_kam
    implicit none
    private
    public :: TMB, tmb_kam_init, c6j, c9j, cg0

    integer, parameter :: NF = 4000
    real(kind=8), save :: lnf(0:NF)          ! ln(n!)
    real(kind=8), save :: ldf(-1:NF)         ! ln(n!!),  (-1)!! = 0!! = 1
    logical, save :: ready = .false.

contains

    subroutine tmb_kam_init()
        integer :: i
        lnf(0) = 0.d0
        do i = 1, NF
            lnf(i) = lnf(i - 1) + log(dble(i))
        end do
        ldf(-1) = 0.d0
        ldf(0) = 0.d0
        do i = 1, NF
            ldf(i) = ldf(i - 2) + log(dble(i))
        end do
        ready = .true.
    end subroutine tmb_kam_init

    ! Triangle test for doubled arguments: 1 if |a-b| <= c <= a+b with
    ! a+b+c even, 0 otherwise.
    integer pure function tri(a, b, c)
        integer, intent(in) :: a, b, c
        integer :: s
        tri = 0
        s = a + b + c
        if (a < 0 .or. b < 0 .or. c < 0) return
        if (mod(s, 2) /= 0) return
        if (c < abs(a - b) .or. c > a + b) return
        tri = 1
    end function tri

    ! ln Delta(a b c) of Varshalovich (8.2.1), doubled arguments.
    real(kind=8) pure function lndelta(a, b, c)
        integer, intent(in) :: a, b, c
        lndelta = 0.5d0*(lnf((a + b - c)/2) + lnf((a - b + c)/2) &
                         + lnf((-a + b + c)/2) - lnf((a + b + c)/2 + 1))
    end function lndelta

    ! 6j symbol, Varshalovich Sec. 9.2.1 (Racah).  Doubled arguments:
    !
    !   { a b c }
    !   { d e f }
    !
    real(kind=8) function c6j(a, b, c, d, e, f) result(ans)
        integer, intent(in) :: a, b, c, d, e, f
        integer :: n1, n2, n3, n4, m1, m2, m3, z, zlo, zhi
        real(kind=8) :: pref, term, s, comp, t
        ans = 0.d0
        if (tri(a, b, c)*tri(a, e, f)*tri(d, b, f)*tri(d, e, c) == 0) return
        n1 = (a + b + c)/2
        n2 = (a + e + f)/2
        n3 = (d + b + f)/2
        n4 = (d + e + c)/2
        m1 = (a + b + d + e)/2
        m2 = (b + c + e + f)/2
        m3 = (c + a + f + d)/2
        zlo = max(n1, n2, n3, n4)
        zhi = min(m1, m2, m3)
        if (zlo > zhi) return
        pref = lndelta(a, b, c) + lndelta(a, e, f) + lndelta(d, b, f) + lndelta(d, e, c)
        s = 0.d0
        comp = 0.d0
        do z = zlo, zhi
            term = exp(pref + lnf(z + 1) - lnf(z - n1) - lnf(z - n2) - lnf(z - n3) &
                       - lnf(z - n4) - lnf(m1 - z) - lnf(m2 - z) - lnf(m3 - z))
            if (mod(z, 2) /= 0) term = -term
            t = s + (term - comp)              ! compensated summation
            comp = (t - s) - (term - comp)
            s = t
        end do
        ans = s
    end function c6j

    ! 9j symbol, Varshalovich Sec. 10.2.4, as a single sum over three 6j.
    ! Doubled arguments, rows (a b c), (d e f), (g h i):
    !
    !   { a b c }              { a d g } { b e h } { c f i }
    !   { d e f } = sum_x (2x+1) { h i x } { d x f } { x a b }
    !   { g h i }
    !
    real(kind=8) function c9j(a, b, c, d, e, f, g, h, i) result(ans)
        integer, intent(in) :: a, b, c, d, e, f, g, h, i
        integer :: x, xlo, xhi
        real(kind=8) :: term, s, comp, t
        ans = 0.d0
        if (tri(a, b, c)*tri(d, e, f)*tri(g, h, i) == 0) return
        if (tri(a, d, g)*tri(b, e, h)*tri(c, f, i) == 0) return
        xlo = max(abs(a - i), abs(h - d), abs(b - f))
        xhi = min(a + i, h + d, b + f)
        s = 0.d0
        comp = 0.d0
        do x = xlo, xhi, 2
            term = dble(x + 1)*c6j(a, d, g, h, i, x)*c6j(b, e, h, d, x, f) &
                   *c6j(c, f, i, x, a, b)
            if (mod(x, 2) /= 0) term = -term         ! (-1)^(2x)
            t = s + (term - comp)
            comp = (t - s) - (term - comp)
            s = t
        end do
        ans = s
    end function c9j

    ! <l1 0, l2 0 | L 0>, Varshalovich Sec. 8.5.1.  UNDOUBLED arguments.
    real(kind=8) function cg0(l1, l2, l3) result(ans)
        integer, intent(in) :: l1, l2, l3
        integer :: j, g
        ans = 0.d0
        if (tri(2*l1, 2*l2, 2*l3) == 0) return
        j = l1 + l2 + l3
        if (mod(j, 2) /= 0) return                   ! odd sum: vanishes
        g = j/2
        ans = exp(lnf(g) - lnf(g - l1) - lnf(g - l2) - lnf(g - l3) &
                  + 0.5d0*(lnf(j - 2*l1) + lnf(j - 2*l2) + lnf(j - 2*l3) - lnf(j + 1))) &
              *sqrt(dble(2*l3 + 1))
        if (mod(g - l3, 2) /= 0) ans = -ans
    end function cg0

    ! The elementary factor G of the compact expression.
    real(kind=8) function gfac(e, l, ea, la, eb, lb) result(ans)
        integer, intent(in) :: e, l, ea, la, eb, lb
        real(kind=8) :: lnt
        ans = 0.d0
        if (tri(2*la, 2*lb, 2*l) == 0) return
        ! ln[ T(e-l, ea-la, eb-lb) T(e+l+1, ea+la+1, eb+lb+1) ]
        lnt = ldf(e - l) - ldf(ea - la) - ldf(eb - lb) &
              + ldf(e + l + 1) - ldf(ea + la + 1) - ldf(eb + lb + 1)
        ans = cg0(la, lb, l)*sqrt(dble((2*la + 1)*(2*lb + 1))*exp(lnt))
    end function gfac

    ! The bracket itself:  <e1 l1, e2 l2 : L | e1s l1s, e2s l2s : L>_d
    !
    ! Argument order and phase convention match CalcTMB of the HOTB tool, so
    ! this is a drop-in replacement for it.  E is accepted and ignored; the
    ! shell is fixed by e1+e2 = e1s+e2s, which is checked here.
    real(kind=8) function TMB(E, L, d, e1, l1, e2, l2, e1s, l1s, e2s, l2s) result(ans)
        integer, intent(in) :: E, L, e1, l1, e2, l2, e1s, l1s, e2s, l2s
        real(kind=8), intent(in) :: d
        integer :: ed, ea, eb, ec, la, lb, lc, ld, edmax, ph
        real(kind=8) :: lnscale, sgn, term, s, comp, t, w
        ans = 0.d0
        if (.not. ready) call tmb_kam_init()
        if (E < 0) continue
        if (e1 + e2 /= e1s + e2s) return
        if (tri(2*l1, 2*l2, 2*L)*tri(2*l1s, 2*l2s, 2*L) == 0) return

        ! sqrt( d^(e1-e2s) / (1+d)^(e1+e2) ), in logarithms
        lnscale = 0.5d0*(dble(e1 - e2s)*log(d) - dble(e1 + e2)*log(1.d0 + d))

        ph = (l1s + l2s - l1 - l2)/2
        sgn = 1.d0
        if (ph/2*2 /= ph) sgn = -1.d0

        edmax = min(e2s, e2)
        s = 0.d0
        comp = 0.d0
        do ed = 0, edmax
            eb = e2s - ed
            ec = e2 - ed
            ea = e1 - e2s + ed
            if (ea < 0) cycle
            ! (-d)^ed folded into the scale, in logarithms
            w = exp(lnscale + dble(ed)*log(d))
            if (mod(ed, 2) /= 0) w = -w
            do ld = ed, 0, -2
                do lb = eb, 0, -2
                    if (tri(2*ld, 2*lb, 2*l2s) == 0) cycle
                    do lc = ec, 0, -2
                        if (tri(2*ld, 2*lc, 2*l2) == 0) cycle
                        do la = ea, 0, -2
                            if (tri(2*la, 2*lb, 2*l1) == 0) cycle
                            if (tri(2*la, 2*l1s, 2*lc) == 0) cycle
                            term = w*c9j(2*la, 2*lb, 2*l1, 2*lc, 2*ld, 2*l2, &
                                         2*l1s, 2*l2s, 2*L) &
                                   *gfac(e1, l1, ea, la, eb, lb) &
                                   *gfac(e2, l2, ec, lc, ed, ld) &
                                   *gfac(e1s, l1s, ea, la, ec, lc) &
                                   *gfac(e2s, l2s, eb, lb, ed, ld)
                            t = s + (term - comp)
                            comp = (t - s) - (term - comp)
                            s = t
                        end do
                    end do
                end do
            end do
        end do
        ans = sgn*s
    end function TMB

end module tmb_kam
