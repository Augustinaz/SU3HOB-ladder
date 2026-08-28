! hob.f90
!
! Command line interface to the SU(3)-basis HOB code (hob_su3.f90),
! arXiv:2501.19131.
!
! Usage:
!   hob E L d                              full HOB matrix of the (E, L) block
!   hob e1 l1 e2 l2 e1' l1' e2' l2' L d    single bracket
!                                          <e1 l1 e2 l2 : L|T(d)|e1'l1'e2'l2' : L>
!
program hob
    use hob_su3
    implicit none

    integer :: nargs, E, L, q(9), i, j
    real(kind=8) :: d, val
    character(len=32) :: arg
    type(hob_solver_t) :: S
    real(kind=8), allocatable :: H(:, :)

    nargs = command_argument_count()
    if (nargs == 3) then
        call get_command_argument(1, arg)
        read (arg, *) E
        call get_command_argument(2, arg)
        read (arg, *) L
        call get_command_argument(3, arg)
        read (arg, *) d
        call hob_init(E)
        call hob_prepare(E, L, S)
        if (S%b%ntot == 0) then
            write (*, '(a)') 'empty (E, L) block'
            stop
        end if
        if (.not. S%ok) write (*, '(a)') 'WARNING: multiplicity check failed'
        call hob_eval(S, d, H)
        write (*, '(a,i3,a,i3,a,f10.6)') 'E =', E, '  L =', L, '  d =', d
        write (*, '(a)') '  e1 l1 e2 l2   e1''l1''e2''l2''     bracket'
        do i = 1, S%b%ntot
            do j = i, S%b%ntot
                write (*, '(2x,4i3,3x,4i3,f16.10)') &
                    S%b%e1(i), S%b%l1(i), E - S%b%e1(i), S%b%l2(i), &
                    S%b%e1(j), S%b%l1(j), E - S%b%e1(j), S%b%l2(j), H(i, j)
            end do
        end do
    else if (nargs == 10) then
        do i = 1, 9
            call get_command_argument(i, arg)
            read (arg, *) q(i)
        end do
        call get_command_argument(10, arg)
        read (arg, *) d
        E = q(1) + q(3)
        L = q(9)
        if (q(5) + q(7) /= E) then
            write (*, '(a)') 'error: e1+e2 /= e1''+e2'''
            stop 1
        end if
        call hob_init(E)
        call hob_prepare(E, L, S)
        val = hob_bracket(S, q(1), q(2), q(3), q(4), q(5), q(6), q(7), q(8), d)
        write (*, '(f18.12)') val
    else
        write (*, '(a)') 'usage: hob E L d'
        write (*, '(a)') '       hob e1 l1 e2 l2 e1'' l1'' e2'' l2'' L d'
        stop 1
    end if

end program hob
