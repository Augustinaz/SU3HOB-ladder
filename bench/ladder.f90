program ladder
    use hob_su3
    implicit none
    type(hob_solver_t) :: S
    real(kind=8), allocatable :: H(:,:), G(:,:)
    real(kind=8) :: t0,t1,tp,te,o
    integer :: E,L,Emax,n,nb,nmax,k
    character(len=32)::arg
    call get_command_argument(1,arg); read(arg,*) Emax
    call hob_init(Emax)
    write(*,'(a)') '#   E  nblk  nmax |  prepare    eval |   max|HH^T-I|'
    do E=0,Emax
      tp=0; te=0; o=0; nb=0; nmax=0
      do L=0,E
        call wall(t0); call hob_prepare(E,L,S); call wall(t1)
        n=S%b%ntot
        if (n==0) cycle
        nb=nb+1; nmax=max(nmax,n); tp=tp+(t1-t0)
        call wall(t0)
        if (allocated(H)) deallocate(H)
        call hob_eval(S,1.d0/3.d0,H); call wall(t1); te=te+(t1-t0)
        if (allocated(G)) deallocate(G); allocate(G(n,n))
        G=matmul(H,transpose(H))
        do k=1,n
          G(k,k)=G(k,k)-1.d0
        end do
        o=max(o,maxval(abs(G)))
      end do
      if (nb==0) cycle
      write(*,'(i5,i6,i6,a,2f9.3,a,es14.3)') E,nb,nmax,' |',tp,te,' |',o
    end do
contains
    subroutine wall(t)
      real(kind=8),intent(out)::t
      integer(kind=8)::c,cr
      call system_clock(c,cr); t=dble(c)/dble(cr)
    end subroutine
end program
