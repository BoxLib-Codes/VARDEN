      subroutine t_dsygv
      implicit none
      character jobz, uplo
      integer, parameter :: dp = kind(1.0d0)
      integer :: n, info, lwork
      real(dp), allocatable ::  a(:,:), w(:), work(:), b(:,:)
      real(dp) :: lw
      integer :: itype, i
      itype = 1
      jobz = 'V'
      uplo = 'L'
      n = 2
      allocate(a(n,n), w(n),b(n,n))
      b = 0
      forall(i=1:n) b(i,i) = 1
      lwork = -1 
      call dsygv(itype, jobz, uplo, n, a, n, b, n, w, 
     &  lw, lwork, info)
      lwork = int(lw)
      allocate(work(lwork))
      print *, 'lw = ', lw
      print *, 'info = ', info
      a = 1
      a(1,1) = 2
      a(2,2) = 3
      call dsygv(itype, jobz, uplo, n, a, n, b, n, w, 
     &  work, lwork, info)
      print *, 'info = ', info
      print *, 'w = ', w
      end
