      subroutine t_dsyev
      implicit none
      character jobz, uplo
      integer, parameter :: dp = kind(1.0d0)
      integer :: n, info, lda, lwork
      real(dp), allocatable ::  a(:,:), w(:), work(:)
      real(dp) :: lw
      jobz = 'V'
      uplo = 'L'
      n = 2
      lda = n
      allocate(a(n,n), w(n))
      lwork = -1 
      call dsyev(jobz, uplo, n, a, lda, w, lw, lwork, info)
      lwork = int(lw)
      allocate(work(lwork))
      print *, 'lw = ', lw
      print *, 'info = ', info
      a = 1
      a(1,1) = 2
      a(2,2) = 3
      call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      print *, 'info = ', info
      print *, 'w = ', w
      end
