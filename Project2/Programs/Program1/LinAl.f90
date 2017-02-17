module LinAl
  implicit none

contains

  real function norm2(vec,vsize)
    !function to compute the (Euclidean)norm of a vector vec of size vsize  !outputs 
    integer :: vsize, i
    real, dimension(vsize) :: vec
    norm2 = 0.
    do i=1,vsize
       norm2 = norm2 + vec(i)**2
       !write(*,*) norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

  real function normF(A,msize,nsize)
    !function to compute the Frobenius norm of a matrix of size dimensions (msize,nsize) 
    integer :: msize, nsize,i,j
    real, dimension(msize,nsize) :: A
    normF = 0.
    do i=1,msize
       do j=1,nsize
          normF = normF + A(i,j)**2
          !write(*,*) norm2
       enddo
    enddo
    normF = sqrt(normF)
  end function normF

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in row-wise split format
    integer :: msize, nsize, i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat

  subroutine Choldec(A,nsize,ipos,ising)

  integer :: nsize
  real, dimension(nsize,nsize) :: A
  integer :: i,j,k
  real :: sum, eps=1e-8
  logical :: ipos, ising

! ***** Performs Cholesky decomposition of a given matrix **! 
! **** Inputs are:
! A = square matrix of size n 
! ising = logical flag for checking if A is singular
! ipos = logical flag for checking if A is positive definite
!***********************************************************
  do j=1,nsize
     sum = A(j,j)
     do k=1,j-1
        sum = sum - A(j,k)**2
     enddo
     if(sum.le.epsilon(0.)) then
        !write(*,*) k,sum
        write (*,*) 'matrix is not positive definite'
        ipos=.TRUE.
        goto 2000
      endif
     A(j,j) = sqrt(sum)
     if(A(j,j).le.epsilon(0.)) then
        write (*,*) 'singular matrix'
        ising =.TRUE.
        goto 2000
      endif
     do i=j+1,nsize
        sum=a(i,j)
        do k=1,j-1
           sum=sum-A(i,k)*A(j,k)
        enddo
        A(i,j) = sum/A(j,j)
     enddo
  enddo
2000 return

end subroutine Choldec

! ************************************************

subroutine Cholsol(A,nsize,b)

  integer :: nsize
  real, dimension(nsize,nsize) :: A
  real, dimension(nsize) :: b

  integer :: i,k
  real :: sum
      
! ************************************************
! Subroutine to solve Ax=b for RHS vector b
! Solution is returned in b
! Inputs:
! Cholesky factorized matrix A (with L 
! in lower triangle) and its size n
! b = RHS of Ax = b
! **************************************************


! *** Forward substitution loop with L 
  do i=1,nsize
     sum=b(i)
     do k=1,i-1
        sum=sum-A(i,k)*b(k)
     enddo
     b(i)=sum/A(i,i)
  enddo

! *** Back substitution loop with L^T
  do i=nsize,1,-1
     sum=b(i)
     do k=i+1,nsize
        sum=sum-A(k,i)*b(k)
     enddo
     b(i)=sum/A(i,i)
  enddo

  return
end subroutine Cholsol

end module LinAl
