module LinAl
  implicit none

contains

  !********************************************************

  subroutine readandallocatemat(mat,msize,nsize,filename)

    character*100 filename
    real, dimension(:,:), allocatable, intent(in out) :: mat
    integer :: msize, nsize

    integer :: i,j

    ! Reads a file containing the matrix mat
    ! Note that the first 2 lines are the matrix dimensions, 
    ! then the next msize lines are the matrix entries
    ! Note that entries must be separated by a tab.
    ! Then allocates an array of size (msize,nsize), populates the matrix,
    ! and returns the array. 

    ! This routine takes as INPUTS:
    ! filename = a character string with the name of file to read
    ! This routine returns as OUTPUTS:
    ! msize = first dimension of read matrix
    ! nsize = second dimension of read matrix
    ! mat = read matrix. Note that mat type is real.

    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) msize,nsize

    ! Allocate matrix
    allocate(mat(msize,nsize))

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)

  end subroutine readandallocatemat


  subroutine invM(mat,msize,r,z)
    ! calculates z = M^-1*r where M is the diagonal matrix containing the diagonal elements of input matrix mat of size msize 
    integer :: msize, i
    real, dimension(msize,msize) :: mat
    real, dimension(msize) :: r,z
    !Intialize z to 0
    z = 0.
    do i=1,msize
       z(i) =  1.0/mat(i,i)*r(i)
    enddo
  end subroutine invM

  real function norm2(vec,vsize)
    !function to compute the (Euclidean)norm of a vector vec of size vsize
    integer :: vsize, i
    real, dimension(vsize) :: vec
    norm2 = 0.
    do i=1,vsize
       norm2 = norm2 + vec(i)**2
       !write(*,*) norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in row-wise split format 
    integer :: msize, nsize, i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat
   
  subroutine ConjGrad(A,b,x0,msize,acc)
  
  real, dimension(msize,msize) :: A
  integer :: msize,i,j,k=1
  real, dimension(msize) :: x0,b,p,r
  real :: acc,alpha,beta,E,E_new
  
! *** Conjugate Gradient algorithm for solving the linear system ****!
! *** Ax = b with A as a matrix of size msize and b is the RHS vector **! 
     
! *** Initialization for p and r ******************!
      
  r = matmul(A,x0)
  p = b - r
  r = p

  alpha = 0.
  beta = 0.
  
! *** Calculate r^T r
  E = dot_product(r,r)
  do while(sqrt(E).gt.acc) 
  !write(*,*) E
! ****** Store Ap_k in b
  b = matmul(A,p)

!****** Calculate alpha_k
  alpha = E/dot_product(b,p)

! ****** Update x and r
     x0 = x0 + alpha*p
     r = r - alpha*b
! ****** Calculate r_k^T r_k
     E_new = dot_product(r,r)
               
! ****** Update beta and err
     beta = E_new/E
     E = E_new
     write(*,*) "Error at iteration",k,"= ",sqrt(E)
     
! ****** Update p_k
     p = r + beta*p
     k = k+1   
  enddo

end subroutine conjgrad

subroutine precondconjgrad(A,b,x0,msize,acc,invM)
  
  external invM
  real, dimension(msize,msize) :: A
  integer :: msize,i,j,k=1
  real, dimension(msize) :: x0,b,p,r,z
  real :: acc,alpha,beta,E,E_new
  
! *** Conjugate Gradient algorithm using diagonal preconditioning ***!! ***  for solving the linear system  Ax = b with A as a matrix *****!
! *** of size msize and b is the RHS vector **!
 
! *** Initialization for p and r ******************!
      
  r = b - matmul(A,x0)
  call invM(A,msize,r,z)
  p = z
  
  alpha = 0.
  beta = 0.
  
  
! *** Calculate r^T r
  E = dot_product(r,z)
  do while(sqrt(E).gt.acc) 
  !write(*,*) E

! ****** Store Ap_k in b
  b = matmul(A,p)

!****** Calculate alpha_k
  alpha = E/dot_product(b,p)

! ****** Update x, r and z
     x0 = x0 + alpha*p
     r = r - alpha*b
     call invM(A,msize,r,z)

! ****** Calculate r_k^T r_k
     E_new = dot_product(r,z)
               
! ****** Update beta and err
     beta = E_new/E
     E = E_new
     write(*,*) "Error at iteration",k,"= ",E
     
! ****** Update p_k
     p = z + beta*p
     k = k+1   
  enddo

end subroutine precondconjgrad

end module LinAl
