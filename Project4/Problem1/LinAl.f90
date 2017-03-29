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

  subroutine GaussJacobi(A,b,msize,D,acc)

    !***** Solves the linear system Ax=b by Gauss-Jacobi method ***!
    !*****  and returns the solution in b**! 
    ! **** Inputs are:
    ! A = matrix of dimensions (msize,msize)
    ! b = vector of length msize 
    ! acc = desired accuracy for the solution of Ax=b

    integer :: msize
    real, dimension(msize,msize) :: A
    real, dimension(msize) :: b, res, x, y
    integer :: i,k,D,iter=1       ! locally defined variables
    real :: sum, acc
    character(len=100) :: D_value, filename, fmt

    fmt = '(I4.4)'

    !***** guess value for x *****!
    x = 0.
    res = 1.
    write(*,*)'using Gauss-Jacobi method'
        !write(*,*)D
    write(D_value,fmt)D
    filename = 'JacobierrD_'//trim(D_value)//'.dat'
    open(20,file=filename)
    do while(norm2(res,msize).gt.acc)
       write(*,*)'at iteration#',iter, "Error is: ", norm2(b-matmul(A,x),msize)

       ! compute y = D^-1*(b-Rx)
       do i=1,msize  
          sum = 0.
          do k=1,msize
             if (k.ne.i) sum =  sum + A(i,k)*x(k) 
             !write(*,*) sum
          end do
          y(i) = 1./A(i,i)*(b(i) - sum)
       end do
       res = y - x
       write(20,*)iter,norm2(b-matmul(A,x),msize)
       x = y
       iter = iter + 1
    enddo
    close(20)
    b = x
    write(*,*)'solution is given by',(b(i),i=1,msize)

  end subroutine GaussJacobi

  subroutine GaussSeidel(A,b,msize,D,acc)

    !***** Solves the linear system Ax=b by Gauss-Seidel method ***!
    !*****  and returns the solution in b**! 
    ! **** Inputs are:
    ! A = matrix of dimensions (msize,msize)
    ! b = RHS vector of length msize 
    ! acc = desired accuracy for the solution of Ax=b

    integer :: msize
    real, dimension(msize,msize) :: A
    real, dimension(msize) :: b, res, x
    integer :: i,k,D,iter=1       ! locally defined variables
    real :: sum, acc   
    character(len=100) :: D_value, filename, fmt

    fmt = '(I4.4)'

    !***** guess value for x *****!
    x = 0.
    res = 1.
    write(*,*)'using Gauss-Seidel method'
    write(D_value,fmt)D
    filename = 'SeidelD_'//trim(D_value)//'.dat'
    open(20,file=filename)
    do while(norm2(res,msize).gt.acc)
       write(*,*)'at iteration#',iter, "Error is: ", norm2(b-matmul(A,x),msize)
       
       ! compute x = D^-1*(b-Rx)
       do i=1,msize
          sum = 0.
          do k=1,msize
             if (k.ne.i) sum =  sum + A(i,k)*x(k) 
             !write(*,*) sum
          end do
          x(i) = 1./A(i,i)*(b(i) - sum)       
       end do       
       res = b - matmul(A,x)
       write(20,*) iter,norm2(res,msize)
       iter = iter + 1
    enddo
    b = x
    write(*,*)'solution is given by',(b(i),i=1,msize)
    close(20)
  end subroutine GaussSeidel

end module LinAl
