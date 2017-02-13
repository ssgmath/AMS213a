program secondprog
  !***** PROGRAM TO PERFORM QR decomposition of a (non-square) matrix ************
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -o prog ***********************
  !***** TO RUN: :/prog atkinson.dat *********************************************

  !***** uses module LinAl that includes QRdec and backsub subroutines ***********

  use LinAl

  real,dimension(:,:), allocatable :: A,As,Q,QT,R
  real,dimension(:,:), allocatable :: I_msize
  real,dimension(:),allocatable :: x, b, bs, error, E ! vectors to store the x and y values to solve the linear problem Ax=b, and find error in solution
  logical :: ising=.false.
  integer :: msize, nsize, i, j
  real :: err, f ! to store error between f value at given x and computed solution value
  character*100 filename,outfilename


  ! **********************************************************
  !
  ! This program calculates the best fit polynomial curve 
  ! f(x) = a_0 + a_1 x + a_2 x^2 + ... a_(n-1) x^(n-1)
  ! to a dataset given in a user defined file (atkinson.dat) 
  ! based on method of Least Squares with QR decomposition
  !
  ! **********************************************************

  if(iargc().ne.1) then 
     write(*,*) 'Wrong number of arguments (need file name)'
     stop
  endif
  call getarg(1,filename)

  ! *** input order of fitting polynomial

  write(*,*) 'order of fitting polynomial' 
  read(*,*) nsize
  nsize = nsize + 1

  ! *** Reads the data from file ******************************

  ! count how many lines the file has
  i = 0
  io = 0
  open(10,file=filename)
  do while(io.eq.0) 
     i = i + 1
     read(10,*,iostat=io) temp,temp
  enddo
  close(10)
  msize = i-1

  ! allocate arrays and vectors for the ls problem
  allocate(A(msize,nsize))
  allocate(As(msize,nsize))
  allocate(I_msize(msize,msize))
  allocate(b(msize))
  allocate(x(msize))

  ! red in the vectors x and b
  open(10,file=filename)
  do i = 1,msize 
     read(10,*) x(i),b(i)
  enddo
  close(10)

  write(*,*) 'The data file has',msize,' points'

  ! *** Constructs the matrix A **************! 

  do i=1,msize
     do j = 1,nsize
        A(i,j) = x(i)**(j-1)
     enddo
  enddo

  ! Copy A into As and b into bs
  As = A
  bs = b

  !write out original form of A and B onto screen
  !write(*,*) 'The matrix A is:'
  !call writemat(As,msize,nsize)
  !write(*,*) 'The matrix B is:'
  !call writemat(B,msize,nsize)
  !call subroutine to reduce A into upper triangular form    
  allocate(Q(msize,msize))
  call QRdec(A,Q,msize,nsize)
  !output A in triangular form
  write(*,*) 'The matrix A after QR decomposition is:'
  call writemat(A,msize,nsize)
  !Construct QT containing first n columns of Q
  allocate(QT(msize,nsize))
  do i=1,msize
     QT(i,1:nsize)=Q(i,1:nsize)
  end do
  !Construct R containing first n rows of A
  !P_A=matmul(Q^,Q^T)
  allocate(R(nsize,nsize))
  do i=1,nsize
     R(i,1:nsize)=A(i,1:nsize)
  end do
  !solution of Rx=Q^T*b by backsustitution
  b = matmul(transpose(QT),b)
  call backsub(R,b,nsize)
  !output the solution onto screen
  write(*,*) 'The solution vector is:'
  !call writemat(B,msize,nsize)
  write(*,*) (b(i), i=1,nsize)
  ! compute error matrix
  allocate(E(nsize))
  E = bs - matmul(As,b)
  write(*,*)'error in solution is',norm2(E,nsize)
  !compute rms error bewtween fitted curve and data points given by E²=sum(b(i)-f(x(i),bs))²
  allocate(error(msize))
  err = 0.0
  open(20,file='fitted.dat')
  do i=1,msize
     f = 0.0
     do j = 1, nsize
        f = f + As(i,j)*b(j)
     end do
     write(20,*) x(i), f
     error(i) = bs(i) - f
     err = err + error(i)**2
  enddo
  close(20)
  write(*,*)'the rms error between fitted curve and data is given by', sqrt(err)

  write(*,*) 'The matrix A-QR is:'
  !call writemat(As-matmul(Q,A),msize,nsize)
  write(*,*)'Frobenius nom of A-QR=',normF(As-matmul(Q,A),msize,nsize)
  do j=1,msize
     do i=1,msize
        I_msize(i,j) = (i/j)*(j/i)    ! integer arithmetic!!! 
     end do
  end do
  write(*,*) 'The matrix Q^T*Q - I is:'
  !call writemat(matmul(transpose(Q),Q)-I_msize,msize,msize)
  write(*,*)'Frobenius norm of Q^T*Q - I =',normF(matmul(transpose(Q),Q)-I_msize,msize,msize)
  deallocate(A,As,Q,QT,R,b,I_msize,E)
end program secondprog
