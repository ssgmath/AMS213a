program firstprog 
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -o prog ***********************
  !***** TO RUN: :/prog Amat.dat Bmat.dat ****************************************

  !****uses module LinAl that includes QRdec and backsub subroutines**

  use LinAl
  real,dimension(:,:), allocatable :: A,As,ATA,Q,QT,R
  real,dimension(:,:), allocatable :: I_msize
  real,dimension(:), allocatable :: x, b, bs, ATb, error, E ! vectors to store the x and y values to solve the linear problem Ax=b, and find error in solution
  logical :: ipos =.FALSE., ising =.FALSE.
  integer :: msize, nsize, i, j
  real :: err, f ! to store error between f value at given x and computed solution value
  character*100 filename,outfilename


  ! ****************************************************************
  !
  ! This program calculates the best fit polynomial curve 
  ! f(x) = a_0 + a_1 x + a_2 x^2 + ... a_(n-1) x^(n-1)
  ! to a dataset given in a user defined file (atkinson.dat) 
  ! based on method of Least Squares using Cholesky decomposition
  ! The output for the solution contains the coefficients a_{j}
  ! the error in solution given by |b-A*x| as well as rms error 
  ! between the points and fitted curve is outputted to screen
  !
  ! ****************************************************************

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
  allocate(b(msize))
  allocate(x(msize))
  allocate(E(nsize))
  allocate(error(msize))

  ! red in the vectors x and b
  open(10,file=filename)
  do i = 1,msize 
     read(10,*) x(i),b(i)
  enddo
  close(10)

  write(*,*) 'The data file has',msize,' points'

  ! *** Constructs the matrix A ************** !

  do i = 1,msize
     do j = 1,nsize
        A(i,j) = x(i)**(j-1)
     enddo
  enddo

  ! Copy A into As and b into bs
  As = A
  bs = b

  ! compute transpose if A^T*A and A^T*b
  allocate(ATA(nsize,nsize))
  allocate(ATb(nsize))
  ATA = matmul(transpose(A),A)
  ATb = matmul(transpose(A),b)
  !call writemat(ATA,nsize,nsize)
  !call subroutine to reduce A into lower triangular form   
  call Choldec(ATA,nsize,ipos,ising)
  !stop if either ipos or ising is true
  if(ipos.or.ising) goto 1000
  !solution of L*L^T*x=b by backsustitution
  call Cholsol(ATA,nsize,ATb)
  !output the solution onto screen
  write(*,*) 'The solution vector is:'
  !call writemat(B,msize,nsize)
  write(*,*) (ATb(i), i=1,nsize)
  ! compute error matrix
  E = bs - matmul(As,Atb)
  write(*,*)'error in solution is',norm2(E,nsize)
  !compute rms error bewtween fitted curve and data points given by E²=1/msize*sum(b(i)-f(x(i),bs))²

  err = 0.0
  open(20,file='fitted.dat')
  do i=1,msize
     f = 0.0
     do j = 1, nsize
     f = f + As(i,j)*ATb(j)
     end do
     write(20,*) x(i), f
     error(i) = bs(i) - f
     err = err + error(i)**2
  enddo
  close(20)
  write(*,*)'the rms error between fitted curve and data is given by', sqrt(err/msize)
1000 write(*,*)'******* END ********'
  deallocate(A,As,ATA,ATb,b,E,error)
end program firstprog
