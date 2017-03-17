program firstprog
  !***** PROGRAM TO solve the linear system Ax = b using Conjugate gradient method without and with (diagonal) preconditioner****
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -lblas -llapack -o prog ***********************
  !***** TO RUN: :./prog 1/2 (1:Conjugate gradient without preconditioning, 2: with diagonal preconditioning)************************!
  use LinAl
  real,dimension(:,:), allocatable :: mat
  real,dimension(:),allocatable :: b, x0
  integer :: msize,D ! size of the matrix A & its diagonal entries
  character*100 filename !filehandle for input file containing matrix
  real :: acc = 1d-6 ! user defined accuracy for solution of Ax = b
  character*100 method ! method to use (with or without preconditioning)
 !**** input the method to use as argument(1/2) *****************
  if(iargc().ne.1) then
  write(*,*) 'Wrong number of arguments (need method to use)'
  stop
  endif
  call getarg(1,method)
  read(method,*)n_method

  ! ask user for size of matrix and diagonal entry
  write(*,*)'Input size of matrix A'
  read(*,*) msize
  write(*,*)'Enter diagonal element D'
  read(*,*) D
  allocate(mat(msize,msize))
  allocate(b(msize))
  allocate(x0(msize))
  
  !Set all of mat to 1 
  mat = 1.
  !Set i-th diagonal element to D or i
  do i=1,msize
     b(i) = real(i)
     do j=1,msize
        if(i.eq.j) mat(i,j) = real(i) !real(D)
     end do
  enddo

  !*********set RHS vector b **********************!
  do i=1,msize
  b(i) = real(i)
  end do
  !*********output the matrix onto screen*****************************
  !write(*,*) 'The input matrix is:'
  !call writemat(mat,msize,msize)

  !******call Conjugate gradient method with or without preconditioning*****!
  ! Initial guess for solution
  x0 = 0.
  if(n_method.eq.1) then
  call ConjGrad(mat,b,x0,msize,acc)
  write (*,*)'the solution using conjugate gradient is given by',(x0(i),i=1,msize)
  else
  call precondConjGrad(mat,b,x0,msize,acc,invM)
  write (*,*)'the solution using preconditioned conjugate gradient is given by',(x0(i),i=1,msize)
  endif
  deallocate(b,x0)
  deallocate(mat)
1000 write(*,*) '**** END OF PROGRAM ************'
end program firstprog
