program firstprog
  !***** PROGRAM TO solve the linear system Ax = b using GaussJacbi & GaussSeidel methods ****
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -lblas -llapack -o prog ***********************
  !***** TO RUN: :./prog 1/2 (1: Gauss-Jacobi, 2: Gauss-Seidel)******!
  use LinAl
  real,dimension(:,:), allocatable :: mat
  real,dimension(:),allocatable :: b
  integer :: msize, i, j, D ! size of the matrix A
  character*100 method ! method to use
  real :: acc = 1d-5 ! user defined accuracy for solution of Ax = b
  !**** input the method to use as argument(1/2) *****************
  if(iargc().ne.1) then
  write(*,*) 'Wrong number of arguments (need method to use)'
  stop
  endif
  call getarg(1,method)
  read(method,*)n_method
  write(*,*)'Input size of matrix A'
  read(*,*) msize
  write(*,*)'Enter diagonal element D'
  read(*,*) D
  !**********call subroutine to read and allocate the matrix**********
  !call readandallocatemat(mat,msize,nsize,filename)
  allocate(mat(msize,msize))
  allocate(b(msize))
  !Set all of mat to 1 
  mat = 1.
  !Set i-th diagonal element to D or i
  do i=1,msize
     b(i) = real(i)
     do j=1,msize
        if(i.eq.j) mat(i,j) = real(i) !real(D)
     end do
     !write(*,*)b(i)
  end do
  
  !*********output the matrix onto screen*****************************
  !write(*,*) 'The input matrix is:'
  !call writemat(mat,msize,msize)
 
 !******call Gauss Jacobi/ Gauss-Seidel  algorithm and iterate till error remains too large *****! 
  if(n_method.eq.1) then
  call GaussJacobi(mat,b,msize,D,acc)
  else
  call GaussSeidel(mat,b,msize,D,acc)
  endif
  deallocate(b)
  deallocate(mat)
1000 write(*,*) '**** END OF PROGRAM ************'
end program firstprog
