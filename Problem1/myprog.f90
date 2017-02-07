program firstprog
  !***** PROGRAM TO calculate trace and norm of columns of a square matrix mat****
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -o prog ***********************
  !***** TO RUN: :/prog Amat.dat**************************************************
  use LinAl
  real,dimension(:,:), allocatable :: mat
  real,dimension(:),allocatable :: colvec !to find norms of columns of mat
  real :: x !to store trace
  integer :: msize, nsize !matrix dimensions
  character*100 filename !filehandle for input file containing matrix
  !**** input the matrix filename as first argument *****************
  if(iargc().ne.1) then
     write(*,*) 'Wrong number of arguments (need file name)'
     stop
  endif
  call getarg(1,filename)
  !**********call subroutine to read and allocate the matrix**********
  call readandallocatemat(mat,msize,nsize,filename)
  allocate(colvec(nsize))
  !*********output the matrix onto screen*****************************
  write(*,*) 'The input matrix is:'
  !*********calculate and output its trace****************************
  call writemat(mat,msize,nsize)
  if(nsize.ne.msize) then
     write(*,*) 'This is not a square matrix, cannot calculate trace'
  else
     x = trace(mat,msize)
     write(*,*) 'The trace of this matrix is ', x
  end if
  !*********calculate norm of the column vectors of input matrix*******
  do j=1,nsize
     do i=1,msize
        colvec(i)= mat(i,j)
        !write(*,*) colvec(i)
     end do
     write(*,*) 'The norm of column',j,'=',norm2(colvec,msize)
  end do
  deallocate(colvec)
  deallocate(mat)
end program firstprog
