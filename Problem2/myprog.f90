program secondprog
  !***** PROGRAM TO SOLVE AX=B using Gaussian elimination with partial pivoting***
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -o prog ***********************
  !***** TO RUN: :/prog Amat.dat Bmat.dat ****************************************

  !****uses module LinAl that includes gausse and backsub subroutines*
  !****to solve the system of linear equation Ax=B ******************* 
  use LinAl
  real,dimension(:,:), allocatable :: A,B,As,Bs,E
  real,dimension(:),allocatable :: colvec !to find norms of columns of A
  logical :: ising
  integer :: msize, nsize
  character*100 filenameA, filenameB
  if(iargc().ne.2) then
     write(*,*) 'Wrong number of arguments (need file names for matrices A & B)'
     stop
  endif
  call getarg(1,filenameA)
  call getarg(2,filenameB)
  call readandallocatemat(A,mize,msize,filenameA)
  call readandallocatemat(B,msize,nsize,filenameB)
  allocate(As(msize,msize))
  allocate(Bs(msize,nsize))
  !save copies of A and B in As and Bs
  As=A
  Bs=B
  !write out original form of A and B onto screen
  write(*,*) 'The matrix A is:'
  call writemat(A,msize,msize)
  write(*,*) 'The matrix B is:'
  call writemat(B,msize,nsize)
  !call subroutine to reduce A into upper triangular form
  call gausse(A,B,msize,nsize,ising)
  !output A in triangular form
  write(*,*) 'The matrix A after gauss elimination is:'
  call writemat(A,msize,msize)
  !solution of Ax=B by backsustitution
  call backsub(A,B,msize,nsize)
  !output the solution onto screen
  write(*,*) 'The solution vector is:'
  call writemat(B,msize,nsize)
  !allocate and compute error matrix s E=As*B - Bs
  allocate(E(msize,nsize))
  E = matmul(As,B) - Bs
  write(*,*) 'The error matrix is:'
  call writemat(E,msize,nsize)
  !find norms of the columns of E
  allocate(colvec(msize))
  do j=1,nsize
     do i=1,msize
        colvec(i)= E(i,j)
        !write(*,*) colvec(i)
     end do
     write(*,*) 'The eucleidean norm of column',j,' is:',norm2(colvec,msize)
  end do
  deallocate(colvec,A,B,As,Bs,E)
end program secondprog
