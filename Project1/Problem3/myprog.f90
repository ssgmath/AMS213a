program thirdprog
!***** PROGRAM TO SOLVE AX=B using LU decomposition with partial pivoting*
!****COMPILE as: gfortran -fdefault-real-8 LinAl.f90 myprog.f90 -o prog***
!****TO RUN: :/prog Amat.dat Bmat.dat ************************************

!****uses module LinAl that includes lufactor and backsub subroutines*
!****to solve the system of linear equation Ax=B *********************
  use LinAl
  double precision,dimension(:,:), allocatable :: A,B,As,Bs,E
  double precision,dimension(:), allocatable :: colvec !to find norms of columns of A
  integer,dimension(:), allocatable :: s
  logical :: ising=.FALSE.
  integer :: msize, nsize, d, int
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
  !allocate permutation vector s
  allocate(s(msize))
  !call subroutine to factor A into L and U
  call lufactor(A,msize,s,ising)
  !write out A in LU form with L-I stored in lower triangular part
  !and U stored in upper triangular part of A
  write(*,*) 'The matrix A after LU decomposition is:'
  call writemat(A,msize,msize)
  !solution of Ax=B by backsustitution
  call backsub(A,B,msize,nsize,s)
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
     end do
     write(*,*) 'The eucleidean norm of column',j,' is:',norm2(colvec,msize)
  end do
  deallocate(colvec,s,A,B,As,Bs)
end program thirdprog
