program firstprog
  !***** PROGRAM TO find eigenvalues and eigenvectors of a real symmetric matrix ****
 !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -lblas -llapack -o prog ***********************
  !***** TO RUN: :./prog A.dat/Emat1.dat/Emat2.dat **************************************************
  use LinAl
  real,dimension(:,:), allocatable :: mat, mmat, Q, V, I_msize
  real,dimension(:),allocatable :: colvec, WORK, eigenvec, error
  real :: x, sum = 0., err=1.0, ABNRM !to store trace and sum of eigenvalues
  integer :: msize, nsize, INFO, iter = 0 !matrix dimensions, LAPACK flag & iteration count
  character*100 filename !filehandle for input file containing matrix
  logical :: flag = .TRUE.
  !**** input the matrix filename as first argument *****************
  if(iargc().ne.1) then
     write(*,*) 'Wrong number of arguments (need file name)'
     stop
  endif
  call getarg(1,filename)
  !**********call subroutine to read and allocate the matrix**********
  call readandallocatemat(mat,msize,nsize,filename)
  allocate(mmat(msize,nsize))
  allocate(Q(msize,msize))
  allocate(V(msize,msize))
  allocate(error(msize))
  allocate(I_msize(msize,msize))
  allocate(colvec(msize))
  allocate(eigenvec(msize))

  !*********output the matrix onto screen*****************************
  write(*,*) 'The input matrix is:'
  !*********calculate and output its trace****************************
  call writemat(mat,msize,nsize)
  if(nsize.ne.msize) then
     write(*,*) 'This is not a square matrix, cannot calulate trace' 
     goto 1000
  else
     x = trace(mat,msize)
     write(*,*) 'The trace of this matrix is ', x
  end if
  !********* create a copy of mat for later use ********!
  mmat = mat
  !**** flag to check for real symmetric matrix *************!
  do i=1,msize
     do j=1,nsize
        if(abs(mat(i,j)-mat(j,i)).gt.epsilon(0.))then
           flag = .FALSE.
           !write(*,*) abs(mat(i,j)-mat(j,i))
        endif
     end do
  end do
  !***** intialise V, I_msize to identity square matrix of size msize **!
  do j=1,msize
     do i=1,msize 
        I_msize(i,j) = (i/j)*(j/i)    ! integer arithmetic!!! 
     end do
  end do
  V = I_msize
  eigenvec = 0
  !******call QR algorithm and iterate till error remains too large *****! 
  if(flag)then
     call QRalgo(mat,V,colvec,msize,nsize)
     write(*,*) "The matrix of eigenvectors is given by"
     call writemat(V,msize,nsize)
  else
     write(*,*)'**** WARNING **** matrix is not symmetric'  
     call QRalgo(mat,V,colvec,msize,nsize)     
  endif

  !*********calculate sum of eigenvalues which should equal trace*****
  do j=1,msize
     do i=1,msize
        !colvec(i)= mmat(i,j)     
        eigenvec(i) = V(i,j)     
     end do
     if(flag)then
     write(*,*)'for lambda=',colvec(j)
     write(*,*)'|Av-lambda*v| =', norm2(matmul(mmat,eigenvec) - colvec(j)*eigenvec,msize)
     endif
     sum = sum + colvec(j)
  end do
  !write(*,*) 'The norm of column',j,'=',norm2(colvec,msize)
  write(*,*)'the sum of eigenvalues =', sum
  deallocate(colvec,eigenvec,error)
  deallocate(mat,mmat,Q,V)
1000 write(*,*) '**** END OF PROGRAM ************'
end program firstprog
