program firstprog
  !***** PROGRAM TO find eigenvalues and eigenvectors of any square matrix using LAPACK routines ****
 !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -lblas -llapack -o prog ***********************
  !***** TO RUN: :./prog Emat1.dat/Emat2.dat **************************************************
  use LinAl
 real,dimension(:,:), allocatable :: mat, mmat, Q, V, I_msize, VL, VR
  real,dimension(:),allocatable :: colvec, WORK, eigenvec, error, WR, WI !to find norms of columns of mat
  real :: x, sum = 0. !to store trace and sum of eigenvalues
  integer :: msize, nsize, INFO !matrix dimensions
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
  allocate(VL(msize,msize))
  allocate(VR(msize,msize))
  allocate(error(msize))
  allocate(I_msize(msize,msize))
  allocate(colvec(msize))
  allocate(eigenvec(msize))
  allocate(WR(msize))
  allocate(WI(msize))
  allocate(WORK(3*msize))
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
  !******** call LAPACK routine to find eigenvalues and eigenvectors ********
  if(flag)then
     call SSYEV('V','U',msize,mat,msize,colvec,WORK,3*msize,INFO)
      write(*,*) "The matrix of eigenvectors is given by"
  call writemat(mat,msize,nsize)  
  else
     write(*,*)'**** WARNING **** matrix is not symmetric'  
     call SGEEV('N', 'V', msize, mat, msize, WR, WI, VL, msize, VR, msize, WORK, 4*msize, INFO)
     write(*,*)'LAPACK D_lambda is given by:'
     call writemat(mat,msize,nsize)
     write(*,*) "The matrix of eigenvectors is given by"
     call writemat(VR,msize,msize) 
  endif
  !write(*,*) INFO
 
  !*********calculate sum of eigenvalues which should equal to trace*******
  do j=1,msize
     do i=1,msize
        if(flag)then
        eigenvec(i) = mat(i,j) 
        else
        colvec(i) = mat(i,j)
        eigenvec(i) = VR(i,j)
        endif
     end do
     write(*,*)'for lambda=',colvec(j)
     write(*,*)'|Av-lambda*v| =', norm2(matmul(mmat,eigenvec) - colvec(j)*eigenvec,msize)
     sum = sum + colvec(j)
  end do
  !write(*,*) 'The norm of column',j,'=',norm2(colvec,msize)
  write(*,*)'the sum of eigenvalues =', sum
  deallocate(colvec,eigenvec,WORK,WR,WI,error)
  deallocate(mat,mmat,Q,V,VL,VR)
1000 write(*,*) '**** END OF PROGRAM ************'
end program firstprog
