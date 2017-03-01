program thirdprog
  !***** PROGRAM TO find eigenvalues of a Lehmer matrix of size msize  ****
  !***** COMPILE as: gfortran LinAl.f90 myprog.f90 -lblas -llapack -o prog ***********************
  !***** TO RUN: :./prog <msize> <0/1> (0:using QR routine, 1: using LAPACK routine) **************************************************
  use LinAl
  real,dimension(:,:), allocatable :: mat, mmat, Q, V, L_msize
  real,dimension(:),allocatable :: colvec, WORK, eigenvec, error
  real :: x, sum = 0., err=1.0, start, finish !to store trace, sum of eigenvalues & measure code-runtime
  integer :: msize, nsize, n_method, INFO, iter = 0 !matrix dimensions, LAPACK flag & iteration count
  logical :: flag =.TRUE.
  character(len=100) :: m_size, method, filename
  !**** input the matrix size as argumentfrom user **********
  if(iargc().ne.2) then
     write(*,*) 'Wrong number of arguments (need size of Lehmer matrix and metid to use:0 for QR/1 for LAPACK)'
     stop
  endif
  call getarg(1,m_size)
  call getarg(2,method)
  read(m_size,*)msize
  read(method,*)n_method
  !********** set flag to switch between user's choice for subroutine or LAPACK ***********
  if(n_method.eq.1)  flag = .FALSE.
  
  !**********call subroutine to read and allocate the matrix**********
  !call readandallocatemat(mat,msize,nsize,filename)
  allocate(Q(msize,msize))
  allocate(V(msize,msize))
  allocate(error(msize))
  allocate(L_msize(msize,msize))
  allocate(colvec(msize))
  allocate(eigenvec(msize))
  allocate(WORK(3*msize))

  !********** call subroutine lehmer_mat *********************!
  call lehmer_mat(L_msize,msize)

  !*********calculate and output its trace****************************
     x = trace(L_msize,msize)
     write(*,*) 'The trace of this matrix is ', x
  
  !********* create a copy of mat for later use ********!
  mmat = L_msize
  !********* initialise V to identity ******************!
  do j=1,msize
     do i=1,msize
        V(i,j) = (i/j)*(j/i) ! integer arithmetic!!! 
     end do
  end do
  
  !****** initialise the eigenvalues ************!
  colvec = 0
  !******call QR algorithm to find eigenvalues and eigenvectors  *****! 
  if(flag)then
     call cpu_time(start)
     call QRalgo(L_msize,colvec,msize,msize)
     call cpu_time(finish)
     !write(*,*) "The matrix of eigenvectors is given by"
     !call writemat(V,msize,msize)
     
  else
     call cpu_time(start)
     !** call LAPACK routine to compare **************!
     call SSYEV('N','U',msize,L_msize,msize,colvec,WORK,3*msize,INFO)
     call cpu_time(finish)
     !write(*,*) 'The matrix of eigenvectors given by LAPACK:'
     !call writemat(L_msize,msize,msize)
  endif
  !********* sort the eigenvalues in ascending order to compare with LAPACK values ***
  call heapsort(colvec, msize)
  !*********calculate sum of eigenvalues which should equal trace*****
  if(flag)then
  filename = 'evalues_code_'//trim(m_size)//'.dat'
  else
  filename = 'evalues_LAPACK_'//trim(m_size)//'.dat'
  endif
  open(10,file=filename)
  do j=1,msize
     do i=1,msize
        if(flag)then    
        eigenvec(i) = V(i,j)
        else
        eigenvec(i) = L_msize(i,j)
        endif
     end do
     write(10,*) j,colvec(j)
     !write(*,*)'|Av-lambda*v| =', norm2(matmul(mmat,eigenvec) - colvec(j)*eigenvec,msize)
     sum = sum + colvec(j)
  end do
  close(10)
  !write(*,*) 'The norm of column',j,'=',norm2(colvec,msize)
  write(*,*)'the sum of eigenvalues =', sum
  deallocate(colvec,eigenvec,WORK,error)
  deallocate(L_msize,mmat,Q,V)
  write (*,*) '("Time = ",f6.3," seconds.")',finish-start
1000 write(*,*) '**** END OF PROGRAM ************'
end program thirdprog
