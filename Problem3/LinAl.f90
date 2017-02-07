module LinAl
  implicit none

contains

  !********************************************************

  subroutine readandallocatemat(mat,msize,nsize,filename)

    character*100 filename
    double precision, dimension(:,:), allocatable, intent(in out) :: mat
    integer :: msize, nsize

    integer :: i,j

    ! Reads a file containing the matrix mat
    ! Note that the first 2 lines are the matrix dimensions, 
    ! then the next msize lines are the matrix entries
    ! Note that entries must be separated by a tab.
    ! Then allocates an array of size (msize,nsize), populates the matrix,
    ! and returns the array. 

    ! This routine takes as INPUTS:
    ! filename = a character string with the name of file to read
    ! This routine returns as OUTPUTS:
    ! msize = first dimension of read matrix
    ! nsize = second dimension of read matrix
    ! mat = read matrix. Note that mat type is real.

    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) msize,nsize

    ! Allocate matrix
    allocate(mat(msize,nsize))

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)

  end subroutine readandallocatemat

  real function trace(mat,msize)
    !function to calculate trace of a square matrix mat of size msize
    integer :: msize, i
    double precision, dimension(msize,msize) :: mat
    trace = 0.
    do i=1,msize
       trace = trace + mat(i,i)
    enddo
  end function trace

  real function norm2(vec,vsize)
    !function to compute the (Euclidean)norm of a vector vec of size vsize
    integer :: vsize, i
    double precision, dimension(vsize) :: vec
    norm2 = 0.
    do i=1,vsize
       norm2 = norm2 + vec(i)**2
       !write(*,*) norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in rowwise split format 
    integer :: msize, nsize, i, j
    double precision, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat

  subroutine lufactor(A,msize,s,ising)
    !****Factors the matrix A of dimension (msize,msize) *******************
    !****into L,U .. On return L-I  is stored in the lower triangular part**
    !****of A and U  is stored in the upper triangular part. The vector s***
    !**** stores the permuted row indices using partial pivoting ***********

    integer :: msize,i,j,k,colmax,srow,index,temp
    double precision, dimension(msize,msize) :: A
    double precision :: p,dum
    integer, dimension(msize) :: s
    logical :: ising

    do j=1,msize
       s(j)=j
    enddo
    !Loop over the columns
    do j=1,msize
       !Find the pivot row index and pivot p=max|A(k,j)|
       p=0.
       do k=j,msize
          srow=dabs(A(k,j))
          if(p.lt.srow)then
             p=srow
             index=k
          end if
       enddo
       !write(*,*)'the pivot is',p,index
       !exchange rows if index is not same as j
       if(index.ne.j) then
          do k=1,msize   ! swap rows in A 
             dum = A(index,k)
             A(index,k) = A(j,k)
             A(j,k) = dum
          enddo
          !exchange entries of s
          temp=s(j)
          s(j)=index
          s(index)=temp
       endif
       !check for zero entry and stop with flag
        if(A(j,j).eq.0.0) then
          write(*,*) 'A is singular'
          ising = .TRUE.
          goto 1000
       endif
       !Calculate the i th column of L
       do i=j+1,msize
          A(i,j)= A(i,j)/A(j,j)
          do k=j+1,msize
             A(i,k)=A(i,k) - A(i,j)*A(j,k)
          end do
       end do
    end do
1000 continue

  end subroutine lufactor

  subroutine backsub(A,B,msize,nsize,s)
    double precision, dimension(msize,msize) :: A
    double precision, dimension(msize,nsize) :: B,Bperm
    integer, dimension(msize) :: s
    real :: sum
    integer :: i,j,k,l,msize,nsize


    !*********** perfoms backsubstitution to solve Ax=B*******
    !*********** with B of dimension (misze,nsize) and square* 
    !*********** matrix A(msize,msize) is in LU form obtained*
    !*********** as output from lufactor subroutine with the**
    !**********  corresponding permution vector s ************

    !*********** y=PB and store in B****!
    do j=1,nsize
       l=0
       do i=1,msize
          Bperm(i,j) = B(s(i),j)
          !write(*,*) s(i)
       end do
       B(:,j) = Bperm(:,j)

       !*********create y=L^-1*P*B and store in B **********!
       do i=1,msize
          sum = B(i,j)
          if (l.ne.0) then
             do k=l,i-1
                sum=sum-A(i,k)*B(k,j) 
             enddo
          else if (sum.ne.0.) then
             l=i
          endif
          B(i,j)=sum         
       enddo

       !*****standard backsubtitution*********!
       do i=msize,1,-1
          sum = 0.0
          if(A(i,i).eq.0.0)then
1000         write(*,*)'singular matrix'
             goto  2000
          end if
          do k=i+1,msize
             sum = sum + A(i,k)*B(k,j)
          end do
          B(i,j) = (B(i,j) - sum)/A(i,i)
       end do
    end do

2000 return
  end subroutine backsub

end module LinAl
