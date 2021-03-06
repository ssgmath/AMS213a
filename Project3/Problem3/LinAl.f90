module LinAl
  implicit none

contains

  !********************************************************

  subroutine readandallocatemat(mat,msize,nsize,filename)

    character*100 filename
    real, dimension(:,:), allocatable, intent(in out) :: mat
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
  
  subroutine lehmer_mat(L_msize,msize)
  real, dimension(msize,msize):: L_msize
  integer :: i,j,msize  

  !***** intialise V, I_msize to identity square matrix of size msize **!
  do j=1,msize
     do i=1,msize
        L_msize(i,j) = real(min(i,j))/real(max(i,j)) ! NOT integer arithmetic!!!
     end do
  end do
  end subroutine lehmer_mat


  real function trace(mat,msize)
    !function to calculate trace of a square matrix mat of size msize
    integer :: msize, i
    real, dimension(msize,msize) :: mat
    trace = 0.
    do i=1,msize
       trace = trace + mat(i,i)
    enddo
  end function trace

  real function norm2(vec,vsize)
    !function to compute the (Euclidean)norm of a vector vec of size vsize
    integer :: vsize, i
    real, dimension(vsize) :: vec
    norm2 = 0.
    do i=1,vsize
       norm2 = norm2 + vec(i)**2
       !write(*,*) norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in row-wise split format 
    integer :: msize, nsize, i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat

  subroutine QRdec(A,Q,msize,nsize)

    !***** Performs  Hausholder QR decomposition of a given matrix**! 
    ! **** Inputs are:
    ! A = matrix of dimensions (msize,nsize)
    ! B = matrix of RHS vectors, dimensions are B(msize,nsize) 
    ! ising = logical flag to check if A is singular

    !**** loop over columns initialising sum and v to 0 *********!
    integer :: msize, nsize
    real, dimension(msize,nsize) :: A
    real, dimension(msize,msize) :: Q, H, I_msize
    real, dimension(msize) :: v
    integer :: i,j,k,imax       ! locally defined variables
    real:: sum,sgn

    !***** intialise Q, I_msize to identity square matrix of size msize ****
    !
    do j=1,msize
       do i=1,msize 
          I_msize(i,j) = (i/j)*(j/i)    ! integer arithmetic!!! 
       end do
    end do
    Q = I_msize

    do j=1,nsize-1
       sum = 0.0
       v = 0.0
       !H = 0.0
       !sgn = 0.0
       !*** loop over rows to the right of diagonal of A to find sum of A(j,k)^2 ****!
       do k=j,msize
          sum = sum + A(k,j)**2
       end do
       !*** sgn = sign(A(j,j)*sqrt(sum) using Fortran function sign ***
       sgn = sign(sqrt(sum),A(j,j))
       !*** compute v vectors as (0,0,.....,A(j,j)+sgn,A(j+1,j),..,A(msize,k) ****!
       do k=j,msize
          if(k.eq.j)then
             v(k) = A(j,j) + sgn
          else
             v(k) = A(k,j)
          endif
       end do
       !write(*,*) (v(i), i=1,msize)
       !*** normalize v ****!
       v = v/norm2(v,msize)

       !******* compute the Housholder matrix H = I - 2*v*v^T ********!
       H = I_msize - 2.0*spread(v(1:msize),dim=2,ncopies=msize)*spread(v(1:msize),dim=1,ncopies=msize)

       !******* store R = H*A in A ****************!
       A = matmul(H,A)
       !call writemat(A,msize,nsize)
       !******* calcuate Q as product of the householder matrices**!
       Q = matmul(Q,H)
    end do
    !****** output Q onto screen ***************!
    !write(*,*)'the Q-matrix is given by:'
    !call writemat(Q,msize,msize)
  end subroutine QRdec

  subroutine QRalgo(A,eigenval,msize,nsize)

    !******call QR algorithm and iterate till error remains too large *****! 
    ! **** Inputs are:
    ! A = matrix of dimensions (msize,nsize)
    ! V = matrix of eigenvectors 
    ! ising = logical flag to check if A is singular

    !**** loop over columns initialising sum and v to 0 *********!
    integer :: msize, nsize
    real, dimension(msize,nsize) :: A
    real, dimension(msize,msize) :: Q !, I_msize, V
    real, dimension(msize) :: eigenval, error
    integer :: i,j,iter = 0   ! locally defined variables
    real:: err = 1.0
    logical :: ising

    !***** intialise V, I_msize to identity square matrix of size msize ****
    !
    !do j=1,msize
    !do i=1,msize 
    !I_msize(i,j) = (i/j)*(j/i)    ! integer arithmetic!!! 
    !end do
    !end do
    !V = I_msize

    do while(err.gt.epsilon(0.))
       iter = iter + 1
       !write(*,*)'**********for iteration #',iter,'*************'
       call QRdec(A,Q,msize,nsize)
       A = matmul(A,Q)
       !write(*,*) 'D_lambda matrix is:'
       !call writemat(A,msize,nsize)
       !V = matmul(V,Q)

       do i=1,msize
          error(i) = (A(i,i) - eigenval(i))/A(i,i)
          eigenval(i) = A(i,i)          
          !write(*,*) 'error in eigenvalue =',error(i)
       enddo
       err = norm2(error,msize)
    enddo
  end subroutine QRalgo
  SUBROUTINE heapsort(a,nsize)
    integer ::  nsize
    real, dimension(nsize) :: a
    !******** Sorts an array a(1:nsize) into ascending order using the Heapsort algorithm. 
    ! nsize is input; a is replaced on output by its sorted form
    integer :: i,ir,j,l
    real :: temp
    if (nsize.lt.2) return
    ! The index l will be decremented from its initial value down to 1 during the ¡°hiring¡± (heap creation) phase. Once it reaches 1, the index ir will be decremented from its initial value down to 1 during the ¡°retirement-and-promotion¡± (heap selection) phase.
    l=nsize/2+1
    ir=nsize
10  continue
    if(l.gt.1)then !Still in hiring phase.
       l=l-1
       temp=a(l)
    else !In retirement-and-promotion phase.
       temp = a(ir) !Clear a space at end of array.
       a(ir)=a(1) !Retire the top of the heap into it.
       ir=ir-1 !Decrease the size of the corporation.
       if(ir.eq.1)then !Done with the last promotion.
          a(1)=temp !The least competent worker of all!
          return
       endif
    endif
    i=l !Whether in the hiring phase or promotion phase, we here set up to sift down element temp to its proper level.
    j=l+l 
20  if(j.le.ir)then !¡°Do while j.le.ir:¡±
       if(j.lt.ir)then
          if(a(j).lt.a(j+1)) j=j+1 !Compare to the better underling.
       endif
       if(temp.lt.a(j))then !Demote temp
          a(i)=a(j)
          i=j
          j=j+j
       else !Set j to terminate the sift-down.
          j=ir+1
       endif
       goto 20
    endif
    a(i)=temp !Put temp into its slot.
    goto 10
  END SUBROUTINE heapsort
end module LinAl
