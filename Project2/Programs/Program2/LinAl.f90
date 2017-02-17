module LinAl
  implicit none

contains

  real function norm2(vec,vsize)
    !function to compute the (Euclidean)norm of a vector vec of size vsize  !outputs 
    integer :: vsize, i
    real, dimension(vsize) :: vec
    norm2 = 0.
    do i=1,vsize
       norm2 = norm2 + vec(i)**2
       !write(*,*) norm2
    enddo
    norm2 = sqrt(norm2)
  end function norm2

  real function normF(A,msize,nsize)
    !function to compute the Frobenius norm of a matrix of size dimensions (msize,nsize) 
    integer :: msize, nsize,i,j
    real, dimension(msize,nsize) :: A
    normF = 0.
    do i=1,msize
       do j=1,nsize
          normF = normF + A(i,j)**2
          !write(*,*) norm2
       enddo
    enddo
    normF = sqrt(normF)
  end function normF

  subroutine writemat(mat,msize,nsize)
    !outputs matrix mat of dimension (msize,nize) in row-wise split format
    integer :: msize, nsize, i, j
    real, dimension(msize,nsize) :: mat
    do i=1,msize
       write(*,*) ( mat(i,j), j=1,nsize )
    enddo
  end subroutine writemat

  subroutine QRdec(A,Q,msize,nsize)

    integer :: msize, nsize
    real, dimension(msize,nsize) :: A
    real, dimension(msize,msize) :: Q, H, I_msize
    real, dimension(msize) :: v
    integer :: i,j,k,imax       ! locally defined variables
    real:: sum,sgn

    !***** intialise Q, I_msize to identity square matrix of size msize ****!
    do j=1,msize 
       do i=1,msize 
          I_msize(i,j) = (i/j)*(j/i)    ! integer arithmetic!!! 
       end do
    end do
    Q = I_msize
    !***** Performs  Hausholder QR decomposition of a given matrix**! 
    ! **** Inputs are:
    ! A = matrix of dimensions (msize,nsize)
    ! B = matrix of RHS vectors, dimensions are B(msize,nsize) 
    ! ising = logical flag to check if A is singular

    !**** loop over columns initialising sum and v to 0 *********!

    do j=1,nsize
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


  subroutine backsub(A,b,nsize)
    integer :: i,j,k,nsize
    real, dimension(nsize,nsize) :: A
    real, dimension(nsize) :: b, x
    real :: sum

    do i=nsize,1,-1
       sum = 0.     
       do k=i+1,nsize
          sum = sum + A(i,k)*b(k)
       end do
       b(i) = (b(i) - sum)/A(i,i)
    end do
  end subroutine backsub

end module LinAl
