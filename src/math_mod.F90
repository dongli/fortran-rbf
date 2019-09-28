module math_utils_mod

  use flogger

  implicit none

  private

  public factorial
  public nchoosek
  public create_polynomial
  public inverse_matrix

contains

  pure recursive integer function factorial(x) result(res)

    integer, intent(in) :: x

    if (x == 1) then
      res = 1
    else
      res = x * factorial(x - 1)
    end if

  end function factorial

  pure integer function nchoosek(m, n) result(res)

    integer, intent(in) :: m
    integer, intent(in) :: n

    if (n == 0) then
      res = 1
    else
      res = factorial(m) / factorial(n) / factorial(m - n)
    end if

  end function nchoosek

  recursive subroutine create_polynomial(num_dim, order, poly, col)

    ! num_dim = 2, order = 0:
    !
    !   0
    !
    ! num_dim = 2, order = 2:
    !
    !   0  1  2  1  1  2
    !   0  0  0  1  2  2
    !
    ! num_dim = 2, order = 3:
    !
    !   0  1  2  1  1  2  1  1  1  2
    !   0  0  0  1  2  2  1  1  2  2
    !   0  0  0  0  0  0  1  2  2  2

    integer, intent(in) :: num_dim
    integer, intent(in) :: order
    integer, intent(inout), allocatable :: poly(:,:)
    integer, intent(inout), optional :: col

    integer i, j
    logical done

    if (.not. present(col)) then
      allocate(poly(max(1, order),nchoosek(order + num_dim, order)))
      poly = 0
      j = 2
      do i = 1, order
        poly(1:i,j) = 1
        j = j + 1
        call create_polynomial(num_dim, i, poly, j)
      end do
      ! ! Debug print
      ! print *, shape(poly)
      ! do i = lbound(poly, 1), ubound(poly, 1)
      !   do j = 1, size(poly, 2)
      !     write(*, '(I3)', advance='no') poly(i,j)
      !   end do
      !   write(*, *)
      ! end do
    else
      do while (.true.)
        j = col - 1
        done = .true.
        do i = order, 1, -1
          if (poly(i,j) /= num_dim) then
            done = .false.
            exit
          end if
        end do
        if (done) return
        poly(:,col) = poly(:,j)
        poly(i:order,col) = poly(i,col) + 1
        col = col + 1
      end do
    end if

  end subroutine create_polynomial

  subroutine inverse_matrix(A, Ai)

    real(8), intent(in) :: A(:,:)
    real(8), intent(out) :: Ai(:,:)

#ifdef USE_LAPACK
    real(8) work(size(A, 1))

    integer ipiv(size(A, 1)), ierr
    
    external dgetrf ! LAPACK factorization subroutine
    external dgetri ! LAPACK inversion subroutine

    Ai = A

    call dgetrf(size(A, 1), size(A, 1), Ai, size(A, 1), ipiv, ierr)
    if (ierr /= 0) call log_error('Matrix is numerically singular!', __FILE__, __LINE__)

    call dgetri(size(A, 1), Ai, size(A, 1), ipiv, work, size(A, 1), ierr)
    if (ierr /= 0) call log_error('math_inv_matrix: Matrix inversion failed!', __FILE__, __LINE__)
#else
    
#endif


  end subroutine inverse_matrix

end module math_utils_mod