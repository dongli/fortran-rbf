module resize_array_mod

  implicit none

  private

  public resize_array

  interface resize_array
    module procedure resize_array_2d_r8
  end interface resize_array

contains

  subroutine resize_array_2d_r8(x, dim, new_size)

    real(8), intent(inout), allocatable :: x(:,:)
    integer, intent(in) :: dim(:)
    integer, intent(in) :: new_size(:)

    real(8), allocatable :: y(:,:)
    integer n(2), i

    n = shape(x)

    do i = 1, size(dim)
      n(dim(i)) = new_size(i)
    end do

    allocate(y(n(1),n(2)))

    y(:n(1),:n(2)) = x(:n(1),:n(2))

    deallocate(x)

    allocate(x(n(1),n(2)))

    x = y

    deallocate(y)

  end subroutine resize_array_2d_r8

end module resize_array_mod