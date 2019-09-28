module distance_mod

  ! This module contains functions that calculate distance between nodes.
  
  use const_mod, only: r8

  implicit none

  private

  public distance
  public distance_choose
  public euclidean_distance

  integer, parameter :: euclidean_distance = 1

  interface
    pure real(r8) function distance_interface(x, y)
      import r8
      real(r8), intent(in) :: x(:)
      real(r8), intent(in) :: y(:)
    end function distance_interface
  end interface

  procedure(distance_interface), pointer :: distance

contains

  subroutine distance_choose(type)

    integer, intent(in) :: type

    select case (type)
    case (euclidean_distance)
      distance => calc_euclidean_distance
    end select

  end subroutine distance_choose

  pure real(r8) function calc_euclidean_distance(x, y) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: y(:)

    res = norm2(x - y)

  end function calc_euclidean_distance

end module distance_mod