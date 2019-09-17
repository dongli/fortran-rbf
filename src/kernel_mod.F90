module kernel_mod

  use const_mod, only: r8

  implicit none

  private

  public kernel_ga_eval
  public kernel_ga_deriv
  public kernel_phs_eval
  public kernel_phs_deriv

  integer, parameter :: kernel_ga  = 1
  integer, parameter :: kernel_phs = 2

contains

  pure real function dist(x, x0) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)

    res = norm2(x - x0)

  end function dist

  pure real function kernel_ga_eval(x, x0, eps) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    real(r8), intent(in) :: eps

    res = exp(-(dist(x, x0) * eps)**2)

  end function kernel_ga_eval

  pure real function kernel_ga_deriv(x, x0, eps, d, ord) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    real(r8), intent(in) :: eps
    integer, intent(in) :: d
    integer, intent(in) :: ord

    !   exp(-(dist(x, x0) * eps)**2)'
    ! = exp(-(dist(x, x0) * eps)**2) * (-(dist(x, x0) * eps)**2)'
    ! = exp(-(dist(x, x0) * eps)**2) * (-2 * (dist(x, x0) * eps)) * (dist(x, x0) * eps)'
    ! = exp(-(dist(x, x0) * eps)**2) * (-2 * (dist(x, x0) * eps)) * eps * dist(x, x0)'
    ! = exp(-(dist(x, x0) * eps)**2) * (-2 * (dist(x, x0) * eps)) * eps * (x(d) - x0(d)) / dist(x, x0)
    ! = exp(-(dist(x, x0) * eps)**2) * (-2 * eps**2 * (x(d) - x0(d)))

    !   exp(-(dist(x, x0) * eps)**2)''
    ! = (exp(-(dist(x, x0) * eps)**2) * (-2 * eps**2 * (x(d) - x0(d))))'
    ! = exp(-(dist(x, x0) * eps)**2)' * (-2 * eps**2 * (x(d) - x0(d))) + exp(-(dist(x, x0) * eps)**2) * (-2 * eps**2 * (x(d) - x0(d)))'
    ! = exp(-(dist(x, x0) * eps)**2) * (-2 * eps**2 * (x(d) - x0(d)))**2 + exp(-(dist(x, x0) * eps)**2) * (-2 * eps**2)
    ! = exp(-(dist(x, x0) * eps)**2) * ((-2 * eps**2 * (x(d) - x0(d)))**2 + (-2 * eps**2))
    ! = exp(-(dist(x, x0) * eps)**2) * (-2 * eps**2 * ((x(d) - x0(d))**2 + 1))

    select case (ord)
    case (1)
      res = kernel_ga_eval(x, x0, eps) * (-2 * eps**2 * (x(d) - x0(d)))
    case (2)
      res = kernel_ga_eval(x, x0, eps) * (-2 * eps**2 * ((x(d) - x0(d))**2 + 1))
    end select

  end function kernel_ga_deriv

  pure real function kernel_phs_eval(x, x0, pow) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    integer, intent(in) :: pow

    res = dist(x, x0)**pow

  end function kernel_phs_eval

  pure real function kernel_phs_deriv(x, x0, pow, d, ord) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    integer, intent(in) :: pow
    integer, intent(in) :: d
    integer, intent(in) :: ord

    !   dist(x, x0)**pow'
    ! = pow * dist(x, x0)**(pow - 1) * dist(x, x0)'
    ! = pow * dist(x, x0)**(pow - 2) * (x(d) - x0(d))

    !   dist(x, x0)**pow''
    ! = (pow * dist(x, x0)**(pow - 2) * (x(d) - x0(d)))'
    ! = pow * dist(x, x0)**(pow - 2)' * (x(d) - x0(d)) + pow * dist(x, x0)**(pow - 2) * (x(d) - x0(d))'
    ! = pow * (pow - 2) * dist(x, x0)**(pow - 3) * dist(x, x0)' * (x(d) - x0(d)) + pow * dist(x, x0)**(pow - 2)
    ! = pow * (pow - 2) * dist(x, x0)**(pow - 4) * (x(d) - x0(d))**2 + pow * dist(x, x0)**(pow - 2)

    select case (ord)
    case (1)
      res = pow * dist(x, x0)**(pow - 2) * (x(d) - x0(d))
    case (2)
      res = pow * (pow - 2) * dist(x, x0)**(pow - 4) * (x(d) - x0(d))**2 + pow * dist(x, x0)**(pow - 2)
    end select

  end function kernel_phs_deriv

end module kernel_mod
