module kernel_mod

  ! This module implements several RBF kernels, such as Gaussian.

  use const_mod, only: r8
  use distance_mod, only: d => distance

  implicit none

  private

  public kernel_ga
  public kernel_phs
  public kernel_evaluate_interface
  public kernel_derivative_interface
  public kernel_ga_evaluate
  public kernel_ga_derivative
  public kernel_phs_evaluate
  public kernel_phs_derivative

  integer, parameter :: kernel_ga  = 1
  integer, parameter :: kernel_phs = 2

  interface
    pure real(r8) function kernel_evaluate_interface(x, x0, opt)
      import r8
      real(r8), intent(in) :: x(:)
      real(r8), intent(in) :: x0(:)
      real(r8), intent(in) :: opt
    end function kernel_evaluate_interface

    pure real(r8) function kernel_derivative_interface(x, x0, opt, dim, ord)
      import r8
      real(r8), intent(in) :: x(:)
      real(r8), intent(in) :: x0(:)
      real(r8), intent(in) :: opt
      integer, intent(in) :: dim
      integer, intent(in) :: ord
    end function kernel_derivative_interface
  end interface

contains

  pure real(r8) function kernel_ga_evaluate(x, x0, eps) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    real(r8), intent(in) :: eps

    res = exp(-(d(x, x0) * eps)**2)

  end function kernel_ga_evaluate

  pure real(r8) function kernel_ga_derivative(x, x0, eps, dim, ord) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    real(r8), intent(in) :: eps
    integer, intent(in) :: dim
    integer, intent(in) :: ord

    select case (ord)
    case (1)
      !   exp(-(d(x, x0) * eps)**2)'
      ! = exp(-(d(x, x0) * eps)**2) * (-(d(x, x0) * eps)**2)'
      ! = exp(-(d(x, x0) * eps)**2) * (-2 * (d(x, x0) * eps)) * (d(x, x0) * eps)'
      ! = exp(-(d(x, x0) * eps)**2) * (-2 * (d(x, x0) * eps)) * eps * d(x, x0)'
      ! = exp(-(d(x, x0) * eps)**2) * (-2 * (d(x, x0) * eps)) * eps * (x(dim) - x0(dim)) / d(x, x0)
      ! = exp(-(d(x, x0) * eps)**2) * (-2 * eps**2 * (x(dim) - x0(dim)))
      res = kernel_ga_evaluate(x, x0, eps) * (-2 * eps**2 * (x(dim) - x0(dim)))
    case (2)
      !   exp(-(d(x, x0) * eps)**2)''
      ! = (exp(-(d(x, x0) * eps)**2) * (-2 * eps**2 * (x(dim) - x0(dim))))'
      ! = exp(-(d(x, x0) * eps)**2)' * (-2 * eps**2 * (x(dim) - x0(dim))) + exp(-(d(x, x0) * eps)**2) * (-2 * eps**2 * (x(dim) - x0(dim)))'
      ! = exp(-(d(x, x0) * eps)**2) * (-2 * eps**2 * (x(dim) - x0(dim)))**2 + exp(-(d(x, x0) * eps)**2) * (-2 * eps**2)
      ! = exp(-(d(x, x0) * eps)**2) * ((-2 * eps**2 * (x(dim) - x0(dim)))**2 + (-2 * eps**2))
      ! = exp(-(d(x, x0) * eps)**2) * (-2 * eps**2 * ((x(dim) - x0(dim))**2 + 1))
      res = kernel_ga_evaluate(x, x0, eps) * (-2 * eps**2 * ((x(dim) - x0(dim))**2 + 1))
    end select

  end function kernel_ga_derivative

  pure real(r8) function kernel_phs_evaluate(x, x0, pow) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    real(r8), intent(in) :: pow

    res = d(x, x0)**pow

  end function kernel_phs_evaluate

  pure real(r8) function kernel_phs_derivative(x, x0, pow, dim, ord) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: x0(:)
    real(r8), intent(in) :: pow
    integer, intent(in) :: dim
    integer, intent(in) :: ord

    select case (ord)
    case (1)
      !   d(x, x0)**pow'
      ! = pow * d(x, x0)**(pow - 1) * d(x, x0)'
      ! = pow * d(x, x0)**(pow - 2) * (x(dim) - x0(dim))
      res = pow * d(x, x0)**(pow - 2) * (x(dim) - x0(dim))
    case (2)
      !   d(x, x0)**pow''
      ! = (pow * d(x, x0)**(pow - 2) * (x(dim) - x0(dim)))'
      ! = pow * d(x, x0)**(pow - 2)' * (x(dim) - x0(dim)) + pow * d(x, x0)**(pow - 2) * (x(dim) - x0(dim))'
      ! = pow * (pow - 2) * d(x, x0)**(pow - 3) * d(x, x0)' * (x(dim) - x0(dim)) + pow * d(x, x0)**(pow - 2)
      ! = pow * (pow - 2) * d(x, x0)**(pow - 4) * (x(dim) - x0(dim))**2 + pow * d(x, x0)**(pow - 2)
      res = pow * (pow - 2) * d(x, x0)**(pow - 4) * (x(dim) - x0(dim))**2 + pow * d(x, x0)**(pow - 2)
    end select

  end function kernel_phs_derivative

end module kernel_mod
