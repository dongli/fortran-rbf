module rbf_mod

  use const_mod, only: r8
  use kernel_mod

  implicit none

  private

  type rbf_interpolant_type
    real(r8), allocatable :: node_wgt_interp(:)
    real(r8), allocatable :: node_wgt_deriv(:,:)
  end type rbf_interpolant_type

contains

  integer function rbf_create(x, kernel, interp) result(ierr)

    real(8), intent(in) :: x(:,:) ! num_dim, num_point
    integer, intent(in) :: kernel ! Kernel function: rbf_kernel_ga, rbf_kernel_phs3, rbf_kernel_phs7
    type(rbf_interpolant_type), intent(inout) :: interp

  end function rbf_create

end module rbf_mod
