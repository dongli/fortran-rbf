module rbf_mod

  ! This module implements radial basis function (RBF).

  use const_mod, only: r8
  use kernel_mod
  use math_utils_mod

  implicit none

  private

  public r8
  public rbf_type
  public kernel_ga
  public kernel_phs

  integer, parameter :: deriv1 = 1
  integer, parameter :: deriv2 = 2

  integer, parameter :: num_deriv_order = 2

  type rbf_type
    procedure(kernel_evaluate_interface  ), pointer, nopass :: phi
    procedure(kernel_derivative_interface), pointer, nopass :: dphi
    integer :: num_ngb    = 0
    integer :: num_dim    = 0
    integer :: num_node   = 0
    integer :: poly_order = -1
    ! Interpolation matrix on each node
    real(r8), allocatable, dimension(:,:,:    ) :: A  ! num_ngb x num_ngb x num_node
    real(r8), allocatable, dimension(:,:,:    ) :: Ai ! Inverse of A
    ! Elemental RBF derivative matrix on each node
    real(r8), allocatable, dimension(:,:,:,:,:) :: B  ! num_ngb x num_ngb x num_deriv_order x num_dim x num_node
    ! Difference matrix on each node
    real(r8), allocatable, dimension(:,:,:,:,:) :: DM ! num_ngb x num_ngb x num_deriv_order x num_dim x num_node
  contains
    procedure :: build => rbf_build
    final :: rbf_final
  end type rbf_type

contains

  subroutine rbf_build(this, x, ngb_idx, kernel_type, kernel_opt, poly_order, y)

    ! The stencil size is implied in the first dimension of ngb_idx.

    class(rbf_type), intent(inout) :: this
    real(r8), intent(in) :: x(:,:)
    integer, intent(in) :: ngb_idx(:,:)
    integer, intent(in) :: kernel_type
    real(r8), intent(in) :: kernel_opt
    integer, intent(in) :: poly_order
    real(r8), intent(in), optional :: y(:)

    integer n
    integer inode, ingb, ingb1, ingb2, dim, ord, iterm, i, j
    integer, allocatable :: poly(:,:)

    select case (kernel_type)
    case (kernel_ga)
      this%phi  => kernel_ga_evaluate
      this%dphi => kernel_ga_derivative
    case (kernel_phs)
      this%phi  => kernel_phs_evaluate
      this%dphi => kernel_phs_derivative
    end select

    if (allocated(this%A )) deallocate(this%A )
    if (allocated(this%Ai)) deallocate(this%Ai)
    if (allocated(this%B )) deallocate(this%B )
    if (allocated(this%DM)) deallocate(this%DM)

    this%num_ngb    = size(ngb_idx, 1)
    this%num_dim    = size(x, 1)
    this%num_node   = size(x, 2)
    this%poly_order = poly_order

    if (poly_order > -1) then
      call create_polynomial(this%num_dim, poly_order, poly)
      n = this%num_ngb + size(poly, 2)
    else
      n = this%num_ngb
    end if

    allocate(this%A (n,n,this%num_node))
    allocate(this%Ai(n,n,this%num_node))
    allocate(this%B (n,n,num_deriv_order,this%num_dim,this%num_node))
    allocate(this%DM(n,n,num_deriv_order,this%num_dim,this%num_node))

    do inode = 1, this%num_node
      do ingb2 = 1, this%num_ngb
        do ingb1 = 1, this%num_ngb
          this%A(ingb1,ingb2,inode) = this%phi(x(:,ngb_idx(ingb1,inode)), &
                                               x(:,ngb_idx(ingb2,inode)), &
                                               kernel_opt)
        end do
      end do
    end do

    if (poly_order > -1) then
      do inode = 1, this%num_node
        this%A(this%num_ngb+1:,:this%num_ngb   ,inode) = 1.0_r8
        this%A(this%num_ngb+1:, this%num_ngb+1:,inode) = 0.0_r8
        do ingb = 1, this%num_ngb
          do iterm = 1, size(poly, 2)
            do j = 1, poly_order
              if (poly(j,iterm) /= 0) then
                i = this%num_ngb+iterm
                this%A(i,ingb,inode) = this%A(i,ingb,inode) * x(poly(j,iterm),ngb_idx(ingb,inode))
              end if
            end do
          end do
        end do
        this%A(:this%num_ngb,this%num_ngb+1:,inode) = transpose(this%A(this%num_ngb+1:,:this%num_ngb,inode))
      end do
      deallocate(poly)
      ! Debug print
      do i = 1, size(this%A, 1)
        do j = 1, size(this%A, 2)
          write(*, '(F5.2)', advance='no') this%A(i,j,10)
          if (j == this%num_ngb) write(*, '(A5)', advance='no') '  |  '
        end do
        write(*, *)
        if (i == this%num_ngb) then
          do j = 1, size(this%A, 2) + 1
            write(*, '(A)', advance='no') '-----'
          end do
          write(*, *)
        end if
      end do
    end if

    do inode = 1, this%num_node
      call inverse_matrix(this%A(:,:,inode), this%Ai(:,:,inode))
    end do
    ! Debug print
    ! do i = 1, size(this%A, 1)
    !   do j = 1, size(this%A, 2)
    !     write(*, '(F8.2)', advance='no') Ai(i,j)
    !     if (j == this%num_ngb) write(*, '(A8)', advance='no') '  |  '
    !   end do
    !   write(*, *)
    !   if (i == this%num_ngb) then
    !     do j = 1, size(this%A, 2) + 1
    !       write(*, '(A)', advance='no') '--------'
    !     end do
    !     write(*, *)
    !   end if
    ! end do

    do inode = 1, this%num_node
      do dim = 1, this%num_dim
        do ord = 1, num_deriv_order
          do ingb2 = 1, this%num_ngb
            do ingb1 = 1, this%num_ngb
              this%B(ingb1,ingb2,ord,dim,inode) = this%dphi(x(:,ngb_idx(ingb1,inode)), &
                                                            x(:,ngb_idx(ingb2,inode)), &
                                                            kernel_opt, dim, ord)
            end do
          end do
          this%DM(:,:,ord,dim,inode) = matmul(transpose(this%B(:,:,ord,dim,inode)), this%Ai(:,:,inode))
        end do
      end do
    end do

  end subroutine rbf_build

  subroutine rbf_final(this)

    type(rbf_type), intent(inout) :: this

    if (allocated(this%A )) deallocate(this%A )
    if (allocated(this%B )) deallocate(this%B )
    if (allocated(this%DM)) deallocate(this%DM)

  end subroutine rbf_final

end module rbf_mod
