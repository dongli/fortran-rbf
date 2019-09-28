program prepare_nodes

  use flogger
  use fiona
  use node_placing_mod
  use kdtree_mod

  implicit none

  real(8), allocatable :: x(:,:)
  integer, allocatable :: ngb_idx(:,:)
  type(kdtree_type) kdtree
  integer i

  call log_notice('Start to run node_placing ...')
  call node_placing([0.0d0,1.0d0,0.0d0,1.0d0], radius, x)

  call kdtree%build(x)

  allocate(ngb_idx(10,size(x, 2)))

  call log_notice('Start to find neighbors of each nodes ...')
  do i = 1, size(x, 2)
    call log_print(to_string(i))
    call kdtree%search(x(:,i), ngb_idx(:,i))
  end do

  call write_nodes(x, ngb_idx)

  deallocate(x)
  deallocate(ngb_idx)

contains

  real(8) function radius(x)

    real(8), intent(in) :: x(2)

    radius = 0.1d0

  end function radius

  subroutine write_nodes(x, ngb_idx)

    real(8), intent(in) :: x(:,:)
    integer, intent(in) :: ngb_idx(:,:)

    call io_init()

    call io_create_dataset('nodes', file_path='nodes.nc', mode='output')
    call io_add_dim('nodes', 'dim' , size=size(x, 1))
    call io_add_dim('nodes', 'node', size=size(x, 2))
    call io_add_dim('nodes', 'ngb' , size=size(ngb_idx, 1))
    call io_add_var('nodes', 'x', 'node coordinates', 'm', ['dim ', 'node'], data_type='real(8)')
    call io_add_var('nodes', 'ngb_idx', 'node neighbor indexing', '', ['ngb ', 'node'], data_type='integer')
    call io_start_output('nodes')
    call io_output('nodes', 'x', x)
    call io_output('nodes', 'ngb_idx', ngb_idx)
    call io_end_output('nodes')

  end subroutine write_nodes

end program prepare_nodes