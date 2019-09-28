program rbf_test

  use fiona
  use distance_mod
  use rbf_mod
  use math_utils_mod

  implicit none

  integer num_dim, num_node, num_ngb
  real(8), allocatable :: x(:,:)
  integer, allocatable :: ngb_idx(:,:)
  type(rbf_type) rbf

  call io_init()

  call io_create_dataset('nodes', file_path='./nodes.nc', mode='input')
  call io_get_dim('nodes', 'dim', size=num_dim)
  call io_get_dim('nodes', 'node', size=num_node)
  call io_get_dim('nodes', 'ngb', size=num_ngb)
  allocate(x(num_dim,num_node))
  allocate(ngb_idx(num_ngb,num_node))
  call io_start_input('nodes')
  call io_input('nodes', 'x', x)
  call io_input('nodes', 'ngb_idx', ngb_idx)
  call io_end_input('nodes')

  call distance_choose(euclidean_distance)
  call rbf%build(x, ngb_idx, kernel_type=kernel_ga, kernel_opt=6.0_r8, poly_order=2)

  deallocate(x)
  deallocate(ngb_idx)

end program rbf_test
