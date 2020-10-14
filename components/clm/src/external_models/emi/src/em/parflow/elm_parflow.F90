module elm_parflow_module
implicit none
!
integer, dimension(:), allocatable :: pf_l2n       !local to natural
integer, dimension(:), allocatable :: pf_l2g       !local to ghosted
integer, dimension(:), allocatable :: pf_g2n       !ghosted to natural
integer, dimension(:), allocatable :: pf_g2l       !ghosted to local
integer, dimension(:), allocatable :: pf_grid_mask !ghosted grid mask, ghosted 0,local 1
double precision, dimension(:), allocatable :: pf_grid_vol  !ghosted grid volume
double precision, dimension(:), allocatable :: pf_grid_area  !ghosted grid top area
integer :: pf_nxdim            !domain size in x
integer :: pf_nydim
integer :: pf_nzdim
integer :: pf_num_loc_nodes    !number of local nodes
integer :: pf_num_nodes        !number of ghosted nodes
double precision,dimension(:),allocatable :: pf_permx
double precision,dimension(:),allocatable :: pf_permy
double precision,dimension(:),allocatable :: pf_permz
double precision,dimension(:),allocatable :: pf_porosity
double precision,dimension(:),allocatable :: pf_por
double precision,dimension(:),allocatable :: pf_press
double precision,dimension(:),allocatable :: pf_sat
double precision,dimension(:),allocatable :: elm_flux
double precision,dimension(:),allocatable :: pf_alpha
double precision,dimension(:),allocatable :: pf_n
double precision,dimension(:),allocatable :: pf_sr
double precision,dimension(:),allocatable :: z_mult
double precision,dimension(:),allocatable :: pf_grid_dz


contains
!
  subroutine parflow_init(filename)
    character *(*), intent (in) :: filename
    !
    call elmparflowinit(filename, pf_num_nodes, pf_num_loc_nodes, pf_nxdim, pf_nydim, pf_nzdim)
    write(*,*) 'in parflow_init',pf_num_loc_nodes,pf_num_nodes,pf_nxdim,pf_nydim,pf_nzdim
    allocate(pf_l2n(pf_num_loc_nodes))
    allocate(pf_l2g(pf_num_loc_nodes))
    allocate(pf_g2n(pf_num_nodes))
    allocate(pf_g2l(pf_num_nodes))
    allocate(pf_grid_mask(pf_num_nodes))
    allocate(pf_grid_vol(pf_num_nodes))
    allocate(pf_grid_area(pf_num_nodes))
    allocate(pf_grid_dz(pf_num_nodes))
    pf_l2n(:) = 0
    pf_l2g(:) = 0
    pf_g2n(:) = 0
    pf_g2l(:) = 0
    pf_grid_mask(:) = 0
    pf_grid_vol(:) = 0.d0
    pf_grid_dz(:) = 0.d0
    pf_grid_area(:) = 0.d0
    !
    call elmparflowgrid(pf_nxdim, pf_nydim, pf_nzdim, pf_l2n, pf_l2g, pf_g2n, pf_g2l, pf_grid_mask, pf_grid_vol,pf_grid_area,pf_grid_dz)   
!    print*,'vol-',pf_grid_vol
!    print *,'l2n-',pf_l2n(1:50)
!    print *,'g2n-',pf_g2n(1:50)
!    stop 
!
  end subroutine parflow_init

  subroutine parflowsoilprop(num_soil_layers)
    integer :: num_soil_layers
    double precision :: dz
    !
    allocate(pf_permx(pf_num_loc_nodes))
    allocate(pf_permy(pf_num_loc_nodes))
    allocate(pf_permz(pf_num_loc_nodes))
    allocate(pf_porosity(pf_num_loc_nodes))
    allocate(pf_alpha(pf_num_loc_nodes))
    allocate(pf_n(pf_num_loc_nodes))
    allocate(pf_sr(pf_num_loc_nodes))
    allocate(z_mult(pf_num_loc_nodes))
    pf_permx(:) = 0.d0
    pf_permy(:) = 0.d0
    pf_permz(:) = 0.d0
    pf_porosity(:) = 0.d0
    pf_alpha(:) = 0.d0
    pf_n(:) = 0.d0
    pf_sr(:) = 0.d0
    z_mult(:) = 0.d0
    call elmparflowsoilprop(pf_permx,pf_permy,pf_permz,pf_porosity,pf_alpha,pf_n,pf_sr,z_mult,num_soil_layers,0,0,0,0)
    pf_alpha(:) = pf_alpha(:) * 1000.d0 !m to mm
    pf_n(:) = 1.d0/pf_n(:)
    pf_grid_vol(:) = pf_grid_vol(:)*z_mult(:)
    pf_grid_dz(:) = pf_grid_dz(:)*z_mult(:)

  end subroutine parflowsoilprop

end module elm_parflow_module

