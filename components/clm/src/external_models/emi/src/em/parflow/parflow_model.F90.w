module parflow_model_module
!#include "definitions.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscviewer.h"
!#include "finclude/petscvec.h90"
  use petscsys
  use Mapping_module

 
  use petscsys
  use petscvec
  use elm_parflow_module  
  use input_read_module
  implicit none


  ! Note:
  !
  ! CLM has the following:
  !   (i) 3D subsurface grid (CLM_SUB);
  !   (ii) 2D surface grid (CLM_SRF).
  ! CLM decomposes the 3D subsurface grid across processors in a 2D (i.e.
  ! cells in Z are not split across processors). Thus, the surface cells of
  ! 3D subsurface grid are on the same processors as the 2D surface grid.
  !
  ! PARFLOW has the following:
  !   (i) 3D subsurface grid (PF_SUB);
  !   (ii) surface control volumes of 3D subsurface grid (PF_2DSUB);
  !   (iii) 2D surface grid (PF_SRF).

  ! map level constants
  PetscInt, parameter, public :: CLM_SUB_TO_PF_SUB           = 1 ! 3D --> 3D
  PetscInt, parameter, public :: CLM_SUB_TO_PF_EXTENDED_SUB  = 2 ! 3D --> extended 3D
  PetscInt, parameter, public :: CLM_SRF_TO_PF_2DSUB         = 3 ! 2D --> SURF of 3D grid
  PetscInt, parameter, public :: CLM_SRF_TO_PF_SRF           = 4 ! 2D --> 2D SURF grid
  PetscInt, parameter, public :: PF_SUB_TO_CLM_SUB           = 5 ! 3D --> 3D
  PetscInt, parameter, public :: PF_SRF_TO_CLM_SRF           = 6 ! 2D SURF grid --> 2D

  ! mesh ids
  PetscInt, parameter, public :: CLM_SUB_MESH   = 1
  PetscInt, parameter, public :: CLM_SRF_MESH   = 2
  PetscInt, parameter, public :: PF_SUB_MESH    = 3
  PetscInt, parameter, public :: PF_2DSUB_MESH  = 4
  PetscInt, parameter, public :: PF_SRF_MESH    = 5

  type, public :: inside_each_overlapped_cell
     PetscInt           :: id
     PetscInt           :: ocell_count
     PetscInt,  pointer :: ocell_id(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_overlapped_cell

  type, public :: parflow_model_type
    
    character(len=MAXSTRINGLENGTH)                         :: map_filename
    type(inside_each_overlapped_cell), pointer :: pf_cells(:)
    type(inside_each_overlapped_cell), pointer :: clm_cells(:)
    type(mapping_type),                pointer :: map_clm_sub_to_pf_sub
    type(mapping_type),                pointer :: map_clm_sub_to_pf_extended_sub
    type(mapping_type),                pointer :: map_clm_srf_to_pf_2dsub
    type(mapping_type),                pointer :: map_clm_srf_to_pf_srf
    type(mapping_type),                pointer :: map_pf_sub_to_clm_sub
    type(mapping_type),                pointer :: map_pf_srf_to_clm_srf
  end type parflow_model_type

  public::parflowModelCreate,               &
       parflowModelInitMapping,             &
       parflowModelSetSoilProp,             &
       parflowModelGetSoilProp,             &
       parflowModelSetICs,                  &
       parflowModelGetSaturation,           &   
       parflowModelDestroy

  private :: &
       parflowModelSetupMappingFiles


contains

! ************************************************************************** !

  function parflowModelCreate(mpicomm, mycommsize,myrank,parflow_file,mapfile)
  ! 
  ! Allocates and initializes the parflowModel object.
  ! 


    implicit none

#include "petsc/finclude/petscsys.h"

    integer, intent(in) :: mpicomm
    integer, intent(in)  :: mycommsize
    integer, intent(in)  :: myrank
    character(len=*), intent(in) :: mapfile 
    character(len=*), intent(in) :: parflow_file

    type(parflow_model_type), pointer :: parflowModelCreate

    type(parflow_model_type),      pointer :: model

    allocate(model)

    call parflow_init(parflow_file)
    !file containing mapfile names
    if (len(trim(mapfile)) > 1) then
      model%map_filename = trim(mapfile)
    else
      print *, 'The external driver must provide the ' // &
           'parflow input file prefix.'
    end if

    PETSC_COMM_SELF = MPI_COMM_SELF
    PETSC_COMM_WORLD = MPI_COMM_WORLD

    call parflowModelSetupMappingFiles(model,mapfile,mycommsize,myrank)
    !assign model to output
    parflowModelCreate => model

  end function parflowModelCreate

! ************************************************************************** !

  subroutine parflowModelSetupMappingFiles(model,mapfile,mycommsize,myrank)
  ! 
  ! parflowModelSetupMappingFiles

    use Mapping_module

    implicit none

#include "petsc/finclude/petscsys.h"

    type(parflow_model_type), pointer, intent(inout) :: model

    type(mapping_type), pointer        :: map

    PetscBool :: clm2pf_flux_file
    PetscBool :: clm2pf_soil_file
    PetscBool :: clm2pf_gflux_file
    PetscBool :: clm2pf_rflux_file
    PetscBool :: pf2clm_flux_file
    PetscBool :: pf2clm_surf_file
    character(len=MAXSTRINGLENGTH) :: string
    character(len=MAXSTRINGLENGTH) :: buf
    character(len=MAXSTRINGLENGTH) :: path
    character(len=MAXSTRINGLENGTH) :: filename
    character(len=MAXSTRINGLENGTH) :: mapfile
    character(len=MAXWORDLENGTH) :: word
    integer                      :: fid
    integer                      :: ierr
    integer                      :: mycommsize
    integer                      :: myrank

    nullify(model%pf_cells)
    nullify(model%clm_cells)
    nullify(model%map_clm_sub_to_pf_sub)
    !nullify(model%map_clm_sub_to_pf_extended_sub)
    !nullify(model%map_clm_srf_to_pf_2dsub)
    !nullify(model%map_clm_srf_to_pf_srf)
    nullify(model%map_pf_sub_to_clm_sub)
    !nullify(model%map_pf_srf_to_clm_srf)

    model%map_clm_sub_to_pf_sub          => MappingCreate()
    !model%map_clm_sub_to_pf_extended_sub => MappingCreate()
    !model%map_clm_srf_to_pf_2dsub        => MappingCreate()
    !model%map_clm_srf_to_pf_srf          => MappingCreate()
    model%map_pf_sub_to_clm_sub          => MappingCreate()
    !model%map_pf_srf_to_clm_srf          => MappingCreate()


    ! Read names of mapping file
    clm2pf_flux_file   = PETSC_FALSE
    clm2pf_soil_file   = PETSC_FALSE
    !clm2pf_gflux_file  = PETSC_FALSE
    !clm2pf_rflux_file  = PETSC_FALSE
    pf2clm_flux_file   = PETSC_FALSE
    !pf2clm_surf_file   = PETSC_FALSE
    
!    string = "MAPPING_FILES"
    fid = 25
    call InputCreate(fid,path,mapfile)
!    call InputFindStringInFile(filename,string,ierr)
!    read(fid,'(a512)',iostat=ierr) buf

    do
      if(clm2pf_flux_file .and. pf2clm_flux_file) exit
      call InputReadString(fid, buf,ierr)
      if (ierr /= 0) exit

      call InputReadWord(buf,word,.true.,ierr)
      call StringToUpper(word)
      
      select case(trim(word))
        case('CLM2PF_FLUX_FILE')
          map => model%map_clm_sub_to_pf_sub
          call InputReadNChars(buf, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE,ierr)
          map%filename     = trim(map%filename)//CHAR(0)
          clm2pf_flux_file = PETSC_TRUE
#if 0
        case('CLM2PF_SOIL_FILE')
          map => model%map_clm_sub_to_pf_extended_sub
          call InputReadNChars(buf, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE,ierr)
          map%filename     = trim(map%filename)//CHAR(0)
          clm2pf_soil_file = PETSC_TRUE
#endif

        case('PF2CLM_FLUX_FILE')
          map => model%map_pf_sub_to_clm_sub
          call InputReadNChars(buf, map%filename, &
                               MAXSTRINGLENGTH, PETSC_TRUE,ierr)
          map%filename     = trim(map%filename)//CHAR(0)
          pf2clm_flux_file = PETSC_TRUE

        case default
!          model%option%io_buffer='Keyword ' // trim(word) // &
!            ' in input file not recognized'
!          call printErrMsg(model%option)
      end select

      ! Read mapping file
!#if 0
!      if (index(map%filename, '.h5') > 0) then
!        call MappingReadHDF5(map, map%filename, mycommsize, myrank)
!      else
        call MappingReadTxtFile(map, map%filename, mycommsize, myrank)
!      endif
!#endif
        
    enddo
    close(fid)
#if 0
    if ((.not. clm2pf_soil_file) .or. (.not. clm2pf_flux_file) .or. &
        (.not. pf2clm_flux_file) ) then
      print *, 'One of the mapping files not found'
    endif
#endif    

  end subroutine parflowModelSetupMappingFiles


! ************************************************************************** !

subroutine parflowModelSetICs(parflow_model)

  ! 
  ! Set initial conditions
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/22/2010
  ! 

    use Mapping_module
    use clm_parflow_interface_data

    implicit none

    type(parflow_model_type), pointer        :: parflow_model

    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: xx_loc_p(:)

    PetscScalar, pointer :: press_pf_loc(:) ! Pressure [Pa]
  ! values passed from parflow
    double precision :: pressure(pf_num_nodes)     ! pressure head, from parflow on grid w/ ghost nodes for current proc


    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_sub, &
                                    clm_pf_idata%press_clm, &
                                    clm_pf_idata%press_pf)

    call VecGetArrayF90(clm_pf_idata%press_pf, press_pf_loc, ierr)

    call VecRestoreArrayF90(clm_pf_idata%press_pf, press_pf_loc, ierr)

end subroutine parflowModelSetICs

! ************************************************************************** !

  subroutine parflowModelSetSoilProp(parflow_model)
  ! 
  ! converts hydraulic properties from clm units
  ! into parflow units.
  ! #ifdef clm_parflow
  ! 
  ! author: gautam bisht
  ! date: 10/22/2010
  ! 

    use clm_parflow_interface_data
    use mapping_module

    implicit none

    type(parflow_model_type), pointer        :: parflow_model

    PetscErrorCode     :: ierr
    PetscInt           :: ghosted_id
    PetscReal          :: den, vis, grav
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: bc_lambda, bc_alpha
    Vec                :: porosity_loc, perm_xx_loc, perm_yy_loc, perm_zz_loc

    PetscScalar, pointer :: hksat_x_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm h2o /s)
    PetscScalar, pointer :: hksat_y_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm h2o /s)
    PetscScalar, pointer :: hksat_z_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm h2o /s)
    PetscScalar, pointer :: watsat_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw_pf_loc(:)     ! clapp and hornberger "b"

    den = 998.2d0       ! [kg/m^3]  @ 20 degc
    vis = 0.001002d0    ! [n s/m^2] @ 20 degc
    grav = 9.81d0       ! [m/s^2]

#if 0
    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_extended_sub, &
                                    clm_pf_idata%hksat_x_clm, &
                                    clm_pf_idata%hksat_x_pf)

    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_extended_sub, &
                                    clm_pf_idata%hksat_y_clm, &
                                    clm_pf_idata%hksat_y_pf)

    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_extended_sub, &
                                    clm_pf_idata%hksat_z_clm, &
                                    clm_pf_idata%hksat_z_pf)

    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_extended_sub, &
                                    clm_pf_idata%sucsat_clm, &
                                    clm_pf_idata%sucsat_pf)

    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_extended_sub, &
                                    clm_pf_idata%bsw_clm, &
                                    clm_pf_idata%bsw_pf)

    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_extended_sub, &
                                    clm_pf_idata%watsat_clm, &
                                    clm_pf_idata%watsat_pf)
    call VecGetArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_pf, hksat_x_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_pf, hksat_y_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_pf, hksat_z_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_pf,  sucsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_pf,  watsat_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_pf,     bsw_pf_loc,     ierr)

#endif

  end subroutine parflowModelSetSoilProp

! ************************************************************************** !

  subroutine parflowmodelgetsoilprop(parflow_model)
  !
  ! converts hydraulic properties from parflow units for clm units
  !
  ! author: gautam bisht
  ! date: 10/27/2014
  !

    use clm_parflow_interface_data
    use mapping_module

    implicit none

    type(parflow_model_type), pointer        :: parflow_model

    PetscErrorCode     :: ierr
    PetscInt           :: ghosted_id, local_id
    PetscInt , parameter :: liq_iphase = 1
    PetscReal          :: den, vis, grav, sr
    PetscReal, pointer :: porosity_loc_p(:), vol_ovlap_arr(:)
    PetscReal, pointer :: perm_xx_loc_p(:), perm_yy_loc_p(:), perm_zz_loc_p(:)
    PetscReal          :: bc_lambda, bc_alpha
    Vec                :: porosity_loc, perm_xx_loc, perm_yy_loc, perm_zz_loc

    PetscScalar, pointer :: hksat_x2_pf_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y2_pf_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z2_pf_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat2_pf_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat2_pf_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw2_pf_loc(:)     ! Clapp and Hornberger "b"
    PetscScalar, pointer :: thetares2_pf_loc(:)! residual soil mosture = sat_res * por

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/S^2]
#if 0

    call VecDuplicate(field%work_loc, porosity_loc, ierr)
    call VecDuplicate(field%work_loc, perm_xx_loc, ierr)
    call VecDuplicate(field%work_loc, perm_yy_loc, ierr)
    call VecDuplicate(field%work_loc, perm_zz_loc, ierr)


    call VecGetArrayF90(clm_pf_idata%hksat_x2_pf, hksat_x2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y2_pf, hksat_y2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z2_pf, hksat_z2_pf_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat2_pf,  sucsat2_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%watsat2_pf,  watsat2_pf_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%bsw2_pf,     bsw2_pf_loc,     ierr)
    call VecGetArrayF90(clm_pf_idata%thetares2_pf,thetares2_pf_loc,ierr)

    call VecGetArrayF90(porosity_loc, porosity_loc_p, ierr)
    call VecGetArrayF90(perm_xx_loc,  perm_xx_loc_p,  ierr)
    call VecGetArrayF90(perm_yy_loc,  perm_yy_loc_p,  ierr)
    call VecGetArrayF90(perm_zz_loc,  perm_zz_loc_p,  ierr)

    hksat_x2_pf_loc(:) = 0.d0
    hksat_y2_pf_loc(:) = 0.d0
    hksat_z2_pf_loc(:) = 0.d0
    sucsat2_pf_loc(:) = 0.d0
    watsat2_pf_loc(:) = 0.d0
    bsw2_pf_loc(:) = 0.d0
    thetares2_pf_loc(:) = 0.d0

    do ghosted_id = 1, pf_num_nodes

      if (patch%sat_func_id(ghosted_id) < 1) cycle

      local_id = pf_g2l(ghosted_id)


      ! bc_alpha [1/Pa]; while sucsat [mm of H20]
      ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
      sucsat2_pf_loc(local_id) = 1.d0/(bc_alpha*grav)

      ! bc_lambda = 1/bsw
      bsw2_pf_loc(local_id) = 1.d0/bc_lambda

      ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
      ! [m^2]          [mm/sec]
      hksat_x2_pf_loc(local_id) = perm_xx_loc_p(ghosted_id)/vis*(den*grav)*1000.d0
      hksat_y2_pf_loc(local_id) = perm_yy_loc_p(ghosted_id)/vis*(den*grav)*1000.d0
      hksat_z2_pf_loc(local_id) = perm_zz_loc_p(ghosted_id)/vis*(den*grav)*1000.d0

      watsat2_pf_loc(local_id) = porosity_loc_p(ghosted_id)


      thetares2_pf_loc(local_id) = porosity_loc_p(ghosted_id)*Sr

   enddo

    call VecRestoreArrayF90(clm_pf_idata%hksat_x2_pf, hksat_x2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y2_pf, hksat_y2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z2_pf, hksat_z2_pf_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat2_pf,  sucsat2_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat2_pf,  watsat2_pf_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw2_pf,     bsw2_pf_loc,     ierr)
    call VecRestoreArrayF90(clm_pf_idata%thetares2_pf,thetares2_pf_loc,ierr)

    call VecRestoreArrayF90(porosity_loc, porosity_loc_p, ierr)
    call VecRestoreArrayF90(perm_xx_loc,  perm_xx_loc_p,  ierr)
    call VecRestoreArrayF90(perm_yy_loc,  perm_yy_loc_p,  ierr)
    call VecRestoreArrayF90(perm_zz_loc,  perm_zz_loc_p,  ierr)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%hksat_x2_pf, &
                                    clm_pf_idata%hksat_x2_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%hksat_y2_pf, &
                                    clm_pf_idata%hksat_y2_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%hksat_z2_pf, &
                                    clm_pf_idata%hksat_z2_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%sucsat2_pf, &
                                    clm_pf_idata%sucsat2_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%bsw2_pf, &
                                    clm_pf_idata%bsw2_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%watsat2_pf, &
                                    clm_pf_idata%watsat2_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%thetares2_pf, &
                                    clm_pf_idata%thetares2_clm)

    call VecDestroy(porosity_loc, ierr)
    call VecDestroy(perm_xx_loc, ierr)
    call VecDestroy(perm_yy_loc, ierr)
    call VecDestroy(perm_zz_loc, ierr)
#endif
  end subroutine parflowModelGetSoilProp

! ************************************************************************** !

  subroutine parflowModelInitMapping(parflow_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id,mycomm,myrank)
  ! 
  ! #endif
  ! Initialize mapping between the two model grid
  ! (CLM and PFLTORAN)
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/24/2011
  ! 

    use Mapping_module
    use elm_parflow_module

    implicit none

    type(parflow_model_type), intent(inout), pointer :: parflow_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename
    integer                                           :: mycomm
    integer                                           :: myrank

    select case (map_id)
      case (CLM_SUB_TO_PF_SUB, CLM_SUB_TO_PF_EXTENDED_SUB, PF_SUB_TO_CLM_SUB)

        call parflowModelInitMappingSub2Sub(parflow_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id,mycomm,myrank)
      case (CLM_SRF_TO_PF_2DSUB)


!        call parflowModelInitMapSrfTo2DSub(parflow_model,  &
!                                            grid_clm_cell_ids_nindex, &
!                                            grid_clm_npts_local, &
!                                            map_id)
      case (CLM_SRF_TO_PF_SRF, PF_SRF_TO_CLM_SRF)
!        call parflowModelInitMapSrfToSrf(parflow_model,  &
!                                            grid_clm_cell_ids_nindex, &
!                                            grid_clm_npts_local, &
!                                            map_id)
      case default
!        parflow_model%option%io_buffer = 'Invalid map_id argument to ' // &
!          'parflowModelInitMapping'
!        call printErrMsg(parflow_model%option)
    end select

  end subroutine parflowModelInitMapping

! ************************************************************************** !

  subroutine parflowModelInitMappingSub2Sub(parflow_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id,mycomm,myrank)
  ! 
  ! Initialize mapping between 3D subsurface
  ! CLM grid and 3D subsurface PARFLOW grid.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/09/2013
  ! 

    use Mapping_module

    implicit none

    type(parflow_model_type), intent(inout), pointer :: parflow_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    character(len=MAXSTRINGLENGTH)                    :: filename

    ! local
    PetscInt                           :: local_id, ghosted_id
    PetscInt                           :: grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)

    type(mapping_type), pointer        :: map
    integer                            :: mycomm
    integer                            :: myrank

    !
    ! Mapping to/from entire PARFLOW 3D subsurface domain
    !

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM_SUB_TO_PF_SUB)
        map => parflow_model%map_clm_sub_to_pf_sub
        source_mesh_id = CLM_SUB_MESH
        dest_mesh_id = PF_SUB_MESH
#if 0
      case(CLM_SUB_TO_PF_EXTENDED_SUB)
        map => parflow_model%map_clm_sub_to_pf_extended_sub
        source_mesh_id = CLM_SUB_MESH
        dest_mesh_id = PF_SUB_MESH
#endif
      case(PF_SUB_TO_CLM_SUB)
        map => parflow_model%map_pf_sub_to_clm_sub
        source_mesh_id = PF_SUB_MESH
      case default
        print *, 'Invalid map_id argument to parflowModelInitMapping'
    end select


    grid_clm_npts_ghost = 0

    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = 1 ! LOCAL
    enddo

    ! Find cell IDs for PARFLOW grid
    grid_pf_npts_local = pf_num_loc_nodes
    grid_pf_npts_ghost = pf_num_nodes - pf_num_loc_nodes

    allocate(grid_pf_local_nindex(pf_num_nodes))
    grid_pf_local_nindex = pf_grid_mask !0-ghost, 1-local

    select case(source_mesh_id)
      case(CLM_SUB_MESH)

        allocate(grid_pf_cell_ids_nindex(pf_num_nodes))
        do ghosted_id = 1, pf_num_nodes
          grid_pf_cell_ids_nindex(ghosted_id) = pf_g2n(ghosted_id)-1
        enddo
     
        call MappingSetSourceMeshCellIds(map, grid_clm_npts_local, &
                                         grid_clm_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, grid_pf_npts_local, &
                                              grid_pf_npts_ghost, &
                                              grid_pf_cell_ids_nindex, &
                                              grid_pf_local_nindex)
      case(PF_SUB_MESH)

        allocate(grid_pf_cell_ids_nindex(pf_num_loc_nodes))
        do local_id = 1, pf_num_loc_nodes
          grid_pf_cell_ids_nindex(local_id) = pf_l2n(local_id)-1
        enddo

        call MappingSetSourceMeshCellIds(map, grid_pf_npts_local, &
                                        grid_pf_cell_ids_nindex)
        call MappingSetDestinationMeshCellIds(map, grid_clm_npts_local, &
                                              grid_clm_npts_ghost, &
                                              grid_clm_cell_ids_nindex, &
                                              grid_clm_local_nindex)
      case default
!        option%io_buffer = 'Invalid argument source_mesh_id passed to parflowModelInitMapping'
!        call printErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_local_nindex)

    call MappingDecompose(map, mycomm)
    call MappingFindDistinctSourceMeshCellIds(map)
    call MappingCreateWeightMatrix(map, myrank)
    call MappingCreateScatterOfSourceMesh(map, mycomm)

  end subroutine parflowModelInitMappingSub2Sub



! ************************************************************************** !

  subroutine parflowModelUpdateSourceSink(parflow_model,pf_flux)
  ! 
  ! Update the source/sink term
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! 

    use clm_parflow_interface_data
    use Mapping_module

    implicit none

    type(parflow_model_type), pointer        :: parflow_model

    PetscScalar, pointer                      :: qflx_pf_loc(:)
    PetscBool                                 :: found
    PetscInt                                  :: iconn
    PetscInt                                  :: local_id
    PetscInt                                  :: ghosted_id
    PetscErrorCode                            :: ierr
    PetscInt                                  :: press_dof
    double precision                          :: pf_flux(pf_num_nodes)

    call MappingSourceToDestination(parflow_model%map_clm_sub_to_pf_sub, &
                                    clm_pf_idata%qflx_clm, &
                                    clm_pf_idata%qflx_pf)


    ! Update the source/sink term to parflow
    call VecGetArrayF90(clm_pf_idata%qflx_pf,qflx_pf_loc,ierr)
    do local_id=1, pf_num_loc_nodes
      ghosted_id=pf_l2g(local_id)
      pf_flux(ghosted_id)=qflx_pf_loc(local_id)
    enddo
    call VecRestoreArrayF90(clm_pf_idata%qflx_pf,qflx_pf_loc,ierr)

  end subroutine parflowModelUpdateSourceSink


! ************************************************************************** !

  subroutine parflowModelGetSaturation(parflow_model,sat,vol,porosity)
  ! 
  ! Extract soil saturation values simulated by
  ! PARFLOW in a PETSc vector.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/22/2011
  ! 

    use clm_parflow_interface_data
    use Mapping_module

    implicit none

    type(parflow_model_type), pointer        :: parflow_model
    PetscErrorCode     :: ierr
    PetscInt           :: local_id, ghosted_id
    PetscReal, pointer :: sat_pf_p(:)
    PetscReal, pointer :: mass_pf_p(:)
    integer            :: nx
    integer            :: ny
    double precision   :: sat(pf_num_nodes)
    double precision   :: vol(pf_num_nodes)
    double precision   :: porosity(pf_num_nodes)
    
    ! Save the saturation values
    call VecGetArrayF90(clm_pf_idata%sat_pf, sat_pf_p, ierr)
    call VecGetArrayF90(clm_pf_idata%mass_pf, mass_pf_p, ierr)
    do local_id=1, pf_num_loc_nodes
      ghosted_id=pf_l2g(local_id)
      sat_pf_p(local_id)=sat(ghosted_id)
      mass_pf_p(local_id)= sat(ghosted_id)*vol(ghosted_id)*porosity(ghosted_id)*1000.d0
    enddo
    call VecRestoreArrayF90(clm_pf_idata%sat_pf, sat_pf_p, ierr)
    call VecRestoreArrayF90(clm_pf_idata%mass_pf, mass_pf_p, ierr)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%sat_pf, &
                                    clm_pf_idata%sat_clm)

    call MappingSourceToDestination(parflow_model%map_pf_sub_to_clm_sub, &
                                    clm_pf_idata%mass_pf, &
                                    clm_pf_idata%mass_clm)


  end subroutine parflowModelGetSaturation

! ************************************************************************** !

  subroutine parflowModelDestroy(model)
  ! 
  ! Deallocates the parflowModel object
  ! 
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! 

    use Mapping_module, only : MappingDestroy

    implicit none

    type(parflow_model_type), pointer :: model

    if (associated(model%map_pf_sub_to_clm_sub)) then
      call MappingDestroy(model%map_pf_sub_to_clm_sub)
      nullify(model%map_pf_sub_to_clm_sub)
    endif

    deallocate(model)

  end subroutine parflowModelDestroy
  
end module parflow_model_module

