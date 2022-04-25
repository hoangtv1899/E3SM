module ExternalModelPARFLOWMod
#ifdef USE_PETSC_LIB

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscviewer.h"

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod                  , only : emi_data_list, emi_data
  use clm_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use MultiPhysicsProbVSFM         , only : mpp_vsfm_type
  use decompMod                    , only : bounds_type
  use parflow_model_module
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use clm_parflow_interface_data
  use elm_parflow_module
  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_parflow_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_gridcell_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_state_wtd
     integer :: index_l2e_init_state_soilp

     integer :: index_l2e_init_h2osoi_liq
     integer :: index_l2e_init_h2osoi_ice

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice
     integer :: index_e2l_init_state_h2osoi_vol
     integer :: index_e2l_init_state_wtd
     integer :: index_e2l_init_parameter_watsatc
     integer :: index_e2l_init_parameter_hksatc
     integer :: index_e2l_init_parameter_bswc
     integer :: index_e2l_init_parameter_sucsatc

     integer :: index_e2l_init_flux_mflx_snowlyr_col
     integer :: index_l2e_init_flux_mflx_snowlyr_col

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_state_tsoil
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_wtd
     integer :: index_e2l_state_soilp

     integer :: index_l2e_flux_infil
     integer :: index_l2e_flux_et
     integer :: index_l2e_flux_dew
     integer :: index_l2e_flux_snow_sub
     integer :: index_l2e_flux_snowlyr
     integer :: index_l2e_flux_drainage

     integer :: index_e2l_flux_qrecharge
     integer :: index_e2l_flux_drain_perched
     integer :: index_e2l_flux_drain
     integer :: index_e2l_flux_qrgwl
     integer :: index_e2l_flux_rsub_sat

     integer :: index_l2e_filter_hydrologyc
     integer :: index_l2e_filter_num_hydrologyc

     integer :: index_l2e_column_zi
     integer :: index_l2e_col_active
     integer :: index_l2e_col_gridcell
     integer :: index_l2e_col_dz

     type(parflow_model_type), pointer :: parflow_m

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_PARFLOW_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_PARFLOW_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_PARFLOW_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_PARFLOW_Populate_E2L_List
     procedure, public :: PreInit                 => EM_PARFLOW_PreInit
     procedure, public :: Init                    => EM_PARFLOW_Init
     procedure, public :: Solve                   => EM_PARFLOW_Solve

  end type em_parflow_type

contains

  !------------------------------------------------------------------------
  subroutine EM_PARFLOW_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PARFLOW from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)             :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_wtd              = index

    id                                         = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_soilp            = index

    id                                         = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_flux_mflx_snowlyr_col  = index

    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index     = index

    id                                         = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_gridcell_index     = index

    id                                         = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi                 = index

    id                                         = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz                 = index

    id                                         = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z                  = index

    id                                         = L2E_COLUMN_AREA
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_area               = index

    id                                         = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type          = index

    id                                         = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint     = index

    id                                         = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint    = index

    id                                         = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsatc      = index

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    id                                         = L2E_PARAMETER_SUCSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_sucsatc      = index

    id                                         = L2E_PARAMETER_EFFPOROSITYC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_effporosityc = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq            = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice            = index

    deallocate(em_stages)

  end subroutine EM_PARFLOW_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PARFLOW_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PARFLOW from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)             :: this
    class(emi_data_list), intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_H2OSOI_VOL_NLEVGRND
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_vol      = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index

    id                                         = E2L_PARAMETER_WATSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_watsatc      = index

    id                                         = E2L_PARAMETER_HKSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_hksatc       = index

    id                                         = E2L_PARAMETER_BSWC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_bswc         = index

    id                                         = E2L_PARAMETER_SUCSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_sucsatc      = index

    deallocate(em_stages)

  end subroutine EM_PARFLOW_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PARFLOW_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PARFLOW from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)             :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    !
    ! !LOCAL VARIABLES:
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PARFLOW_SOIL_HYDRO_STAGE

    id                                   = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoil           = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_INFIL_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infil            = index

    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    id                                   = L2E_FLUX_DEW_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew              = index

    id                                   = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow_sub         = index

    id                                   = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snowlyr          = index

    id                                   = L2E_FLUX_DRAINAGE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drainage         = index

    id                                   = L2E_FILTER_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_hydrologyc     = index

    id                                   = L2E_FILTER_NUM_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_hydrologyc = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    id                                   = L2E_COLUMN_ACTIVE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_active            = index

    id                                   = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_gridcell          = index

    id                                    = L2E_COLUMN_DZ
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_dz                 = index

    deallocate(em_stages)

  end subroutine EM_PARFLOW_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_PARFLOW_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PARFLOW from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)             :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PARFLOW_SOIL_HYDRO_STAGE

    id                              = E2L_STATE_H2OSOI_LIQ
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_liq = index

    id                              = E2L_STATE_H2OSOI_ICE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_ice = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_wtd        = index

    id                              = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_soilp      = index

    id                              = E2L_FLUX_AQUIFER_RECHARGE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrecharge   = index

    id                              = E2L_FLUX_DRAIN_PERCHED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_drain_perched   = index

    id                              = E2L_FLUX_DRAIN
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_drain       = index

    id                              = E2L_FLUX_QRGWL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrgwl   = index

    id                              = E2L_FLUX_RSUB_SAT
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_rsub_sat   = index

    deallocate(em_stages)

  end subroutine EM_PARFLOW_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_PARFLOW_PreInit(this, prefix,mapfile, restart_stamp)
    !
    use spmdMod     , only : mpicom, npes, iam
    !
    implicit none
    !
    class(em_parflow_type)      :: this
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: mapfile
    character(len=*), intent(in) :: restart_stamp
    
    this%parflow_m => parflowModelCreate(mpicom, npes,iam,prefix,mapfile)
!    call pflotranModelSetupRestart(this%pflotran_m, restart_stamp)

  end subroutine EM_PARFLOW_PreInit

  !------------------------------------------------------------------------
  subroutine EM_PARFLOW_Init(this, l2e_init_list, e2l_init_list, iam,bounds_clump)

    !
    !use spmdMod         , only : mpicom, masterproc, npes, iam
    use elm_parflow_module
    !
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)              :: this

    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    integer              , intent(in)    :: iam
    !
    ! LOCAL VARAIBLES:
    integer :: clm_npts
    integer :: clm_surf_npts

    ! Initialize parflow model, num_nodes,num_loc_nodes,
    ! pf_l2n,pf_gln,pf_grid_mask,pf_grid_vol
    ! then Create mapping model
    !parflow_m => parflowModelCreate(mpicom,npes,iam,prefix,mapfile)

    ! Initialize size and vector for clm and parflow
    call CreateCLMPARFLOWInterfaceDate(this, bounds_clump, clm_npts, clm_surf_npts)

    ! Create CLM-PARFLOW mapping files
    call CreateCLMPARFLOWMaps(this, bounds_clump, clm_npts, clm_surf_npts)

    ! Initialize PARFLOW states

    !call parflowModelStepperRunInit(this%parflow_m)

    ! Get top surface area
    !call parflowModelGetTopFaceArea(this%parflow_m)

    ! Get PARFLOW states
    !call parflowModelGetUpdatedData(this%parflow_m)
    ! Get parflow soil properties
    call parflowsoilprop(this%parflow_m%map_clm_sub_to_pf_sub%clm_nlevsoi)
    call parflowModelGetTopFaceArea(this%parflow_m)

    ! Save the data need by ELM
    call extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)

  end subroutine EM_PARFLOW_Init

  !-----------------------------------------------------------------------
  subroutine CreateCLMPARFLOWInterfaceDate(this, bounds, clm_npts, clm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Allocates memory for CLM-PARFLOW interface
    !
    ! !USES:
    use decompMod                     , only : bounds_type
    use clm_varpar                    , only : nlevsoi, nlevgrnd
    use shr_kind_mod                  , only: r8 => shr_kind_r8
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(inout) :: clm_npts
    integer                 , intent(inout) :: clm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer                                      :: nlevmapped
    character(len= 128)                          :: subname = 'CreateCLMPARFLOWInterfaceDate' ! subroutine name

    ! Initialize PETSc vector for data transfer between CLM and PARFLOW, defined in
    ! parflow_dir/elm
    call CLMPARFLOWIDataInit()


    ! Compute number of cells in CLM domain.
    ! Assumption-1: One column per CLM grid cell.

    ! Check if the number of CLM vertical soil layers defined in the mapping
    ! file read by PARFLOW matches either nlevsoi or nlevgrnd
    clm_pf_idata%nzclm_mapped = this%parflow_m%map_clm_sub_to_pf_sub%clm_nlevsoi
    nlevmapped                = clm_pf_idata%nzclm_mapped
    if ( (nlevmapped /= nlevsoi) .and. (nlevmapped /= nlevgrnd) ) then
       call endrun(trim(subname)//' ERROR: Number of layers PARFLOW thinks CLM should '// &
            'have do not match either nlevsoi or nlevgrnd. Abortting' )
    end if

    clm_npts = (bounds%endg - bounds%begg + 1)*nlevmapped
    clm_surf_npts = (bounds%endg - bounds%begg + 1)

    ! CLM: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_sub = clm_npts
    clm_pf_idata%ngclm_sub = clm_npts

    ! CLM: Surface of subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_2dsub = (bounds%endg - bounds%begg + 1)
    clm_pf_idata%ngclm_2dsub = (bounds%endg - bounds%begg + 1)
    ! For CLM: Same as surface of subsurface domain
    clm_pf_idata%nlclm_srf = clm_surf_npts
    clm_pf_idata%ngclm_srf = clm_surf_npts

    ! PARFLOW: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlpf_sub = pf_num_loc_nodes
    clm_pf_idata%ngpf_sub = pf_num_nodes

    ! PARFLOW: Surface of subsurface domain (local and ghosted cells)
    !clm_pf_idata%nlpf_2dsub = 0
    !clm_pf_idata%ngpf_2dsub = 0

    ! PARFLOW: Surface domain (local and ghosted cells)
    !clm_pf_idata%nlpf_srf = 0
    !clm_pf_idata%ngpf_srf = 0

    ! Allocate vectors for data transfer between CLM and PARFLOW. defined in
    ! parflow_dir/elm
    call CLMPARFLOWIDataCreateVec(MPI_COMM_WORLD)

  end subroutine CreateCLMPARFLOWInterfaceDate

  !-----------------------------------------------------------------------
  subroutine CreateCLMPARFLOWMaps(this, bounds, clm_npts, clm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Creates maps to transfer data between CLM and PARFLOW
    !
    ! !USES:
    use decompMod       , only : bounds_type, ldecomp
    use spmdMod     , only : mpicom, npes, iam
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(inout) :: clm_npts
    integer                 , intent(inout) :: clm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer             :: g,j
    integer             :: nlevmapped
    integer, pointer    :: clm_cell_ids_nindex(:)
    integer, pointer    :: clm_surf_cell_ids_nindex(:)

    ! Save cell IDs of CLM grid
    allocate(clm_cell_ids_nindex(     1:clm_npts     ))
    allocate(clm_surf_cell_ids_nindex(1:clm_surf_npts))

    nlevmapped     = clm_pf_idata%nzclm_mapped
    clm_npts       = 0
    clm_surf_npts  = 0
    do g = bounds%begg, bounds%endg
       do j = 1,nlevmapped
          clm_npts = clm_npts + 1
          clm_cell_ids_nindex(clm_npts) = (ldecomp%gdc2glo_rc(g)-1)*nlevmapped + j - 1
       enddo
       clm_surf_npts = clm_surf_npts + 1
       clm_surf_cell_ids_nindex(clm_surf_npts) = (ldecomp%gdc2glo_rc(g)-1)*nlevmapped
    enddo

    ! Initialize maps for transferring data between CLM and PARFLOW. Defined in
    ! parflow_dir/elm/parflow_model.F90
    call parflowModelInitMapping(this%parflow_m, clm_cell_ids_nindex, &
                                  clm_npts,CLM_SUB_TO_PF_SUB,mpicom,iam)
!    call parflowModelInitMapping(this%parflow_m, clm_cell_ids_nindex, &
!                                  clm_npts, nlmax_pf,CLM_SUB_TO_PF_EXTENDED_SUB)
    call parflowModelInitMapping(this%parflow_m, clm_cell_ids_nindex, &
                                  clm_npts, PF_SUB_TO_CLM_SUB,mpicom,iam)
    deallocate(clm_cell_ids_nindex)
    deallocate(clm_surf_cell_ids_nindex)

  end subroutine CreateCLMPARFLOWMaps


  !-----------------------------------------------------------------------
  subroutine extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    use landunit_varcon           , only : istsoil
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use clm_varcon                , only : denice, denh2o
    use clm_varpar                , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_parflow_type)              :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer               :: c,j,g,l,pf_j  ! do loop indices
    integer               :: nlevmapped
    integer               :: gcount
    integer               :: bounds_proc_begc, bounds_proc_endc

    integer     , pointer :: col_active(:)
    integer     , pointer :: col_landunit(:)
    integer     , pointer :: col_gridcell(:)
    integer     , pointer :: lun_type(:)

    real(r8)    , pointer :: l2e_h2osoi_ice(:,:)

    real(r8)    , pointer :: e2l_h2osoi_liq(:,:)
    real(r8)    , pointer :: e2l_h2osoi_ice(:,:)
    real(r8)    , pointer :: e2l_h2osoi_vol(:,:)
    real(r8)    , pointer :: e2l_zwt(:)
    real(r8)    , pointer :: e2l_mflx_snowlyr_col(:)
    real(r8)    , pointer :: e2l_watsatc(:,:)
    real(r8)    , pointer :: e2l_hksatc(:,:)
    real(r8)    , pointer :: e2l_bswc(:,:)
    real(r8)    , pointer :: e2l_sucsatc(:,:)

    real(r8)    , pointer :: dz(:,:)

    PetscScalar , pointer :: sat_clm_loc(:)
    PetscScalar , pointer :: watsat_clm_loc(:)
    PetscScalar , pointer :: hksat_clm_loc(:)
    PetscScalar , pointer :: bsw_clm_loc(:)
    PetscScalar , pointer :: sucsat_clm_loc(:)
    PetscErrorCode        :: ierr

    character(len= 128)   :: subname = 'extract_data_for_elm'
    !-----------------------------------------------------------------------

    bounds_proc_begc = bounds_clump%begc
    bounds_proc_endc = bounds_clump%endc
    nlevmapped       = clm_pf_idata%nzclm_mapped

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active             , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index     , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_gridcell_index     , col_gridcell )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz                , dz           )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type          , lun_type     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_vol      , e2l_h2osoi_vol       )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_watsatc     , e2l_watsatc )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_hksatc      , e2l_hksatc  )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_bswc        , e2l_bswc    )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_sucsatc     , e2l_sucsatc )
    
    call parflowModelGetSoilProp(this%parflow_m)
    call parflowModelGetSaturation(this%parflow_m,pf_sat,pf_grid_vol,pf_por)

    ! Set initial value of for ELM
    e2l_mflx_snowlyr_col(:) = 0._r8
    e2l_zwt(:)              = 0._r8

    ! Initialize soil moisture
    call VecGetArrayF90(clm_pf_idata%sat_clm      , sat_clm_loc    , ierr)
    call VecGetArrayF90(clm_pf_idata%watsat2_clm  , watsat_clm_loc , ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_x2_clm , hksat_clm_loc  , ierr)
    call VecGetArrayF90(clm_pf_idata%bsw2_clm     , bsw_clm_loc    , ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat2_clm  , sucsat_clm_loc , ierr)

    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          l = col_landunit(c)
          if (lun_type(l) == istsoil) then
             g = col_gridcell(c)
             gcount = g - bounds_clump%begg
             do j = 1, nlevgrnd
                pf_j = gcount*nlevmapped + j

                if (j <= nlevmapped) then
                   e2l_h2osoi_liq(c,j) = sat_clm_loc(pf_j)*watsat_clm_loc(pf_j)*dz(c,j)*1.e3_r8 !1000kg/m2 = m3/m2

                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                        l2e_h2osoi_ice(c,j)/dz(c,j)/denice
                   e2l_h2osoi_vol(c,j) = min(e2l_h2osoi_vol(c,j),watsat_clm_loc(pf_j))
                   e2l_h2osoi_ice(c,j) = 0._r8

                   e2l_watsatc(c,j) = watsat_clm_loc(pf_j)
                   e2l_hksatc(c,j)  = hksat_clm_loc(pf_j)
                   e2l_bswc(c,j)    = bsw_clm_loc(pf_j)
                   e2l_sucsatc(c,j) = sucsat_clm_loc(pf_j)
                else
                   e2l_h2osoi_liq(c,j) = e2l_h2osoi_liq(c,nlevmapped)
                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_vol(c,nlevmapped)
                   e2l_h2osoi_ice(c,j) = 0._r8
                   e2l_watsatc(c,j)    = e2l_watsatc(c,nlevmapped)
                   e2l_hksatc(c,j)     = e2l_hksatc(c,nlevmapped)
                   e2l_bswc(c,j)       = e2l_bswc(c,nlevmapped)
                   e2l_sucsatc(c,j)    = e2l_sucsatc(c,nlevmapped)
                end if

             enddo
          else
             write(iulog,*)'WARNING: Land Unit type other than soil type is present within the domain'
             call endrun( trim(subname)//' ERROR: Land Unit type not supported' )             
          endif
       endif
    enddo

    call VecRestoreArrayF90(clm_pf_idata%sat_clm      , sat_clm_loc    , ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat2_clm  , watsat_clm_loc , ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_x2_clm , hksat_clm_loc  , ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw2_clm     , bsw_clm_loc    , ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat2_clm  , sucsat_clm_loc , ierr)

   end subroutine extract_data_for_elm

   !------------------------------------------------------------------------
   subroutine EM_PARFLOW_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The VSFM dirver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    select case (em_stage)
    case (EM_PARFLOW_SOIL_HYDRO_STAGE)

       call EM_PARFLOW_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
            bounds_clump)

    case default
       write(iulog,*)'EM_FATES_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_PARFLOW_Solve

    !------------------------------------------------------------------------
  subroutine EM_PARFLOW_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! Solve the Variably Saturated Flow Model (VSFM) in soil columns.
    !
#include <petsc/finclude/petsc.h>
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use mpp_varpar                , only : nlevgrnd
    use clm_varcon                , only : denice, denh2o
    use petscsnes
    use elm_parflow_module
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_parflow_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g                                                       ! do loop indices
    integer                              :: pi                                                               ! pft index
    real(r8)                             :: dzsum                                                            ! summation of dzmm of layers below water table (mm)
    real(r8)                             :: dtime

    real(r8)  , pointer                  :: mflx_et_col_1d         (:)
    real(r8)  , pointer                  :: mflx_infl_col_1d       (:)
    real(r8)  , pointer                  :: mflx_dew_col_1d        (:)
    real(r8)  , pointer                  :: mflx_drain_col_1d      (:)
    real(r8)  , pointer                  :: mflx_sub_snow_col_1d   (:)
    real(r8)  , pointer                  :: mflx_snowlyr_col_1d    (:)
    real(r8)  , pointer                  :: t_soil_col_1d          (:)
    integer    , pointer                 :: col_active             (:)
    integer    , pointer                 :: col_gridcell           (:)

    real(r8)  , pointer                  :: fliq_col_1d       (:)
    real(r8)  , pointer                  :: mass_col_1d       (:)
    real(r8)  , pointer                  :: smpl_col_1d       (:)
    real(r8)  , pointer                  :: soilp_col_1d      (:)
    real(r8)  , pointer                  :: sat_col_1d        (:)

    real(r8)  , pointer                  :: frac_ice                    (:,:) ! fraction of ice
    real(r8)  , pointer                  :: total_mass_flux_col         (:)            ! Sum of all source-sinks conditions for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_et_col      (:)            ! ET sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_infl_col    (:)            ! Infiltration source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_dew_col     (:)            ! Dew source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_drain_col   (:)            ! Drainage sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_snowlyr_col (:)            ! Flux due to disappearance of snow for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_sub_col     (:)            ! Sublimation sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_lateral_col (:)            ! Lateral flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_seepage_col (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: qflx_seepage                (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: mass_prev_col          (:,:) ! Mass of water before a VSFM solve
    real(r8)  , pointer                  :: dmass_col              (:)            ! Change in mass of water after a VSFM solve
    real(r8)  , pointer                  :: mass_beg_col                (:)            ! Total mass before a VSFM solve
    real(r8)  , pointer                  :: mass_end_col                (:)            ! Total mass after a VSFM solve
    integer                              :: ier                                                              ! error status

    integer                              :: begc, endc
    integer                              :: g_idx, c_idx

    PetscInt                             :: soe_auxvar_id                                                    ! Index of system-of-equation's (SoE's) auxvar
    PetscErrorCode                       :: ierr                                                             ! PETSc return error code

    PetscBool                            :: converged                                                        ! Did VSFM solver converge to a solution with given PETSc SNES tolerances
    PetscInt                             :: converged_reason                                                 ! SNES converged due to which criteria
    PetscReal                            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal                            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal                            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt                             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt                             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal                            :: stol                                                             ! solution convergance tolerance
    PetscReal                            :: rtol                                                             ! relative convergance tolerance
    PetscReal,parameter                  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal                            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal                            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal                            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_lateral                                          ! Sum of lateral mass flux for all active soil columns
    PetscReal                            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt                             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter                  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt                             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt                             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal                            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter                 :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscBool                            :: successful_step                                                  ! Is the solution return by VSFM acceptable
    PetscReal , pointer                  :: soilp_col_ghosted_1d(:)
    PetscReal , pointer                  :: fliq_col_ghosted_1d(:)
    PetscReal , pointer                  :: mflx_lateral_col_1d(:)
    PetscReal , pointer                  :: lat_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_press_1d(:)

    integer                              :: jwt
    real(r8)                             :: z_dn, z_up

    real(r8)  , pointer                  :: l2e_mflux_infil(:)
    real(r8)  , pointer                  :: l2e_mflux_dew(:)
    real(r8)  , pointer                  :: l2e_mflux_sub_snow(:)
    real(r8)  , pointer                  :: l2e_mflux_snowlyr(:)
    real(r8)  , pointer                  :: l2e_mflux_et(:,:)
    real(r8)  , pointer                  :: l2e_mflux_drain(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: l2e_zi(:,:)
    real(r8)  , pointer                  :: col_dz(:,:)
    integer   , pointer                  :: l2e_filter_hydrologyc(:)
    integer                              :: l2e_num_hydrologyc
 
    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp(:,:)
    real(r8)  , pointer                  :: e2l_wtd(:)
    real(r8)  , pointer                  :: e2l_soilp(:,:)
    real(r8)  , pointer                  :: e2l_qrecharge(:)

    PetscScalar, pointer :: qflx_clm_loc(:)
    PetscScalar, pointer :: area_clm_loc(:)
    PetscScalar, pointer :: thetares2_clm_loc(:)
    PetscScalar, pointer :: watsat_clm_loc(:)
    PetscScalar, pointer :: sat_clm_loc(:)
    PetscScalar, pointer :: mass_clm_loc(:)
    PetscScalar, pointer :: e2l_drain_perched(:)
    PetscScalar, pointer :: e2l_drain(:)
    PetscScalar, pointer :: e2l_qrgwl(:)
    PetscScalar, pointer :: e2l_rsub_sat(:)

    integer :: bounds_proc_begc, bounds_proc_endc
    integer :: nlevmapped
    real(r8) :: col_wtgcell
    real(r8) :: curr_secs
    real(r8) :: pftime
    real(r8) :: pfdt
    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    ! Get time step

    dtime = dt

    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_infil       , l2e_mflux_infil       )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_dew         , l2e_mflux_dew         )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snow_sub    , l2e_mflux_sub_snow    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snowlyr     , l2e_mflux_snowlyr     )

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_et          , l2e_mflux_et          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_drainage    , l2e_mflux_drain       )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq , l2e_h2osoi_liq        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice , l2e_h2osoi_ice        )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_hydrologyc , l2e_filter_hydrologyc )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_hydrologyc   , l2e_num_hydrologyc    )

    call l2e_list%GetPointerToReal2D(this%index_l2e_column_zi        , l2e_zi                )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_active        , col_active            )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_gridcell      , col_gridcell          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_col_dz           , col_dz                )

    call e2l_list%GetPointerToReal1D(this%index_e2l_state_wtd        , e2l_wtd               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_liq , e2l_h2osoi_liq        )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_ice , e2l_h2osoi_ice        )
    !call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp        , e2l_smp               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_soilp      , e2l_soilp             )

    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrecharge    , e2l_qrecharge        )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_drain_perched, e2l_drain_perched    )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_drain        , e2l_drain            )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrgwl        , e2l_qrgwl            )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_rsub_sat     , e2l_rsub_sat         )

    begc = bounds_proc_begc
    endc = bounds_proc_endc

    allocate(frac_ice                    (begc:endc,1:nlevgrnd))
    allocate(total_mass_flux_col         (begc:endc))
    allocate(total_mass_flux_et_col      (begc:endc))
    allocate(total_mass_flux_infl_col    (begc:endc))
    allocate(total_mass_flux_dew_col     (begc:endc))
    allocate(total_mass_flux_drain_col   (begc:endc))
    allocate(total_mass_flux_snowlyr_col (begc:endc))
    allocate(total_mass_flux_sub_col     (begc:endc))
    allocate(total_mass_flux_lateral_col (begc:endc))
    allocate(total_mass_flux_seepage_col (begc:endc))
    allocate(qflx_seepage                (begc:endc))
    allocate(mass_prev_col          (begc:endc,1:nlevgrnd))
    allocate(dmass_col              (begc:endc))
    allocate(mass_beg_col                (begc:endc))
    allocate(mass_end_col                (begc:endc))

    allocate(mflx_et_col_1d              ((endc-begc+1)*nlevgrnd))
    allocate(mflx_drain_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(mflx_infl_col_1d            (endc-begc+1))
    allocate(mflx_dew_col_1d             (endc-begc+1))
    allocate(mflx_sub_snow_col_1d        (endc-begc+1))
    allocate(mflx_snowlyr_col_1d         (endc-begc+1))
    allocate(t_soil_col_1d               ((endc-begc+1)*nlevgrnd))

    allocate(mass_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(fliq_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(smpl_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(soilp_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(sat_col_1d             ((endc-begc+1)*nlevgrnd))

    ! initialize
    mflx_et_col_1d(:)                = 0.d0
    mflx_infl_col_1d(:)              = 0.d0
    mflx_dew_col_1d(:)               = 0.d0
    mflx_drain_col_1d(:)             = 0.d0
    mflx_sub_snow_col_1d(:)          = 0.d0
    mflx_snowlyr_col_1d(:)           = 0.d0
    t_soil_col_1d(:)                 = 298.15d0

    mass_beg                         = 0.d0
    mass_end                         = 0.d0
    total_mass_flux                  = 0.d0
    total_mass_flux_et               = 0.d0
    total_mass_flux_infl             = 0.d0
    total_mass_flux_dew              = 0.d0
    total_mass_flux_drain            = 0.d0
    total_mass_flux_snowlyr          = 0.d0
    total_mass_flux_sub              = 0.d0
    total_mass_flux_lateral          = 0.d0

    mass_beg_col(:)                  = 0.d0
    mass_end_col(:)                  = 0.d0
    total_mass_flux_col(:)           = 0.d0
    total_mass_flux_et_col(:)        = 0.d0
    total_mass_flux_infl_col(:)      = 0.d0
    total_mass_flux_dew_col(:)       = 0.d0
    total_mass_flux_drain_col(:)     = 0.d0
    total_mass_flux_snowlyr_col(:)   = 0.d0
    total_mass_flux_sub_col(:)       = 0.d0
    total_mass_flux_lateral_col(:)   = 0.d0

    mass_prev_col(:,:)          = 0.d0
    dmass_col(:)                = 0.d0

    nlevmapped = clm_pf_idata%nzclm_mapped

    call parflowModelGetSaturation(this%parflow_m,pf_sat,pf_grid_vol,pf_por)

    ! Get total mass
!    call parflowModelGetUpdatedData( this%parflow_m )

    call VecGetArrayF90(clm_pf_idata%mass_clm  , mass_clm_loc  , ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%area_top_face_clm, area_clm_loc, ierr); CHKERRQ(ierr)

    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       g = col_gridcell(c)

       do j = 1, nlevmapped

          c_idx = (c - begc)*nlevgrnd + j
          g_idx = (g - bounds_clump%begg)*nlevmapped + j

          mflx_et_col_1d(c_idx)          = l2e_mflux_et(c,j)
          mflx_drain_col_1d(c_idx)       = l2e_mflux_drain(c,j)

          total_mass_flux_et           = total_mass_flux_et           + mflx_et_col_1d(c_idx)
          total_mass_flux_et_col(c)    = total_mass_flux_et_col(c)    + mflx_et_col_1d(c_idx)

          total_mass_flux_drain        = total_mass_flux_drain        + mflx_drain_col_1d(c_idx)
          total_mass_flux_drain_col(c) = total_mass_flux_drain_col(c) + mflx_drain_col_1d(c_idx)

          mass_beg                     = mass_beg                     + mass_clm_loc(g_idx)/area_clm_loc(g_idx)
          mass_beg_col(c)              = mass_beg_col(c)              + mass_clm_loc(g_idx)/area_clm_loc(g_idx)
          mass_prev_col(c,j)      = mass_col_1d(c_idx)
       end do

       c_idx = c - begc+1

       mflx_dew_col_1d(c_idx)        = l2e_mflux_dew(c)
       mflx_infl_col_1d(c_idx)       = l2e_mflux_infil(c)
       mflx_snowlyr_col_1d(c_idx)    = l2e_mflux_snowlyr(c)
       mflx_sub_snow_col_1d(c_idx)   = l2e_mflux_sub_snow(c)

       total_mass_flux_dew            = total_mass_flux_dew            + mflx_dew_col_1d(c_idx)
       total_mass_flux_dew_col(c)     = total_mass_flux_dew_col(c)     + mflx_dew_col_1d(c_idx)

       total_mass_flux_infl           = total_mass_flux_infl           + mflx_infl_col_1d(c_idx)
       total_mass_flux_infl_col(c)    = total_mass_flux_infl_col(c)    + mflx_infl_col_1d(c_idx)

       total_mass_flux_snowlyr        = total_mass_flux_snowlyr        + mflx_snowlyr_col_1d(c_idx)
       total_mass_flux_snowlyr_col(c) = total_mass_flux_snowlyr_col(c) + mflx_snowlyr_col_1d(c_idx)

       total_mass_flux_sub            = total_mass_flux_sub            + mflx_sub_snow_col_1d(c_idx)
       total_mass_flux_sub_col(c)     = total_mass_flux_sub_col(c)     + mflx_sub_snow_col_1d(c_idx)

       total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
            total_mass_flux_infl_col(c)    + &
            total_mass_flux_dew_col(c)     + &
            total_mass_flux_drain_col(c)   + &
            total_mass_flux_snowlyr_col(c) + &
            total_mass_flux_sub_col(c)     + &
            total_mass_flux_lateral_col(c)
    end do
    total_mass_flux        = &
         total_mass_flux_et        + &
         total_mass_flux_infl      + &
         total_mass_flux_dew       + &
         total_mass_flux_drain     + &
         total_mass_flux_snowlyr   + &
         total_mass_flux_sub       + &
         total_mass_flux_lateral
    call VecRestoreArrayF90(clm_pf_idata%mass_clm  , mass_clm_loc  , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clm, area_clm_loc, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(clm_pf_idata%qflx_clm, qflx_clm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%area_top_face_clm, area_clm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%thetares2_clm, thetares2_clm_loc, ierr); CHKERRQ(ierr)

    frac_ice(:,:)       = 0.d0
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       do j = 1, nlevmapped
          frac_ice(c,j) = l2e_h2osoi_ice(c,j)/(l2e_h2osoi_liq(c,j) + l2e_h2osoi_ice(c,j))
       end do
    end do

    ! Initialize ET sink to ZERO
    do g = bounds_clump%begg, bounds_clump%endg
       do j = 1,nlevmapped
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          qflx_clm_loc(g_idx) = 0.0_r8
       end do
    end do

    ! Account for following fluxes in the top soil layer:
    ! - infiltration
    ! - dew
    ! - disapperance of snow layer
    ! - sublimation of snow
    j = 1
    col_wtgcell = 1._r8
    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          ! Set gridcell indices
          g = col_gridcell(c)
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          c_idx = c - begc+1
          qflx_clm_loc(g_idx) = qflx_clm_loc(g_idx) + &
               (&
               mflx_infl_col_1d(c_idx)    + &
               mflx_dew_col_1d(c_idx)     + &
               mflx_snowlyr_col_1d(c_idx) + &
               mflx_sub_snow_col_1d(c_idx)  &
               )*col_wtgcell*area_clm_loc(g_idx)
       end if
       total_mass_flux_col(c) = 0.d0
    enddo

    ! Account for evapotranspiration flux
    do c = bounds_proc_begc, bounds_proc_endc
       do j = 1,nlevmapped
          g = col_gridcell(c)
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          c_idx = (c - bounds_proc_begc)*nlevgrnd+j
          if (col_active(c) == 1) then
             qflx_clm_loc(g_idx) = qflx_clm_loc(g_idx) + &
                  mflx_et_col_1d(c_idx)*area_clm_loc(g_idx) !kg/s*m2
             total_mass_flux_col(c) = total_mass_flux_col(c) + qflx_clm_loc(g_idx)/area_clm_loc(g_idx) !kg/s
 !  convert unit to 1/hr, already done in soilwatermovement.F90
 !  PARFLOW needs units of 1/h for PF_FLUX so divide by dz  
 !             qflx_clm_loc(g_idx) = qflx_clm_loc(g_idx)*3.6d0/col_dz(c,j)/(denh2o * 1.0d-3)/area_clm_loc(g_idx) !1/hr
             qflx_clm_loc(g_idx) = qflx_clm_loc(g_idx)/area_clm_loc(g_idx) !kg/s
          end if
       end do
    end do
    call VecRestoreArrayF90(clm_pf_idata%qflx_clm, qflx_clm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%thetares2_clm, thetares2_clm_loc, ierr); CHKERRQ(ierr)

 !   call parflowModelUpdateFlowConds( this%parflow_m )
 !   call parflowModelStepperRunTillPauseTime( this%parflow_m, (nstep+1.0d0)*dtime )
    call parflowModelUpdateSourceSink(this%parflow_m,elm_flux)
 !
    curr_secs = nstep*dtime
    pftime = curr_secs/3600.d0
    pfdt = dtime/3600.d0
    ! convert unit from kg/s (mm/s) to 1/hr
!    if(sum(abs(elm_flux)) .ne. 0.d0) then
!     print*,'dz-',pf_grid_dz(1:nlevmapped)
!     print *,'elm-flux-',elm_flux(1:nlevmapped)
!stop
!    endif
    elm_flux(:) = elm_flux(:) * 3600.d0 * 1.0d-3 / pf_grid_dz(:) 
    call elmparflowadvance(pftime,pfdt,elm_flux,pf_press,pf_porosity,pf_sat,nlevmapped, &
                           0,0,0,0)


    call VecGetArrayF90(clm_pf_idata%sat_clm   , sat_clm_loc   , ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%mass_clm  , mass_clm_loc  , ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr); CHKERRQ(ierr)

    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       g = col_gridcell(c)

       ! initialization
       jwt = -1

       ! Loops in decreasing j so WTD can be computed in the same loop
       e2l_h2osoi_liq(c,:) = 0._r8
       e2l_h2osoi_ice(c,:) = 0._r8
       do j = nlevmapped, 1, -1
          g_idx = (g - bounds_clump%begg)*nlevmapped + j

          e2l_h2osoi_liq(c,j) =  (1.d0 - frac_ice(c,j))*mass_clm_loc(g_idx)/area_clm_loc(g_idx) !mm = kg/m2
          e2l_h2osoi_ice(c,j) =  frac_ice(c,j)         *mass_clm_loc(g_idx)/area_clm_loc(g_idx)

          mass_end        = mass_end        + mass_clm_loc(g_idx)/area_clm_loc(g_idx)
          mass_end_col(c) = mass_end_col(c) + mass_clm_loc(g_idx)/area_clm_loc(g_idx)

       end do

       ! Find maximum water balance error over the column
       !abs_mass_error_col = max(abs_mass_error_col,                     &
       !     abs(mass_beg_col(c) - mass_end_col(c) + &
       !     total_mass_flux_col(c)*dt))
       e2l_qrecharge     (c) = 0._r8

       e2l_wtd(c) = l2e_zi(c,nlevmapped)
    end do

    ! Save soil liquid pressure from VSFM for all (active+nonactive) cells.
    ! soilp_col is used for restarting VSFM.
    do c = begc, endc
       do j = 1, nlevgrnd
          c_idx = (c - begc)*nlevgrnd + j
          e2l_soilp(c,j) = 0._r8!soilp_col_1d(c_idx)
       end do
    end do

    e2l_drain_perched (:) = 0._r8
    e2l_drain         (:) = 0._r8
    e2l_qrgwl         (:) = 0._r8
    e2l_rsub_sat      (:) = 0._r8
    
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clm, area_clm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%sat_clm   , sat_clm_loc   , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%mass_clm  , mass_clm_loc  , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr); CHKERRQ(ierr)


    deallocate(frac_ice                    )
    deallocate(total_mass_flux_col         )
    deallocate(total_mass_flux_et_col      )
    deallocate(total_mass_flux_infl_col    )
    deallocate(total_mass_flux_dew_col     )
    deallocate(total_mass_flux_drain_col   )
    deallocate(total_mass_flux_snowlyr_col )
    deallocate(total_mass_flux_sub_col     )
    deallocate(total_mass_flux_lateral_col )
    deallocate(total_mass_flux_seepage_col )
    deallocate(qflx_seepage                )
    deallocate(mass_prev_col          )
    deallocate(dmass_col              )
    deallocate(mass_beg_col                )
    deallocate(mass_end_col                )

    deallocate(mflx_et_col_1d              )
    deallocate(mflx_drain_col_1d           )
    deallocate(mflx_infl_col_1d            )
    deallocate(mflx_dew_col_1d             )
    deallocate(mflx_sub_snow_col_1d        )
    deallocate(mflx_snowlyr_col_1d         )
    deallocate(t_soil_col_1d               )

    deallocate(mass_col_1d            )
    deallocate(fliq_col_1d            )
    deallocate(smpl_col_1d            )
    deallocate(soilp_col_1d           )
    deallocate(sat_col_1d             )

  end subroutine EM_PARFLOW_Solve_Soil_Hydro
!..................................................................
#endif
end module ExternalModelPARFLOWMod
