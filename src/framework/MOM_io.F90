module MOM_io

! This file is part of MOM6. See LICENSE.md for the license.


use MOM_error_handler,    only : MOM_error, NOTE, FATAL, WARNING
use MOM_domains,          only : MOM_domain_type, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,          only : get_simple_array_i_ind, get_simple_array_j_ind
use MOM_file_parser,      only : log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_dyn_horgrid,      only : dyn_horgrid_type
use MOM_string_functions, only : lowercase, slasher
use MOM_string_functions, only : append_substring
use MOM_time_manager,     only : time_type, time_type_to_real
use MOM_verticalGrid,     only : verticalGrid_type
use ensemble_manager_mod, only : get_ensemble_id
use fms_mod,              only : write_version_number, open_namelist_file, check_nml_error
use MOM_string_functions, only : extract_word
use mpp_mod,              only : mpp_max, mpp_pe, mpp_npes
use mpp_domains_mod,      only : domain1d, domain2d, domainug, mpp_get_domain_components, &
                                 mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
use mpp_domains_mod,      only : mpp_get_domain_npes, mpp_define_io_domain, mpp_get_io_domain, mpp_get_io_domain_layout
use mpp_domains_mod,      only : CENTER, CORNER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_io_mod,           only : open_file => mpp_open, close_file => mpp_close
use mpp_io_mod,           only : mpp_write_meta
use mpp_io_mod,           only : axistype
use mpp_io_mod,           only : fieldtype, axistype, flush_file => mpp_flush
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, NETCDF_FILE=>MPP_NETCDF
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : SINGLE_FILE=>MPP_SINGLE, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod,           only : MPP_APPEND, MPP_MULTI, MPP_OVERWR, MPP_NETCDF, MPP_RDONLY
use mpp_io_mod,           only : io_infra_init=>mpp_io_init
use fms2_io_mod,          only : check_if_open, &
                                get_dimension_names, &
                                get_dimension_size, &
                                get_compute_domain_dimension_indices, &
                                get_global_attribute, &
                                get_global_io_domain_indices, &
                                get_num_dimensions, &
                                get_num_variables, &
                                get_unlimited_dimension_name, &
                                get_variable_dimension_names, &
                                get_variable_names, &
                                get_variable_num_dimensions, &
                                get_variable_size, &
                                get_variable_units, &
                                get_variable_unlimited_dimension_index, &
                                global_att_exists, &
                                is_dimension_unlimited, &
                                is_dimension_registered, &
                                read_data, &
                                read_restart, &
                                register_restart_field, &
                                register_axis, &
                                register_field, &
                                register_variable_attribute, &
                                fms2_open_file => open_file, &
                                fms2_close_file => close_file, &
                                write_data, &
                                write_restart, &
                                attribute_exists => variable_att_exists, &
                                variable_exists, &
                                dimension_exists, &
                                file_exists, &
                                FmsNetcdfDomainFile_t, &
                                FmsNetcdfFile_t, &
                                FmsNetcdfUnstructuredDomainFile_t, &
                                unlimited

use netcdf

implicit none ; private

public :: close_file, open_file, fieldtype, flush_file, write_version_number
public :: num_timelevels, MOM_read_vector, ensembler, create_file
public :: slasher, MOM_io_init, field_exists, field_size, read_axis_data
public :: open_namelist_file, check_nml_error, io_infra_init
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE
public :: var_desc, modify_vardesc, query_vardesc, cmor_long_std
public :: scale_data, convert_checksum_to_string
! new FMS-IO routines and wrappers
public :: attribute_exists
public :: check_if_open
public :: dimension_exists
public :: file_exists
public :: fms2_open_file
public :: fms2_close_file
public :: FmsNetcdfFile_t
public :: FmsNetcdfDomainFile_t
public :: FmsNetcdfUnstructuredDomainFile_t
public :: get_compute_domain_dimension_indices
public :: get_dimension_names
public :: get_dimension_size
public :: get_global_io_domain_indices
public :: get_global_attribute
public :: get_horizontal_grid_logic
public :: get_num_dimensions
public :: get_num_variables
public :: get_time_units
public :: get_unlimited_dimension_name
public :: get_var_dimension_features
public :: get_variable_dimension_names
public :: get_variable_byte_size
public :: get_variable_num_dimensions
public :: get_variable_size
public :: get_variable_units
public :: get_variable_unlimited_dimension_index
public :: global_att_exists
public :: is_dimension_registered
public :: is_dimension_unlimited
public :: MOM_get_diagnostic_axis_data
public :: MOM_read_data
public :: MOM_register_diagnostic_axis
public :: MOM_register_variable_axes
public :: read_data
public :: read_restart
public :: register_axis
public :: register_field
public :: register_restart_field
public :: register_variable_attribute
public :: variable_exists
public :: write_data
public :: write_field
public :: write_restart
public :: unlimited

!> Type for describing a variable, typically a tracer
type, public :: vardesc
  character(len=64)  :: name               !< Variable name in a NetCDF file
  character(len=48)  :: units              !< Physical dimensions of the variable
  character(len=240) :: longname           !< Long name of the variable
  character(len=8)   :: hor_grid           !< Horizontal grid:  u, v, h, q, Cu, Cv, T, Bu, or 1
  character(len=8)   :: z_grid             !< Vertical grid:  L, i, or 1
  character(len=8)   :: t_grid             !< Time description: s, p, or 1
  character(len=64)  :: cmor_field_name    !< CMOR name
  character(len=64)  :: cmor_units         !< CMOR physical dimensions of the variable
  character(len=240) :: cmor_longname      !< CMOR long name of the variable
  real               :: conversion         !< for unit conversions, such as needed to
                                           !! convert from intensive to extensive
end type vardesc

!> A type for making arrays of pointers to real 1-d arrays
type p1d
  real, dimension(:), pointer :: p => NULL() !< A pointer to a 1d array
end type p1d

!> A structure with information about a single axis variable
type axis_atts
  character(len=64)  :: name                    !< Names of the axis
  character(len=48)  :: units                   !< Physical dimensions of the axis
  character(len=240) :: longname                !< Long name of the axis
  character(len=8)   :: positive                !< Positive-definite direction: up, down, east, west, north, south
  integer            :: horgrid_position        !< Horizontal grid position
  logical            :: is_domain_decomposed    !< if .true. the axis data are domain-decomposed
                                                !! and need to be indexed by the compute domain
                                                !! before passing to write_data
end type axis_atts

!> Type for describing an axis variable (e.g., lath, lonh, Time)
type, public :: axis_data_type
  !> An array of descriptions of the registered axes
  type(axis_atts), pointer :: axis(:) => NULL()  !< structure with axis attributes
  type(p1d), pointer       :: data(:) => NULL()  !< pointer to the axis data
end type axis_data_type

!> interface to read data from a netcdf file
interface MOM_read_data
  module procedure MOM_read_data_4d_DD
  module procedure MOM_read_data_3d_DD
  module procedure MOM_read_data_2d_DD
  module procedure MOM_read_data_1d_DD
  module procedure MOM_read_data_scalar
  module procedure MOM_read_data_4d_noDD
  module procedure MOM_read_data_3d_noDD
  module procedure MOM_read_data_2d_noDD
  module procedure MOM_read_data_1d_noDD
  module procedure MOM_read_data_2d_noDD_diag_axes
  module procedure MOM_read_data_2d_supergrid
end interface

!> Read a pair of data fields representing the two components of a vector from a netcdf file
interface MOM_read_vector
  module procedure MOM_read_vector_3d
  module procedure MOM_read_vector_2d
end interface

!> interface to write data to a netcdf file generated by create_file
interface write_field
  module procedure write_field_4d_DD
  module procedure write_field_3d_DD
  module procedure write_field_2d_DD
  module procedure write_field_1d_DD
  module procedure write_scalar
  module procedure write_field_4d_noDD
  module procedure write_field_3d_noDD
  module procedure write_field_2d_noDD
  module procedure write_field_1d_noDD
end interface

!> interface to scale data after reading in a field
interface scale_data
  module procedure scale_data_4d
  module procedure scale_data_3d
  module procedure scale_data_2d
  module procedure scale_data_1d
end interface

!> interface for creating files to write/overwrite/append to
interface create_file
  module procedure create_file_filename
  module procedure create_file_fileobj_dd
end interface create_file

!> interface for registering axes associated with a variable to a netCDF file object 
interface MOM_register_variable_axes
  module procedure MOM_register_variable_axes_subdomain
  module procedure MOM_register_variable_axes_full
end interface MOM_register_variable_axes

!>\note CAUTION: The following variables are saved by default, and are only necessary for consecutive calls to
!! MOM_read_data or write_field with the same file name. The user should ensure that fms2_close_file on
!! the fileobj_write_field and fileobj_read structures are called at every requisite time step at after the last
!! variable is written to the file by omitting the optional leave_file_open argument, or setting it to .false.

type(FmsNetcdfFile_t), private :: fileobj_write_field !> netCDF non-domain-decomposed file object returned by call to
                                                      !! open_file in write_field calls.
type(FmsNetcdfDomainFile_t), private :: fileobj_write_field_dd !> netCDF domain-decomposed file object returned by call
                                                               !! to open_file in write_field calls.
type(FmsNetcdfDomainFile_t), private :: fileobj_read_dd !> netCDF domain-decomposed file object returned by call to
                                                        !! open_file in MOM_read_data_DD calls
type(FmsNetcdfFile_t), private :: fileobj_read !> netCDF domain-decomposed file object returned by call to
                                               !! open_file in MOM_read_data_noDD calls
!> Type with variable metadata for a netCDF file opened to read domain-decomposed data
type file_variable_meta_DD
  integer :: nvars = 0!> number of variables in a netCDF file opened to read domain-decomposed data
  character(len=96), allocatable, dimension(:) :: var_names !> array for names of variables in a netCDF
                                                            !! file opened to read domain-decomposed data
end type file_variable_meta_DD

!> Type with variable metadata for a netCDF file opened to read non-domain-decomposed data
type file_variable_meta_noDD
  integer :: nvars = 0 !> number of variables in a netCDF file opened to read non-domain-decomposed data
  character(len=96), allocatable, dimension(:) :: var_names !> array for names of variables in a netCDF
                                                            !! file opened to read non-domain-decomposed data
end type file_variable_meta_noDD

type (file_variable_meta_DD), private :: file_var_meta_DD
type (file_variable_meta_noDD), private :: file_var_meta_noDD
integer, private :: write_field_time_index !> index of the time_level value that is written to netCDF file by the
                                           !! write_field routines.

contains

!> This routine opens a netcdf file in "write" or "overwrite" mode, registers the global diagnostic axes, and writes
!! the axis data and metadata to the file
subroutine create_file_filename(filename, vars, numVariables, register_time, G, DG, GV, checksums, is_restart)
  character(len=*),      intent(in)               :: filename !< full path to the netcdf file
  type(vardesc), dimension(:), intent(in)         :: vars !< structures describing the output
  integer,               intent(in)               :: numVariables !< number of variables to write to the file
  logical, optional, intent(in) :: register_time !< if .true., register a time dimension to the file
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums(:,:)  !< checksums of the variables
  logical, optional, intent(in) :: is_restart !< indicates whether file is a restart file

  ! local
  type(FmsNetcdfFile_t) :: fileObjNoDD ! non-domain-decomposed netcdf file object returned by open_file
  type(FmsNetcdfDomainFile_t) :: fileObjDD ! domain-decomposed netcdf file object returned by open_file
  type(axis_data_type) :: axis_data_CS ! structure for coordinate variable metadata
  type(MOM_domain_type), pointer :: Domain => NULL()
  logical :: file_open_successDD, file_open_successNoDD ! true if netcdf file is opened
  logical :: one_file, domain_set ! indicates whether the file will be domain-decomposed or not
  logical :: reg_time ! register the time if .true.
  logical :: is_restart_file
  character(len=10) :: nc_mode
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=1024) :: filename_temp
  character(len=48), allocatable, dimension(:,:) :: dim_names ! variable dimension names
  integer :: i, is, ie, j, substring_index, total_axes
  integer :: num_dims ! number of dimensions
  integer :: thread ! indicates whether threading is used
  integer, dimension(4) :: dim_lengths ! variable dimension lengths
  integer, allocatable :: pelist(:) ! list of pes associated with the file
  real :: time

  ! determine whether the file will be domain-decomposed or not
  domain_set=.false.
  if (present(G)) then
    domain_set = .true. ; Domain => G%Domain
  elseif (present(dG)) then
    domain_set = .true. ; Domain => dG%Domain
  endif

  is_restart_file = .false.
  if (present(is_restart)) is_restart_file = is_restart
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index < 1) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  nc_mode = ""
  if (file_exists(trim(filename_temp))) then
    nc_mode = "overwrite"
  else
    nc_mode = "write"
  endif

  reg_time = .false.
  if (present(register_time)) reg_time = register_time

  ! open the file
  file_open_successNoDD=.false.
  file_open_successDD=.false.

  if (domain_set) then
    ! define the io domain if on one pe and the io domain is not set
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif

    if (.not. check_if_open(fileObjDD)) &
      file_open_successDD=fms2_open_file(fileObjDD, filename_temp, trim(nc_mode), Domain%mpp_domain, &
                                          is_restart=is_restart_file)
  else
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    allocate(pelist(mpp_npes()))
    pelist(:) = 0
    do i=1,size(pelist)
      pelist(i) = i-1
    enddo

    if (.not. check_if_open(fileObjNoDD)) &
      file_open_successNoDD=fms2_open_file(fileObjNoDD, filename_temp, trim(nc_mode), &
                                           is_restart=is_restart_file, pelist=pelist)
  endif
  ! allocate the output data variable dimension attributes
  allocate(dim_names(numVariables,4))
  dim_names(:,:) = ""
  ! allocate the axis data and attribute types for the file
  !> \note The user should increase the sizes of the axis and data attributes to accommodate more axes if necessary.
  allocate(axis_data_CS%axis(7))
  allocate(axis_data_CS%data(7))
  ! axis registration procedure for the domain-decomposed case
  if (file_open_successDD) then
    do i=1,numVariables
      num_dims=0
      dim_lengths(:) = 0

      !> \note The time dimension is registered separately at the end of the procedure if reg_time = .true.
      !! so the t_grid argument in get_var_dimension_features is set to '1' (do nothing)
      if (present(G)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                          dim_lengths, num_dims, G=G)
      elseif(present(dG)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                          dim_lengths, num_dims, dG=dG)
      endif

      if(present(GV)) &
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                        dim_lengths, num_dims, GV=GV)
      !> \note num_dims will be 0 for scalar values
      if (num_dims .le. 0) cycle

      do j=1,num_dims
        ! register the variable axes to the file if they are not already registered
        if (dim_lengths(j) .gt. 0) then
          if (.not.(dimension_exists(fileObjDD, dim_names(i,j)))) then

            if (present(G)) then
              if (present(GV)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, GV=GV)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G)
              endif   
            elseif (present(dG)) then
              if (present(GV)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, GV=GV)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG)
              endif
            elseif (present(GV)) then
               call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV)
            endif
              call MOM_register_diagnostic_axis(fileObjDD, trim(dim_names(i,j)), dim_lengths(j))
          endif
          ! register the axis attributes and write the axis data to the file
          if (.not.(variable_exists(fileObjDD, trim(axis_data_CS%axis(j)%name)))) then
            if (associated(axis_data_CS%data(j)%p)) then

                call register_field(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                  "double", dimensions=(/trim(axis_data_CS%axis(j)%name)/))

                call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                               'long_name', axis_data_CS%axis(j)%longname)

                call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                               'units', trim(axis_data_CS%axis(j)%units))

                if (len_trim(axis_data_CS%axis(j)%positive)>1) &
                  call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                                   'positive', trim(axis_data_CS%axis(j)%positive))

                if (axis_data_CS%axis(j)%is_domain_decomposed) then
                  call get_global_io_domain_indices(fileObjDD, trim(axis_data_CS%axis(j)%name), is, ie)
                  call write_data(fileObjDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p(is:ie))
                else
                  call write_data(fileObjDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p)
                endif
            endif
          endif
        endif
      enddo
    enddo

    if (reg_time) then
     if (.not.(dimension_exists(fileObjDD,"Time"))) &
       call register_axis(fileObjDD, "Time", unlimited)
    endif

    if (check_if_open(fileObjDD)) call fms2_close_file(fileObjDD)
  ! axis registration and write procedure for the non-domain-decomposed case
  elseif (file_open_successNoDD) then
    do i=1,numVariables
      num_dims=0
      dim_lengths(:) = 0

      !> \note The time dimension is registered separately at the end of the procedure if reg_time = .true.
      !! so the t_grid argument in get_var_dimension_features is set to '1' (do nothing)
      if (present(G)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                          dim_lengths, num_dims, G=G)
      elseif(present(dG)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                          dim_lengths, num_dims, dG=dG)
      endif

      if(present(GV)) &
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                        dim_lengths, num_dims, GV=GV)
      !> \note num_dims will be 0 for scalar variables
      if (num_dims .le. 0) cycle

      do j=1,num_dims
        ! register the variable axes to the file if they are not already registered
        if (dim_lengths(j) .gt. 0) then
          if (.not.(dimension_exists(fileObjNoDD, dim_names(i,j)))) then
            if (present(G)) then
              if (present(GV)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, GV=GV)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G)
              endif   
            elseif (present(dG)) then
              if (present(GV)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, GV=GV)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG)
              endif
            elseif (present(GV)) then
               call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV)
            endif
            call register_axis(fileObjNoDD, trim(dim_names(i,j)), dim_lengths(j))
          endif
          ! register the axis attributes and write the axis data to the file
          if (.not.(variable_exists(fileObjNoDD, trim(axis_data_CS%axis(j)%name)))) then
            if (associated(axis_data_CS%data(j)%p)) then
                call register_field(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                  "double", dimensions=(/trim(axis_data_CS%axis(j)%name)/))

                call register_variable_attribute(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                               'long_name', axis_data_CS%axis(j)%longname)

                call register_variable_attribute(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                               'units', trim(axis_data_CS%axis(j)%units))

                 if (len_trim(axis_data_CS%axis(j)%positive)>1) &
                   call register_variable_attribute(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                                    'positive', trim(axis_data_CS%axis(j)%positive))

                if (lowercase(trim(axis_data_CS%axis(j)%name)) .ne. 'time') then
                  call write_data(fileObjNoDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p)
                endif
            endif
          endif
        endif
      enddo
    enddo

    if (reg_time) then
     if (.not.(dimension_exists(fileObjNoDD,"Time"))) &
       call register_axis(fileObjNoDD, "Time" , unlimited)
    endif

    if (check_if_open(fileObjNoDD)) call fms2_close_file(fileObjNoDD)
  endif

  deallocate(dim_names)
  deallocate(axis_data_CS%axis)
  deallocate(axis_data_CS%data)
  if (allocated(pelist)) deallocate(pelist)
  nullify(Domain)

end subroutine create_file_filename

!> This routine opens a netcdf file in "write" or "overwrite" mode, registers the global diagnostic axes, and writes
!! the axis data and metadata to the file. It returns the netcdf file object for additional writing.
subroutine create_file_fileobj_dd(filename, fileObjDD, vars, numVariables, register_time, G, DG, GV, checksums, is_restart)
  character(len=*),      intent(in)               :: filename !< full path to the netcdf file
  type(FmsNetcdfDomainFile_t), intent(inout)      :: fileObjDD !< domain-decomposed netcdf file object
                                                               !! returned by open_file
  type(vardesc), dimension(:), intent(in)         :: vars !< structures describing the output
  integer,               intent(in)               :: numVariables !< number of variables to write to the file
  logical, optional, intent(in) :: register_time !< if .true., register a time dimension to the file
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums(:,:)  !< checksums of the variables
  logical, optional, intent(in) :: is_restart !< indicates whether file is a restart file

  ! local
  type(axis_data_type) :: axis_data_CS ! structure for coordinate variable metadata
  type(MOM_domain_type), pointer :: Domain => NULL()
  logical :: file_open_successDD ! true if netcdf file is opened
  logical :: one_file, domain_set ! indicates whether the file will be domain-decomposed or not
  logical :: reg_time ! register the time if .true.
  logical :: is_restart_file ! .true. if the file is a restart file
  character(len=10) :: nc_mode
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=1024) :: filename_temp
  character(len=48), allocatable, dimension(:,:) :: dim_names ! variable dimension names
  integer :: i, is, ie, j, substring_index, total_axes
  integer :: num_dims ! number of dimensions
  integer :: thread ! indicates whether threading is used
  integer, dimension(4) :: dim_lengths ! variable dimension lengths
  integer, allocatable :: pelist(:) ! list of pes associated with the file
  real :: time

  ! determine whether the file will be domain-decomposed or not
  domain_set=.false.
  if (present(G)) then
    domain_set = .true. ; Domain => G%Domain
  elseif (present(dG)) then
    domain_set = .true. ; Domain => dG%Domain
  endif

  is_restart_file = .false.
  if (present(is_restart)) is_restart_file = is_restart
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index < 1) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  nc_mode = ""
  if (file_exists(trim(filename_temp))) then
    nc_mode = "overwrite"
  else
    nc_mode = "write"
  endif

  reg_time = .false.
  if (present(register_time)) reg_time = register_time
  ! open the file
  file_open_successDD=.false.
  ! define the io domain if on one pe and the io domain is not set
  if (domain_set) then
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
  else
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    allocate(pelist(mpp_npes()))
    pelist(:) = 0
    do i=1,size(pelist)
      pelist(i) = i-1
    enddo
  endif
  if (.not. check_if_open(fileObjDD)) &
    file_open_successDD=fms2_open_file(fileObjDD, filename_temp, trim(nc_mode), Domain%mpp_domain, &
                                       is_restart=is_restart_file)
  ! allocate the output data variable dimension attributes
  allocate(dim_names(numVariables,4))
  dim_names(:,:) = ""
  ! allocate the axis data and attribute types for the file
  !> \note The user should increase the sizes of the axis and data attributes to accommodate more axes if necessary.
  allocate(axis_data_CS%axis(7))
  allocate(axis_data_CS%data(7))
  ! axis registration procedure for the domain-decomposed case
  if (file_open_successDD) then
    do i=1,numVariables
      num_dims=0
      dim_lengths(:) = 0
      !> \note The time dimension is registered separately at the end of the procedure if reg_time = .true.
      !! so the t_grid argument in get_var_dimension_features is set to '1' (do nothing)
      if (present(G)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                          dim_lengths, num_dims, G=G)
      elseif(present(dG)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                          dim_lengths, num_dims, dG=dG)
      endif

      if(present(GV)) &
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, '1', dim_names(i,:), &
                                        dim_lengths, num_dims, GV=GV)
      !> \note num_dims will be 0 for scalar values
      if (num_dims .le. 0) cycle

      do j=1,num_dims
        ! register the variable axes to the file if they are not already registered
        if (dim_lengths(j) .gt. 0) then
          if (.not.(dimension_exists(fileObjDD, dim_names(i,j)))) then

            if (present(G)) then
              if (present(GV)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, GV=GV)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G)
              endif   
            elseif (present(dG)) then
              if (present(GV)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, GV=GV)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG)
              endif
            elseif (present(GV)) then
               call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV)
            endif
              call MOM_register_diagnostic_axis(fileObjDD, trim(dim_names(i,j)), dim_lengths(j))
          endif
          ! register the axis attributes and write the axis data to the file
          if (.not.(variable_exists(fileObjDD, trim(axis_data_CS%axis(j)%name)))) then
            if (associated(axis_data_CS%data(j)%p)) then

                call register_field(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                  "double", dimensions=(/trim(axis_data_CS%axis(j)%name)/))

                call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                               'long_name', axis_data_CS%axis(j)%longname)

                call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                               'units', trim(axis_data_CS%axis(j)%units))

                if (len_trim(axis_data_CS%axis(j)%positive)>1) &
                  call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                                   'positive', trim(axis_data_CS%axis(j)%positive))

                if (axis_data_CS%axis(j)%is_domain_decomposed) then
                  call get_global_io_domain_indices(fileObjDD, trim(axis_data_CS%axis(j)%name), is, ie)
                  call write_data(fileObjDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p(is:ie))
                else
                  call write_data(fileObjDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p)
                endif
            endif
          endif
        endif
      enddo
    enddo

    if (reg_time) then
     if (.not.(dimension_exists(fileObjDD,"Time"))) &
       call register_axis(fileObjDD, "Time", unlimited)
    endif
  else
    call MOM_error(FATAL, "MOM_io::create_file_filobj_dd: unable to open file "//trim(filename))
  endif

  deallocate(dim_names)
  deallocate(axis_data_CS%axis)
  deallocate(axis_data_CS%data)
  if (allocated(pelist)) deallocate(pelist)
  nullify(Domain)

end subroutine create_file_fileobj_dd


!> This function uses the fms_io function write_data to write a 1-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_1d_DD(filename, fieldname, data, mode, domain, var_desc, start_index, edge_lengths, time_level, &
                             time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:), intent(in) :: data !< The 1-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type),intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave the file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer, dimension(:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  integer :: num_dims, substring_index
  integer :: dim_unlim_size! size of the unlimited dimension
  integer, dimension(1) :: start, nwrite ! indices for first data value and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=1024) :: filename_temp
  character(len=48), dimension(2) :: dim_names !< variable dimension names (or name, in the 1-D case); 1 extra
                                               !! dimension in case appending along the time axis
  integer, dimension(2) :: dim_lengths !< variable dimension lengths (or length, in the 1-D case)

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_unlim_size=0
  dim_unlim_name=""
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! get the dimension names and lengths
  ! NOTE: the t_grid argument is set to '1' (do nothing) because the presence of a time dimension is user-specified
  ! and not assumed from the var_desc%t_grid value
  if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                      dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                      dim_lengths, num_dims, dG=dG)
  endif

  if (present(GV)) &
  call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                  dim_lengths, num_dims, GV=GV)
  ! define the start and edge_length arguments
  start(:) = 1
  nwrite(:) = dim_lengths(1)
  if (present(start_index)) then
    start(1) = max(1, start_index(1))
  endif

  if (present(edge_lengths)) then
    nwrite(1) = max(dim_lengths(1),edge_lengths(1))
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  if (.not.(check_if_open(fileobj_write_field_dd))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_1d_DD:mode argument must be write, overwrite, or append")
    ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    ! get the time_level index
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field_dd, trim(filename_temp), lowercase(trim(mode)), &
                                       domain%mpp_domain, is_restart=.false.)
    ! register the diagnostic axis associated with the variable
    call MOM_register_diagnostic_axis(fileobj_write_field_dd, trim(dim_names(1)), dim_lengths(1))
  endif
  ! register and write the time_level
  if (present(time_level)) then
    if (.not. (variable_exists(fileobj_write_field_dd, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field_dd, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field_dd, trim(dim_unlim_name), 'units', trim(t_units))
       call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time_level if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), &
                        corner=(/write_field_time_index/), edge_lengths=(/1/))
    endif
  endif
  ! register the field if it is not already in the file
  if (.not.(variable_exists(fileobj_write_field_dd, trim(fieldname)))) then
    call register_field(fileobj_write_field_dd, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), "checksum", checksum_char)
    endif
  endif
  ! write the variable
  if (present(time_level)) then
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field_dd)) call fms2_close_file(fileobj_write_field_dd)
    write_field_time_index = 0
  endif
  nullify(data_tmp)
end subroutine write_field_1d_DD

!> This function uses the fms_io function write_data to write a 2-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_2d_DD(filename, fieldname, data, mode, domain, var_desc, start_index, edge_lengths, time_level, &
                             time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:,:), intent(in) :: data !< The 2-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(2), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave the file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer :: data_tmp(:,:) => null()
  integer :: i, is, ie, js, je, j, ndims, num_dims, substring_index
  integer, allocatable, dimension(:) :: x_inds, y_inds
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer :: file_dim_length
  integer, dimension(2) :: start, nwrite ! indices for starting points and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(3) :: dim_names ! variable dimension names; 1 extra dimension in case appending
                                               ! along the time axis
  character(len=48), allocatable, dimension(:) :: file_dim_names
  integer, dimension(3) :: dim_lengths ! variable dimension lengths

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_lengths(:) = 0
  dim_names(:) = ""
  dim_unlim_size = 0
  dim_unlim_name = ""
  ndims = 2
  num_dims = 0
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! get the dimension names and lengths
  ! NOTE: the t_grid argument is set to '1' (do nothing) because the presence of a time dimension
  ! is user-specified rather than derived from the var_desc%t_grid value
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif
  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, GV=GV)
  ! set the start (start_index) and nwrite (edge_lengths) values
  start(:) = 1
  nwrite(:) = dim_lengths(1:ndims)

  if (present(start_index)) then
    do i=1,ndims
      start(i) = max(1,start_index(i))
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      nwrite(i) = max(dim_lengths(i),edge_lengths(i))
    enddo
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  if (.not.(check_if_open(fileobj_write_field_dd))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
       call MOM_error(FATAL,"MOM_io:write_2d_DD:mode argument must be write, overwrite, or append")
    ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    ! get the time_level index
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field_dd, trim(filename_temp), lowercase(trim(mode)), &
                                       domain%mpp_domain, is_restart=.false.)
  endif
  ! register the horizontal diagnostic axes associated with the variable
  do i=1,num_dims
    if (.not.(is_dimension_registered(fileobj_write_field_dd, trim(dim_names(i))))) &
      call MOM_register_diagnostic_axis(fileobj_write_field_dd, trim(dim_names(i)), dim_lengths(i))
  enddo
  ! register and write the time_level
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field_dd, dim_unlim_name)
    call get_dimension_size(fileobj_write_field_dd, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)

    if (.not. (variable_exists(fileobj_write_field_dd, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field_dd, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field_dd, trim(dim_unlim_name), 'units', trim(t_units))
       call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time_level if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), &
                        corner=(/write_field_time_index/), edge_lengths=(/1/))
    endif
  endif
  ! register the variable if it is not already in the file
  if (.not.(variable_exists(fileobj_write_field_dd, trim(fieldname)))) then
    call register_field(fileobj_write_field_dd, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), "checksum", checksum_char)
    endif
  endif
  ! write the variable
  if (present(time_level)) then
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field_dd)) call fms2_close_file(fileobj_write_field_dd)
    write_field_time_index=0
    if (allocated(file_dim_names)) deallocate(file_dim_names)
  endif
  nullify(data_tmp)
end subroutine write_field_2d_DD

!> This function uses the fms_io function write_data to write a 3-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_3d_DD(filename, fieldname, data, mode, domain, var_desc, start_index, edge_lengths, time_level, &
                             time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:,:,:), intent(in) :: data !< The 3-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type),  intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(3), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave the file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer, dimension(:,:,:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  integer :: i, is, ie, js, je, ndims, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer, dimension(3) :: start, nwrite ! indices for first data value and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names !< variable dimension names; 1 extra dimension in case appending
                                               !! along the time axis
  integer, dimension(4) :: dim_lengths !< variable dimension lengths

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_unlim_size = 0
  dim_unlim_name = ""
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! get the dimension names and lengths
  ! NOTE: the t_grid argument is set to '1' (do nothing) because the presence of a time dimension is user-specified
  ! and not assumed from the var_desc%t_grid value
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif

  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, GV=GV)
  ! set the start (start_index) and nwrite (edge_lengths) values
  ndims = 3
  start(:) = 1
  nwrite(:) = dim_lengths(1:3)
  if (present(start_index)) then
    do i=1,ndims
      start(i) = max(1,start_index(i))
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      nwrite(i) = max(dim_lengths(i), edge_lengths(i))
    enddo
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif
  ! open the file
  if (.not.(check_if_open(fileobj_write_field_dd))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_3d_DD:mode argument must be write, overwrite, or append")
    ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    ! get the time_level index
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field_dd, trim(filename_temp), lowercase(trim(mode)), &
                                       domain%mpp_domain, is_restart=.false.)
    ! register the horizontal and vertical diagnostic axes associated with the variable
    do i=1,ndims
      call MOM_register_diagnostic_axis(fileobj_write_field_dd, trim(dim_names(i)), dim_lengths(i))
    enddo
  endif
  ! register and write the time_level
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field_dd ,dim_unlim_name)
    call get_dimension_size(fileobj_write_field_dd, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)

    if (.not. (variable_exists(fileobj_write_field_dd, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field_dd, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field_dd, trim(dim_unlim_name), 'units', trim(t_units))
       call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time_level if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size ) &
        call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), &
                         corner=(/write_field_time_index/), edge_lengths=(/1/))
    endif
  endif
  ! register the field if it is not already in the file
  if (.not.(variable_exists(fileobj_write_field_dd, trim(fieldname)))) then
    call register_field(fileobj_write_field_dd, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), "checksum", checksum_char)
    endif
  endif
  ! write the data
  if (present(time_level)) then
    call get_dimension_size(fileobj_write_field_dd, trim(dim_unlim_name), dim_unlim_size)
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field_dd)) call fms2_close_file(fileobj_write_field_dd)
    write_field_time_index=0
  endif
  nullify(data_tmp)

end subroutine write_field_3d_DD

!> This function uses the fms_io function write_data to write a 4-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_4d_DD(filename, fieldname, data, mode, domain, var_desc, start_index, edge_lengths, time_level, &
                             time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:,:,:,:), intent(in) :: data !< The 4-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(4), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(4), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave the file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer, dimension(:,:,:,:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  real :: file_time ! most recent time currently written to file
  integer :: i, ndims, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer, dimension(4) :: start, nwrite ! indices for first data value and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names ! variable dimension names
  integer, dimension(4) :: dim_lengths ! variable dimension lengths

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  num_dims = 0
  dim_unlim_size = 0
  dim_unlim_name = ""
  dim_names(:) = ""
  dim_lengths(:) = 0
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! get the dimension names and lengths
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif

  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                    dim_lengths, num_dims, GV=GV)
  ! set the start (start_index) and nwrite (edge_lengths) values
  ndims = 4
  start(:) = 1
  nwrite(:) = dim_lengths(:)
  if (present(start_index)) then
    do i=1,ndims
      start(i) = max(1,start_index(i))
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      nwrite(i) = max(dim_lengths(i), edge_lengths(i))
    enddo
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif
  ! open the file
  if (.not.(check_if_open(fileobj_write_field_dd))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_4d_DD:mode argument must be write, overwrite, or append")
    ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    ! get the index of the corresponding time_level the first time the file is opened
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field_dd, trim(filename_temp), lowercase(trim(mode)), &
                                       domain%mpp_domain, is_restart=.false.)
    ! register the horizontal and vertical diagnostic axes associated with the variable
    do i=1,ndims
      call MOM_register_diagnostic_axis(fileobj_write_field_dd, trim(dim_names(i)), dim_lengths(i))
    enddo
  endif
  ! register the time dimension and write the time_level
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field_dd, dim_unlim_name)
    call get_dimension_size(fileobj_write_field_dd, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)

    if (.not. (variable_exists(fileobj_write_field_dd, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field_dd, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field_dd, trim(dim_unlim_name), 'units', trim(t_units))
       call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time_level if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field_dd, trim(dim_unlim_name), (/time_level/), &
                        corner=(/write_field_time_index/), edge_lengths=(/1/))
    endif
  endif
  ! register the variable if it is not already in the file
  if (.not.(variable_exists(fileobj_write_field_dd, trim(fieldname)))) then
    call register_field(fileobj_write_field_dd, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field_dd, trim(fieldname), "checksum", checksum_char)
    endif
  endif
  ! write the data
  if (present(time_level)) then
    call get_dimension_size(fileobj_write_field_dd, trim(dim_unlim_name), dim_unlim_size)
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field_dd, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field_dd)) call fms2_close_file(fileobj_write_field_dd)
    write_field_time_index=0
  endif
  nullify(data_tmp)
end subroutine write_field_4d_DD

!> This routine uses the fms_io function write_data to write a scalar variable named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_scalar(filename, fieldname, data, mode, var_desc, time_level, time_units, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, intent(in) :: data !< The 1-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave the file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=48), dimension(1) :: dim_names ! variable dimension names
  integer :: i, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  real, allocatable, dimension(:) :: file_times
  integer, dimension(1) :: dim_lengths ! variable dimension lengths
  integer, allocatable, dimension(:) :: pelist ! list of pes associated with the netCDF file

  dim_unlim_size = 0
  dim_unlim_name= ""
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  if (.not.(check_if_open(fileobj_write_field))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_scaler:mode argument must be write, overwrite, or append")
    ! get the index of the corresponding time_level the first time the file is opened 
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    if (.not.(allocated(pelist))) then
      allocate(pelist(mpp_npes()))
      pelist(:) = 0
      do i=1,size(pelist)
        pelist(i) = i-1
      enddo
    endif
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field, trim(filename_temp), trim(mode), is_restart=.false., &
                                       pelist=pelist)
  endif
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field, dim_unlim_name)
    call get_dimension_size(fileobj_write_field, trim(dim_unlim_name), dim_unlim_size)
    ! write the time value if it is not already written to the file
    if (.not.(variable_exists(fileobj_write_field, trim(dim_unlim_name)))) then
      ! set the time units
      t_units = ""
      if (present(time_units)) then
        t_units = get_time_units(time_units)
      else
        t_units = "days"
      endif

      call register_field(fileobj_write_field, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
      call register_variable_attribute(fileobj_write_field, trim(dim_unlim_name), 'units', trim(t_units))
      call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/))
    else
      ! write the next time value if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/), &
                        edge_lengths=(/1/))
    endif
  endif
  ! register the variable
  if (.not.(variable_exists(fileobj_write_field, trim(fieldname)))) then
    if (present(time_level)) then
      call register_field(fileobj_write_field, trim(fieldname), "double", dimensions=(/trim(dim_unlim_name)/))
    else
      call register_field(fileobj_write_field, trim(fieldname), "double")
    endif
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'long_name', trim(var_desc%longname))
  endif
  ! write the data
  if (present(time_level)) then
    call write_data(fileobj_write_field, trim(fieldname), data, unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field, trim(fieldname), data)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field)) call fms2_close_file(fileobj_write_field)
    if (allocated(pelist)) deallocate(pelist)
    write_field_time_index=0
  endif
end subroutine write_scalar

!> This function uses the fms_io function write_data to write a 1-D non-domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_1d_noDD(filename, fieldname, data, mode, var_desc, start_index, edge_lengths, time_level, &
                               time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:), intent(in) :: data !< The 1-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  real, pointer, dimension(:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  integer :: i, ndims, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer, dimension(1) :: start, nwrite ! indices for first data value and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(2) :: dim_names ! variable dimension names (up to 2 if appended at time level)
  integer, dimension(2) :: dim_lengths ! variable dimension lengths
  integer, allocatable, dimension(:) :: pelist ! list of pes associated with the netCDF file
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_unlim_size = 0
  dim_unlim_name= "Time"
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! get the dimension names and lengths
  ! NOTE: the t_grid argument is set to '1' (do nothing) because the presence of a time dimension is user-specified
  ! and not assumed from the var_desc%t_grid value.
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif

  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, GV=GV)
  ! set the start (start_index) and nwrite (edge_lengths) values
  start(:) = 1
  nwrite(:) = dim_lengths(1)
  if (present(start_index)) then
    start(1) = max(1,start_index(1))
  endif

  if (present(edge_lengths)) then
    nwrite(1) = max(dim_lengths(1),edge_lengths(1))
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  if (.not.(check_if_open(fileobj_write_field))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
       (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_1d_noDD:mode argument must be write, overwrite, or append")
    ! get the index of the corresponding time_level the first time the file is opened 
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    if (.not.(allocated(pelist))) then 
      allocate(pelist(mpp_npes()))
      pelist(:) = 0
      do i=1,size(pelist)
        pelist(i) = i-1
      enddo
    endif
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field, trim(filename_temp), lowercase(trim(mode)), &
                                       is_restart=.false., pelist=pelist)
  endif
  ! write the data, and the time value if it is not already written to the file
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field,dim_unlim_name)
    call get_dimension_size(fileobj_write_field, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)

    if (.not. (variable_exists(fileobj_write_field, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field, trim(dim_unlim_name), 'units', trim(t_units))
      call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time value if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/), &
                        edge_lengths=(/1/))
    endif
  endif
  ! register the field if it is not already in the file
  if (.not.(variable_exists(fileobj_write_field, trim(fieldname)))) then
    call register_field(fileobj_write_field, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
    ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field, trim(fieldname), "checksum", checksum_char)
    endif
  endif
  ! write the variable to the file
  if (present(time_level)) then
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field)) call fms2_close_file(fileobj_write_field)
    if (allocated(pelist)) deallocate(pelist)
    write_field_time_index = 0
  endif
  nullify(data_tmp)

end subroutine write_field_1d_noDD

!> This function uses the fms_io function write_data to write a scalar variable named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_2d_noDD(filename, fieldname, data, mode, var_desc, start_index, edge_lengths, time_level, &
                               time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:,:), intent(in) :: data !< The 2-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(2), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success ! .true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer, dimension (:,:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  integer :: i, ndims, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer, dimension(2) :: start, nwrite ! indices for starting points and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(3) :: dim_names ! variable dimension names
  integer, dimension(3) :: dim_lengths ! variable dimension lengths
  integer, allocatable, dimension(:) :: pelist ! list of pes associated with the netCDF file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_unlim_size = 0
  dim_unlim_name = ""
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  ! get the dimension names and lengths
  ! NOTE: the t_grid argument is set to '1' (do nothing) because the presence of a time dimension is user-specified
  ! and not assumed from the var_desc%t_grid value
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif

  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                      dim_lengths, num_dims, GV=GV)

  ! set the start (start_index) and nwrite (edge_lengths) values
  ndims=2
  start(:) = 1
  nwrite(:) = dim_lengths(1:2)
  if (present(start_index)) then
    do i=1,ndims
      start(i) = max(1,start_index(i))
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      nwrite(i) = max(dim_lengths(i),edge_lengths(i))
    enddo
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  if (.not.(check_if_open(fileobj_write_field))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_2d_noDD:mode argument must be write, overwrite, or append")
    ! get the time_level index
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    if(.not.(allocated(pelist))) then
      allocate(pelist(mpp_npes()))
      pelist(:) = 0
      do i=1,size(pelist)
        pelist(i) = i-1
      enddo
    endif
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field, trim(filename_temp), lowercase(trim(mode)), &
                                       is_restart=.false., pelist=pelist)
  endif

  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field,dim_unlim_name)
    call get_dimension_size(fileobj_write_field, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)

    if (.not. (variable_exists(fileobj_write_field, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field, trim(dim_unlim_name), 'units', trim(t_units))
      call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time value if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/), &
                        edge_lengths=(/1/))
    endif
  endif

   ! register the variable to the file
  if (.not.(variable_exists(fileobj_write_field, trim(fieldname)))) then
    call register_field(fileobj_write_field, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field, trim(fieldname), "checksum", checksum_char)
    endif
  endif
  ! write the variable to the file
  if (present(time_level)) then
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field)) call fms2_close_file(fileobj_write_field)
    if (allocated(pelist)) deallocate(pelist)
    write_field_time_index=0
  endif
  nullify(data_tmp)
end subroutine write_field_2d_noDD

!> This function uses the fms_io function write_data to write a 3-D non-domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_3d_noDD(filename, fieldname, data, mode, var_desc, start_index, edge_lengths, time_level, &
                               time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:,:,:), intent(in) :: data !< The 3-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(4), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(4), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer, dimension(:,:,:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  integer :: i, ndims, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer, dimension(3) :: start, nwrite ! indices for first data value and number of values to write
  character(len=20) :: t_units ! time_units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names ! variable dimension names
  integer, dimension(4) :: dim_lengths ! variable dimension lengths
  integer, allocatable, dimension(:) :: pelist ! list of pes associated with the netCDF file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_unlim_size = 0
  dim_unlim_name = ""
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! get the dimension names and lengths
  ! NOTE: the t_grid argument is set to '1' (do nothing) because the presence of a time dimension is user-specified
  ! and not assumed from the var_desc%t_grid value
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif

  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, '1', dim_names, &
                                    dim_lengths, num_dims, GV=GV)
  ! set the start (start_index) and nwrite (edge_lengths) values
  ndims = 3
  start(:) = 1
  nwrite(:) = dim_lengths(1:3)
  if (present(start_index)) then
    do i=1,ndims
      start(i) = max(1,start_index(i))
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      nwrite(i) = max(dim_lengths(i), edge_lengths(i))
    enddo
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif
  ! open the file
  if (.not.(check_if_open(fileobj_write_field))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_3d_noDD:mode argument must be write, overwrite, or append")
    ! get the time_level index
    if (present(time_level))  write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    if (.not.(allocated(pelist))) then
      allocate(pelist(mpp_npes()))
      pelist(:) = 0
      do i=1,size(pelist)
        pelist(i) = i-1
      enddo
    endif
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field, trim(filename_temp), lowercase(trim(mode)), &
                                       is_restart=.false., pelist=pelist)
  endif
  ! register and write the time_level
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field,dim_unlim_name)
    call get_dimension_size(fileobj_write_field, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)

    if (.not. (variable_exists(fileobj_write_field, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field, trim(dim_unlim_name), 'units', trim(t_units))
       call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      ! write the time_level if it is larger than the most recent file time
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/), &
                        edge_lengths=(/1/))
    endif
  endif
  ! register the field if it is not already in the file
  if (.not.(variable_exists(fileobj_write_field, trim(fieldname)))) then
    call register_field(fileobj_write_field, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileobj_write_field, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  if (present(time_level)) then
    call get_dimension_size(fileobj_write_field, trim(dim_unlim_name), dim_unlim_size)
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field)) call fms2_close_file(fileobj_write_field)
    if (allocated(pelist)) deallocate(pelist)
    write_field_time_index=0
  endif
  nullify(data_tmp)
end subroutine write_field_3d_noDD

!> This function uses the fms_io function write_data to write a 4-D non-domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_4d_noDD(filename, fieldname, data, mode, var_desc, start_index, edge_lengths, time_level, &
                               time_units, scale, checksums, G, dG, GV, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, target, dimension(:,:,:,:), intent(in) :: data !< The 4-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(4), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(4), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: time_units !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  real, pointer, dimension(:,:,:,:) :: data_tmp => null() ! enables data to be passed to functions as intent(inout)
  integer :: i, ndims, num_dims, substring_index
  integer :: dim_unlim_size ! size of the unlimited dimension
  integer, dimension(4) :: start, nwrite ! indices for first data value and number of values to write
  character(len=20) :: t_units ! time units
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names ! variable dimension names
  integer, dimension(4) :: dim_lengths ! variable dimension lengths
  integer, allocatable, dimension(:) :: pelist ! list of pes associated with the netCDF file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_unlim_size = 0
  dim_unlim_name = ""
  dim_names(:) = ""
  dim_lengths(:) = 0
  ndims = 4
  num_dims = 0
  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  ! get the dimension names and lengths
  if (present(G)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                    dim_lengths, num_dims, G=G)
  elseif(present(dG)) then
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                    dim_lengths, num_dims, dG=dG)
  endif
  if (present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                    dim_lengths, num_dims, GV=GV)
  ! set the start (start_index) and nwrite (edge_lengths) values
  start(:) = 1
  nwrite(:) = dim_lengths(:)
  if (present(start_index)) then
    do i=1,ndims
      start(i) = max(1, start_index(i))
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      nwrite(i) = max(dim_lengths(i), edge_lengths(i))
    enddo
  endif

  data_tmp => data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif
  ! open the file
  if (.not.(check_if_open(fileobj_write_field))) then
    if ((lowercase(trim(mode)) .ne. "write") .and. (lowercase(trim(mode)) .ne. "append") .and. &
        (lowercase(trim(mode)) .ne. "overwrite")) &
      call MOM_error(FATAL,"MOM_io:write_4d_noDD:mode argument must be write, overwrite, or append")
    ! get the time_level index
    if (present(time_level)) write_field_time_index = get_time_index(trim(filename_temp), time_level)
    ! get the pes associated with the file.
    !>\note this is required so that only pe(1) is identified as the root pe to create the file
    !! Otherwise, multiple pes may try to open the file in write (NC_NOCLOBBER) mode, leading to failure
    if (.not.(allocated(pelist))) then
      allocate(pelist(mpp_npes()))
      pelist(:) = 0
      do i=1,size(pelist)
        pelist(i) = i-1
      enddo
    endif
    ! open the file in write or append mode
    file_open_success = fms2_open_file(fileobj_write_field, trim(filename_temp), lowercase(trim(mode)), &
                                       is_restart=.false., pelist=pelist)
  endif
  ! register and write the time_level
  if (present(time_level)) then
    call get_unlimited_dimension_name(fileobj_write_field,dim_unlim_name)
    call get_dimension_size(fileobj_write_field, trim(dim_unlim_name), dim_unlim_size)
    num_dims=num_dims+1
    dim_names(num_dims) = trim(dim_unlim_name)
    ! write the time value if it is not already written to the file
    if (.not. (variable_exists(fileobj_write_field, trim(dim_unlim_name)))) then
       ! set the time units
       t_units = ""
       if (present(time_units)) then
         t_units = get_time_units(time_units)
       else
         t_units = "days"
       endif

       call register_field(fileobj_write_field, trim(dim_unlim_name), "double", dimensions=(/trim(dim_unlim_name)/))
       call register_variable_attribute(fileobj_write_field, trim(dim_unlim_name), 'units', trim(t_units))
       call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/))
    else
      if (write_field_time_index .gt. dim_unlim_size) &
        call write_data(fileobj_write_field, trim(dim_unlim_name), (/time_level/), corner=(/write_field_time_index/), &
                        edge_lengths=(/1/))
    endif
  endif
  ! register the variable  
  if (.not.(variable_exists(fileobj_write_field, trim(fieldname)))) then
     call register_field(fileobj_write_field, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
     call register_variable_attribute(fileobj_write_field, trim(fieldname), 'units', trim(var_desc%units))
     call register_variable_attribute(fileobj_write_field, trim(fieldname), 'long_name', trim(var_desc%longname))
     ! write the checksum attribute
     if (present(checksums)) then
       ! convert the checksum to a string
       checksum_char = ""
       checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileobj_write_field, trim(fieldname), "checksum", checksum_char)
     endif
  endif
  ! write the variable to the file
  if (present(time_level)) then
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                    unlim_dim_level=write_field_time_index)
  else
    call write_data(fileobj_write_field, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_write_field)) call fms2_close_file(fileobj_write_field)
    deallocate(pelist)
    write_field_time_index=0
  endif
  nullify(data_tmp)
end subroutine write_field_4d_noDD

!> This routine calls the fms_io read_data subroutine to read 1-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_1d_DD(filename, fieldname, data, domain, start_index, edge_lengths, timelevel, scale, &
                               x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:), intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                             !! variable size
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, num_var_dims, dim_unlim_size
  integer, dimension(1) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos ! x and y domain positions

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    ! define the io domain for 1-pe jobs because it is required to read domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_1d_DD: "//trim(fieldname)//" not found in"//&
                                            trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)
  ! set the start and nread values that will be passed as the read_data start_index and edge_lengths arguments
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)

  start(1)=1
  if (present(timelevel)) then
    if (is_dimension_unlimited(fileobj_read_dd, dim_names(1))) start(1) = timelevel
  elseif (present(start_index)) then
    start(1) = start_index(1)
  endif

  if (present(edge_lengths)) then
    nread(1) = edge_lengths(1)
  else
    call get_dimension_size(fileobj_read_dd, trim(dim_names(1)), nread(1))
  endif
  ! read the data
  dim_unlim_size = 0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
        exit
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_io::MOM_read_data_1d_DD: time level specified, but the variable "//&
                      trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
     if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
     file_var_meta_DD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_read_data_1d_DD

!> This routine calls the fms_io read_data subroutine to read 2-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_2d_DD(filename, fieldname, data, domain, start_index, edge_lengths, timelevel, scale, &
                               x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename  !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:), intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(2), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims, first(2), last(2)
  integer :: start(2), nread(2) ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg
  type(domain2D), pointer :: io_domain => NULL()

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    ! define the io domain for 1-pe jobs because it is required to read domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_2d_DD: "//trim(fieldname)//" not found in "//&
                                            trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)

  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read_dd, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)
  ! get the IO domain
  io_domain => mpp_get_io_domain(domain%mpp_domain)
  ! Get the global indicies
  call mpp_get_global_domain(io_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, position=pos)
  ! Get the compute indicies
  call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=pos)
  last(1) = iec - isg + 1 ! get array indices for the axis data
  last(2) = jec - jsg + 1
  first(1) = isc - isg + 1
  first(2) = jsc - jsg + 1

  start(:) = 1
  if (present(start_index)) then
    start = start_index
  else
    start(:) = first(:)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    !do i=1,num_var_dims
    !  call get_dimension_size(fileobj_read_dd, trim(dim_names(i)), nread(i))
    !enddo
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
    if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
    file_var_meta_DD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)
  if (associated(io_domain)) nullify(io_domain)
end subroutine MOM_read_data_2d_DD

!> This routine calls the fms_io read_data subroutine to read 3-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_3d_DD(filename, fieldname, data, domain, start_index, edge_lengths, timelevel, scale, &
                               x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename  !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(3), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                             !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! if .true., the variable was found in the netCDF file
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(3) :: start, nread, first, last ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg
  type(domain2D), pointer :: io_domain => NULL()

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    ! define the io domain for 1-pe jobs because it is required to read domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_3d_DD: "//trim(fieldname)//" not found in"//&
                                            trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)
  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read_dd, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)
  ! get the IO domain
  io_domain => mpp_get_io_domain(domain%mpp_domain)
  ! Get the global indicies
  call mpp_get_global_domain(io_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, position=pos)
  ! Get the compute indicies
  call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=pos)
  last(1) = iec - isg + 1 ! get array indices for the axis data
  last(2) = jec - jsg + 1
  first(1) = isc - isg + 1
  first(2) = jsc - jsg + 1

  call mpp_get_global_domain(domain%mpp_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg,    position=pos) ! Get the global indicies
  call mpp_get_compute_domain(domain%mpp_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=pos) ! Get the compute indicies
  last(1) = iec - isg + 1 ! get array indices for the axis data
  last(2) = jec - jsg + 1
  first(1) = isc - isg + 1
  first(2) = jsc - jsg + 1

  start(:) = 1
  if (present(start_index)) then
    start = start_index
  else
    start(1:2) = first(1:2)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread(1) = last(1) - first(1) + 1
    nread(2) = last(2) - first(2) + 1
    call get_dimension_size(fileobj_read_dd, trim(dim_names(3)), nread(3))
    !nread = shape(data)
  endif
 ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_io::MOM_read_data_3d_DD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
     if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
    file_var_meta_DD%nvars = 0
  endif

  if (allocated(dim_names)) deallocate(dim_names)
  if (associated(io_domain)) nullify(io_domain)
end subroutine MOM_read_data_3d_DD

!> This routine calls the fms_io read_data subroutine to read 4-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_4d_DD(filename, fieldname, data, domain, start_index, edge_lengths, timelevel, scale, &
                               x_position, y_position, leave_file_open)
  character(len=*),       intent(in) :: filename  !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data !< The 4-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(4),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(4) :: start, nread, first, last ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg
  type(domain2D), pointer :: io_domain => NULL()

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    ! define the io domain for 1-pe jobs because it is required to read domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_4d_DD: "//trim(fieldname)//" not found in"//&
                                            trim(filename))
   ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)
  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read_dd, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)
  ! get the IO domain
  io_domain => mpp_get_io_domain(domain%mpp_domain)
  ! Get the global indicies
  call mpp_get_global_domain(domain%mpp_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, position=pos)
  ! Get the compute indicies
  call mpp_get_compute_domain(domain%mpp_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=pos)
  last(1) = iec - isg + 1 ! get array indices for the axis data
  last(2) = jec - jsg + 1
  first(1) = isc - isg + 1
  first(2) = jsc - jsg + 1

  call mpp_get_global_domain(domain%mpp_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg,    position=pos) ! Get the global indicies
  call mpp_get_compute_domain(domain%mpp_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=pos) ! Get the compute indicies
  last(1) = iec - isg + 1 ! get array indices for the axis data
  last(2) = jec - jsg + 1
  first(1) = isc - isg + 1
  first(2) = jsc - jsg + 1

  start(:) = 1
  if (present(start_index)) then
    start(:) = start_index(:)
  else
    start(1:2) = first(1:2)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread(1) = last(1) - first(1) + 1
    nread(2) = last(2) - first(2) + 1
    do i=3,num_var_dims
      call get_dimension_size(fileobj_read_dd, trim(dim_names(i)), nread(i))
    enddo
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1, num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
      endif
      if (i .eq. 4) then
        nread(i) = 1
        start(i) = timelevel
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_io::MOM_read_data_4d_DD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
     if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
    file_var_meta_DD%nvars = 0
  endif
  if (associated(io_domain)) nullify(io_domain)
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_read_data_4d_DD

!!> This routine calls the fms_io read_data subroutine to read a scalar (0-D) field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_scalar(filename, fieldname, data, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, intent(inout) :: data !< data buffer to pass to read_data
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  integer :: i
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_scalar: "//trim(fieldname)// &
                                            " not found in"//trim(filename))
  ! read the data
  call read_data(fileobj_read, trim(fieldname), data)
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
end subroutine MOM_read_data_scalar

!> This routine calls the fms_io read_data subroutine to read 1-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_1d_noDD(filename, fieldname, data, start_index, edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is
                                                             !! the variable size
  integer,      optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  integer :: i, num_var_dims, dim_unlim_size
  integer, dimension(1) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable:: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_1d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))

  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)

  ! set the start and nread values that will be passed as the read_data start_index and edge_lengths arguments
  start(1)=1
  if (present(timelevel)) then
    if (is_dimension_unlimited(fileobj_read, dim_names(1))) start(1) = timelevel
  elseif (present(start_index)) then
    start(1) = start_index(1)
  endif

  if (present(edge_lengths)) then
    nread(1) = edge_lengths(1)
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size = 0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
        exit
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_io::MOM_read_data_1d_noDD: time level specified, but the variable "//&
                      trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_read_data_1d_noDD

!> This routine calls the fms_io read_data subroutine to read 2-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_2d_noDD(filename, fieldname, data, start_index, edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  integer, dimension(2),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(2) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) &
    call MOM_error(FATAL, "MOM_io:MOM_read_data_2d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
    allocate(dim_names(num_var_dims))
    call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
    dim_names(:) = ""
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .LE. 0) then
      call MOM_error(WARNING, "MOM_io::MOM_read_data_2d_noDD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if(allocated(dim_names)) deallocate(dim_names)

end subroutine MOM_read_data_2d_noDD

!> This routine calls the fms_io read_data subroutine to read 3-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_3d_noDD(filename, fieldname, data, start_index, edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  integer, dimension(3), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(3) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_3d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))
  ! get the variable dimensions
  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .LE. 0) then
      call MOM_error(WARNING, "MOM_io::MOM_read_data_3d_noDD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_read_data_3d_noDD

!> This routine calls the fms_io read_data subroutine to read 4-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_4d_noDD(filename, fieldname, data, start_index, edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:),   intent(inout) :: data !< The 4-dimensional array to pass to read_data
  integer, dimension(4),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after read_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(4) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_4d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))
  ! get the variable dimensions
  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1, num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
      if (i .eq. 4) then
        nread(i) = 1
        start(i) = timelevel
      endif
    enddo
    if (dim_unlim_size .LE. 0) then
      call MOM_error(WARNING, "MOM_io::MOM_read_data_4d_noDD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
 ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_read_data_4d_noDD

!> This routine calls the fms_io read_data subroutine to read 2-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call. The supergrid variable axis lengths are determined from compute domain lengths, and 
!! the domain indices are computed from the difference between the global and compute domain indices
subroutine MOM_read_data_2d_supergrid(filename, fieldname, data, domain, is_supergrid, start_index, edge_lengths, &
                              timelevel, scale, x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename  !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:), intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  logical, intent(in) :: is_supergrid !< flag indicating whether to use supergrid
  integer, dimension(2), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims, first(2), last(2)
  integer :: start(2), nread(2) ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg
  type(domain2D), pointer :: io_domain => NULL()

  if (.not.(is_supergrid)) call MOM_read_data(filename, fieldname, data, domain, start_index, edge_lengths, &
                                              timelevel, scale, x_position, y_position, leave_file_open)
  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    ! define the io domain for 1-pe jobs because it is required to read domain-decomposed files
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_2d_supergrid: "//&
      trim(fieldname)//" not found in "//trim(filename))

  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
  ! get the IO domain
  io_domain => mpp_get_io_domain(domain%mpp_domain)
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read, trim(variable_to_read), io_domain, xPosition=xpos, yPosition=ypos) 
  ! Get the global indicies
  call mpp_get_global_domain(io_domain, xbegin=isg, xend=ieg, ybegin=jsg, yend=jeg, position=pos)
  ! Get the compute indicies
  call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec, position=pos)
  ! get array indices for the axis data
  last(1) = iec - isg + 1
  last(2) = jec - jsg + 1
  first(1) = isc - isg + 1
  first(2) = jsc - jsg + 1

  start(:) = 1
  if (present(start_index)) then
    start = start_index
  else
    start(:) = first(:)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread(1) = last(1) - first(1) + 1
    nread(2) = last(2) - first(2) + 1
    !do i=1,num_var_dims
    !  call get_dimension_size(fileobj_read, trim(dim_names(i)), nread(i))
    !enddo
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)
  if (associated(io_domain)) nullify(io_domain)
end subroutine MOM_read_data_2d_supergrid


!> This is a kluge interface to obtain the correct x and y values used to define the diagnostic axes 
!! lath, latq, lonh, and lonq in MOM_grid_intialize. This routine allocates a buffer of the size of the entire field in
!! the input file. If reading in 'y' (latitude), a mask is created to read in a subset of supergrid points using the
!! indexing conventions from set_grid_metrics_from_mosaic, and is applied to the data buffer.
!! Next, a search is performed on the masked data for the column (longitude) that contains the maximum number of non-NAN
!! latitude points, and the resulting index is used to define output data array that sets the "lath" and "latq" values.
!! Finally, the indexed data is multipled by the optional scale argument if needed.
subroutine MOM_read_data_2d_noDD_diag_axes(filename, fieldname, data, define_diagnostic_axes, G, start_index, &
                                           edge_lengths, timelevel, scale, grid_type, leave_file_open)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  logical, intent(in) :: define_diagnostic_axes !< if .true., read in the full data array, search for the
                                                   !! full ranges of x/lon and y/lat values needed to define
                                                   !! the diagnostic axes
  type(dyn_horgrid_type), optional, intent(in) :: G !< The dynamic horizontal grid type; required if
                                                    !! define_diagnostic_axes=.true.
  integer, dimension(2),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edge_lengths !< number of data values to read in.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  character(len=*), optional, intent(in) :: grid_type !< grid type; either 'b' or 't'
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after read_data is called; default is .true.
  logical, dimension(:), allocatable :: yuse ! Mask for geoLatT and geoLatB values
  integer :: i, j, dim_unlim_index, substring_index, xidx, yidx1, yidx2
  integer, dimension(2) :: start, nread, dimSizes ! indices for first data value and number of values to read,
                                                  ! variable dimension sizes
  character(len=40), dimension(2) :: dim_names ! variable dimension names
  character(len=40) :: units ! variable units
  character(len=1) :: gtype ! grid type
  real :: ymax1, ymax2 ! max values from latitude array
  real, allocatable, dimension(:,:) :: tmpGlbl ! temporary array to hold data
  
  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  dim_names(:) = ""
  gtype=""
  units=""
  substring_index = 0
  ! open the file
  if (.not.(check_if_open(fileobj_read))) &
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)

  call get_variable_dimension_names(fileobj_read, trim(fieldname), dim_names)

  do i=1,2
    call get_dimension_size(fileobj_read, trim(dim_names(i)), nread(i))
  enddo

  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) nread = edge_lengths

  if (present(timelevel)) then
    dim_unlim_index=0
    do i=1,2
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        dim_unlim_index=i
        start(i)=timelevel
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD_diag_axes: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
  !> \note: the following indexing procedure searches for global x-dimension indices that correspond to the
  !! full ranges of y (geoLatT and geoLatB) values needed for the diagnostic indices (lath and latq).
  !! Defining the gridLatT and gridLatB values using the previous values in tmpGlbl(1,:) did not result in the correct
  !! diagnostic index values due to different indexing conventions for non-domain-decomposed IO in the new and
  !! old procedures.
  if (define_diagnostic_axes) then
    allocate(tmpGlbl(nread(1),nread(2)))
    !tmpGlbl(:,:) = -999.0
    ! read the data into the temporary array
    call read_data(fileobj_read, trim(fieldname), tmpGlbl, corner=start, edge_lengths=nread)

    call get_variable_units(fileobj_read, fieldname, units)
    substring_index = index(lowercase(trim(units)), "north")

    if ((substring_index .gt. 0) .or. (index(lowercase(trim(fieldname)), 'y') .gt. 0)) then
      if (.not.(present(grid_type))) call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD_diag_axes: grid_type"// &
         " argument must be defined if define_diag_axes=.true. and reading in a y-axis/latitude variable")

      ! create a mask for the T-grid latitude values in the tmpGlbl array
      allocate(yuse(nread(2)))
      yuse(:) = .FALSE.

      gtype = lowercase(trim(grid_type))
      select case(gtype)
        case ('t')
          do j=G%jsg,G%jeg
            yuse(2*(j-G%jsg)+2) = .TRUE.
          enddo
        case ('b')
          do j=G%jsg-1,G%jeg
            yuse(2*(j-G%jsg)+3) = .TRUE.
          enddo
        case default
          call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD_diag_axes: grid_type must be either t or b")
      end select
      ! get the index for the x-dimension with the maximum T-grid latitude (geoLatT) value in the tmpGlbl array
      xidx = 0
      yidx1 = 0
      ymax1 = -999.0
      yidx2 = 0
      ymax2 = -999.0

      do i=1,nread(1)
        ! find index of the maximum T-grid latitude value for the ith x (longitude) dimension
        yidx1 = MAXLOC(tmpGlbl(i,:), 1, MASK=yuse)
        ymax1 = tmpGLbl(i,yidx1)
        ! if the new max is greater than the current max, set the current max to the new max value
        if ( MAX(ymax1,ymax2) .EQ. ymax1) then
          yidx2 = yidx1
          ymax2 = ymax1
          xidx = i
        endif
      enddo

      if (xidx .LT. 1) call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_nodDD_diag_axes: xidx is less than 1")
      if ((ymax2 .LT. -90.0) .or. (ymax2 .GT. 90.0)) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_nodDD_diag_axes: ymax2 is not between -90.0 and 90.0")

      data(1,1:nread(2)) = tmpGlbl(xidx,:)

      deallocate(yuse)
    else
      substring_index = index(lowercase(trim(units)), "east")
      if ((substring_index .gt. 0) .or. (index(lowercase(trim(fieldname)), 'x') .gt. 0)) then
        ! all latitude indices contain the full range of x/longitude values required for the diagnostic axes
        data(1:nread(1),:) = tmpGlbl(:,1:size(data,2))
      else
        call MOM_error(FATAL, "MOM_io:MOM_read_data_2d_noDD_diag_axes: axis units "//trim(units)// " are invalid."// &
                       "lowercase(units) must contain 'north' or 'east' to indicate the axis orientation.")
      endif
    endif

    deallocate(tmpGlbl)
  else ! just read the data into the user-specified data buffer
    call read_data(fileobj_read, trim(fieldname), data, corner=start, edge_Lengths=nread)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
  endif

end subroutine MOM_read_data_2d_noDD_diag_axes

!> This routine uses the fms_io read_data interface to read a pair of distributed
!! 2-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_2d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scale, leave_file_open)
  character(len=*),       intent(in)    :: filename !< name of the netcdf file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  integer :: is, ie, js, je, i, ndims, dim_unlim_index
  integer :: u_pos, v_pos
  integer, allocatable :: dim_sizes_u(:), dim_sizes_v(:)
  character(len=32), allocatable :: dim_names_u(:), dim_names_v(:), units_u(:), units_v(:)
  character(len=1) :: x_or_y ! orientation of cartesian coordinate axis
  logical :: is_valid
  logical :: file_open_success ! .true. if open file is successful
  logical :: close_the_file ! indicates whether to close the file after MOM_read_vector is called; default is .true.

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) &
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
  if (.not. file_open_success) call MOM_error(FATAL, "MOM_read_vector_2d: netcdf file "//trim(filename)//" not opened.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE .or. stagger == BGRID_NE ) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  ndims = get_variable_num_dimensions(fileobj_read_dd, u_fieldname)
  allocate(dim_sizes_u(ndims))
  allocate(dim_sizes_v(ndims))
  allocate(dim_names_u(ndims))
  allocate(dim_names_v(ndims))
  allocate(units_u(ndims))
  allocate(units_v(ndims))

  units_u(:) = ""
  units_v(:) = ""
  dim_names_u(:) = ""
  dim_names_v(:) = ""
  dim_sizes_u(:) = 0
  dim_sizes_v(:) = 0

  call get_variable_size(fileobj_read_dd, u_fieldname, dim_sizes_u)
  call get_variable_size(fileobj_read_dd, v_fieldname, dim_sizes_v)
  call get_variable_dimension_names(fileobj_read_dd, u_fieldname, dim_names_u)
  call get_variable_dimension_names(fileobj_read_dd, v_fieldname, dim_names_v)

  do i=1,ndims
    ! register the u axes
    if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_u(i)))) then
      call get_variable_units(fileobj_read_dd, dim_names_u(i), units_u(i))
      call validate_lat_lon_units(units_u(i), x_or_y, is_valid) 
      if (is_valid) then
        call register_axis(fileobj_read_dd, dim_names_u(i), x_or_y, domain_position=u_pos)
      else
        call register_axis(fileobj_read_dd, dim_names_u(i), dim_sizes_u(i))
      endif
    endif
    ! Register the v axes if they differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then
      if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_v(i)))) then
        call get_variable_units(fileobj_read_dd, dim_names_v(i), units_v(i))
        call validate_lat_lon_units(units_v(i), x_or_y, is_valid) 
        if (is_valid) then
          call register_axis(fileobj_read_dd, dim_names_v(i), x_or_y, domain_position=v_pos)
        else
          call register_axis(fileobj_read_dd, dim_names_v(i), dim_sizes_v(i))
        endif
      endif
    endif
  enddo
  ! read the data
  dim_unlim_index = 0
  if (present(timelevel)) then
    do i=1,ndims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names_u(i))) then
        dim_unlim_index = i
        exit
      endif
    enddo
    if (dim_unlim_index .gt. 0) then
      call read_data(fileobj_read_dd, u_fieldname,u_data, unlim_dim_level=timelevel)
      call read_data(fileobj_read_dd, v_fieldname, v_data, unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read_dd, u_fieldname, u_data)
      call read_data(fileobj_read_dd, v_fieldname, v_data) 
    endif
  else
    call read_data(fileobj_read_dd, u_fieldname, u_data)
    call read_data(fileobj_read_dd, v_fieldname, v_data)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(u_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(u_data,2), js, je)
    u_data(is:ie,js:je) = scale*u_data(is:ie,js:je)
    call get_simple_array_i_ind(MOM_Domain, size(v_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(v_data,2), js, je)
    v_data(is:ie,js:je) = scale*v_data(is:ie,js:je)
  endif ; endif
  if (allocated(dim_names_u)) deallocate(dim_names_u)
  if (allocated(dim_names_v)) deallocate(dim_names_v)
  if (allocated(dim_sizes_u)) deallocate(dim_sizes_u)
  if (allocated(dim_sizes_v)) deallocate(dim_sizes_v)
  if (allocated(units_u)) deallocate(units_u)
  if (allocated(units_v)) deallocate(units_v)
end subroutine MOM_read_vector_2d

!> This routine uses the fms_io read_data interface to read a pair of distributed
!! 3-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_3d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scale, leave_file_open)
  character(len=*),       intent(in)    :: filename !< name of the netcdf file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:,:), intent(inout) :: u_data    !< The 3 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:,:), intent(inout) :: v_data    !< The 3 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  integer :: is, ie, js, je, i, dim_unlim, ndims
  integer :: u_pos, v_pos
  integer, allocatable :: dim_sizes_u(:), dim_sizes_v(:)
  character(len=32), allocatable :: dim_names_u(:), dim_names_v(:), units_u(:), units_v(:)
  character(len=1) :: x_or_y
  logical :: is_valid
  logical :: file_open_success ! .true. if open file is successful
  logical :: close_the_file ! indicates whether to close the file after MOM_read_vector is called; default is .true.

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
    if (.not. file_open_success) &
      call MOM_error(FATAL, "MOM_read_vector_3d: netcdf file "//trim(filename)//" not opened.")
  endif

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  ndims = get_variable_num_dimensions(fileobj_read_dd, u_fieldname)
  allocate(dim_sizes_u(ndims))
  allocate(dim_sizes_v(ndims))
  allocate(dim_names_u(ndims))
  allocate(dim_names_v(ndims))
  allocate(units_u(ndims))
  allocate(units_v(ndims))

  units_u(:) = ""
  units_v(:) = ""
  dim_names_u(:) = ""
  dim_names_v(:) = ""

  call get_variable_size(fileobj_read_dd, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileobj_read_dd, v_fieldname, dim_sizes_v, broadcast=.true.)
  call get_variable_dimension_names(fileobj_read_dd, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileobj_read_dd, v_fieldname, dim_names_v, broadcast=.true.)

  do i=1,ndims
    ! register the u axes
    if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_u(i)))) then
      call get_variable_units(fileobj_read_dd, dim_names_u(i), units_u(i))
      call validate_lat_lon_units(units_u(i), x_or_y, is_valid) 
      if (is_valid) then
        call register_axis(fileobj_read_dd, dim_names_u(i), x_or_y, domain_position=u_pos)
      else
        call register_axis(fileobj_read_dd, dim_names_u(i), dim_sizes_u(i))
      endif
    endif
    ! Register the v axes if they differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then
      if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_v(i)))) then
        call get_variable_units(fileobj_read_dd, dim_names_v(i), units_v(i))
        call validate_lat_lon_units(units_v(i), x_or_y, is_valid) 
        if (is_valid) then
          call register_axis(fileobj_read_dd, dim_names_v(i), x_or_y, domain_position=v_pos)
        else
          call register_axis(fileobj_read_dd, dim_names_v(i), dim_sizes_v(i))
        endif
      endif
    endif
  enddo
  ! read the data
  dim_unlim = 0
  if (present(timelevel)) then
    do i=1,ndims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names_u(i))) then
        dim_unlim = i
        exit
      endif
    enddo
    if (dim_unlim .gt. 0) then
      call read_data(fileobj_read_dd, u_fieldname, u_data, unlim_dim_level=timelevel)
      call read_data(fileobj_read_dd, v_fieldname, v_data, unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read_dd, u_fieldname, u_data, edge_lengths=dim_sizes_u)
      call read_data(fileobj_read_dd, v_fieldname, v_data, edge_lengths=dim_sizes_v) 
    endif
  else
    call read_data(fileobj_read_dd, u_fieldname, u_data, edge_lengths=dim_sizes_u)
    call read_data(fileobj_read_dd, v_fieldname, v_data, edge_lengths=dim_sizes_v)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
  endif
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(u_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(u_data,2), js, je)
    u_data(is:ie,js:je,:) = scale*u_data(is:ie,js:je,:)
    call get_simple_array_i_ind(MOM_Domain, size(v_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(v_data,2), js, je)
    v_data(is:ie,js:je,:) = scale*v_data(is:ie,js:je,:)
  endif ; endif
  if (allocated(dim_names_u)) deallocate(dim_names_u)
  if (allocated(dim_names_v)) deallocate(dim_names_v)
  if (allocated(dim_sizes_u)) deallocate(dim_sizes_u)
  if (allocated(dim_sizes_v)) deallocate(dim_sizes_v)
  if (allocated(units_u)) deallocate(units_u)
  if (allocated(units_v)) deallocate(units_v)
end subroutine MOM_read_vector_3d

!> register a MOM diagnostic axis to a domain-decomposed file
subroutine MOM_register_diagnostic_axis(fileObj, axisName, axisLength)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in) :: axisName !< name of the axis to register to file
  integer, intent(in), optional :: axisLength !< length of axis/dimension ;only needed for Layer, Interface, Time,
                                              !! Period
  select case (trim(axisName))
    case ('latq'); call register_axis(fileObj,'latq','y', domain_position=NORTH_FACE)
    case ('lath'); call register_axis(fileObj,'lath','y', domain_position=CENTER)
    case ('lonq'); call register_axis(fileObj,'lonq','x', domain_position=EAST_FACE)
    case ('lonh'); call register_axis(fileObj,'lonh','x', domain_position=CENTER)
    case default
      if (.not. present(axisLength)) call MOM_error(FATAL,"MOM_io:register_diagnostic_axis: "//&
                        "An axis_length argument is required to register the axis "//trim(axisName))
      call register_axis(fileObj, trim(axisName), axisLength)
  end select
end subroutine MOM_register_diagnostic_axis

!> register axes associated with a variable from a domain-decomposed netCDF file that are mapped to
!! a sub-domain (e.g., a supergrid).
!> \note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes to obtain
!! the correct domain decomposition for the data buffer.
subroutine MOM_register_variable_axes_subdomain(fileObj, variableName, io_domain, xPosition, yPosition)
  type(FmsNetcdfFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in) :: variableName !< name of the variable
  type(domain2D), intent(in) :: io_domain !< type that contains the mpp io domain
  integer, intent(in), optional :: xPosition !< domain position of the x-axis
  integer, intent(in), optional :: yPosition !< domain position of the y-axi
  ! local
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i, isg, ieg, isc, iec, jsg, jeg, jsc, jec, xlen, ylen
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos, pos ! domain positions for x and y axes. Default is CENTER
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:register_variable_axes_subdomain: The fileObj "// &
                                                  " has not been opened. Call fms2_open_file(fileObj,...) "// &
                                                  "before passing the fileObj argument to this function.")
  xPos=CENTER
  yPos=CENTER
  if (present(xPosition)) xPos=xPosition
  if (present(yPosition)) yPos=yPosition
  ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  call get_variable_size(fileObj, trim(variableName), dimSizes, broadcast=.true.)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)
  ! determine the position to pass to the mpp domain calls
  if (xPos .eq. EAST_FACE) then
    if (yPos .eq. NORTH_FACE) then
      pos = CORNER
    else
      pos = EAST_FACE
    endif
  elseif (yPos .eq. NORTH_FACE) then
    pos = NORTH_FACE
  endif
  ! Get the lengths of the global indicies   
  call mpp_get_compute_domain(io_domain, xsize=xlen, ysize=ylen, position=pos)
  ! register the axes
  !>\note: This is not a comprehensive check for all possible supported horizontal axes associated with variables
  !! read from netCDF files. Developers should add/remove cases as needed.
  do i=1,ndims
    !if (.not.(is_dimension_registered(fileObj, trim(dim_names(i))))) then
      select case(trim(lowercase(dim_names(i))))
        case ("grid_x_t")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case ("nx")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("nxp")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("longitude")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("long")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("lon")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("lonh")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("lonq")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("xh")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case ("grid_y_t")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case ("ny")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("nyp")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("latitude")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("lat")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("lath")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("latq")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("yh")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case default ! assumes that the axis is not domain-decomposed
          if (.not. is_dimension_unlimited(fileObj, trim(dim_names(i)))) &
            call MOM_error(WARNING,"MOM_register_variable_axes_subdomain: the axis "//trim(dim_names(i))//&
              "is not included in the valid x and y dimension cases. If the code hangs, check the whether "//&
              "an x or y axis is being registered as a non-domain-decomposed variable, "//&
              "and add it to the accepted cases if necessary.")
          call register_axis(fileObj, trim(dim_names(i)), dimSizes(i))
      end select
    endif
  enddo

  if (allocated(dimSizes)) deallocate(dimSizes)
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_register_variable_axes_subdomain


!> register axes associated with a variable from a domain-decomposed netCDF file
!> \note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes to obtain
!! the correct domain decomposition for the data buffer.
subroutine MOM_register_variable_axes_full(fileObj, variableName, xPosition, yPosition)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in) :: variableName !< name of the variable
  integer, intent(in), optional :: xPosition !< domain position of the x-axis
  integer, intent(in), optional :: yPosition !< domain position of the y-axis
  ! local
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos ! domain positions for x and y axes. Default is CENTER
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:register_variable_axes: The fileObj has "// &
                                                  "not been opened. Call fms2_open_file(fileObj,...) before "// &
                                                  "passing the fileObj argument to this function.")
  xPos=CENTER
  yPos=CENTER
  if (present(xPosition)) xPos=xPosition
  if (present(yPosition)) yPos=yPosition

  ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  call get_variable_size(fileObj, trim(variableName), dimSizes, broadcast=.true.)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)
  ! register the axes
  !>\note: This is not a comprehensive check for all possible supported horizontal axes associated with variables
  !! read from netCDF files. Developers should add/remove cases as needed.
  do i=1,ndims
    if (.not.(is_dimension_registered(fileobj, trim(dim_names(i))))) then
      select case(trim(lowercase(dim_names(i))))
        case ("grid_x_t")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case ("nx")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("nxp")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("longitude")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("long")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("lon")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("lonh")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("lonq")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("xh")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case ("grid_y_t")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case ("ny")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("nyp")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("latitude")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("lat")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("lath")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("latq")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("yh")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=xPos)
        case default ! assumes that the axis is not domain-decomposed
          if (.not. is_dimension_unlimited(fileObj, trim(dim_names(i)))) &
            call MOM_error(WARNING,"MOM_register_variable_axes_full: the axis "//trim(dim_names(i))//" is not "//&
              "included in the valid x and y dimension cases. If the code hangs, check the whether "//&
              "an x or y axis is being registered as a non-domain-decomposed variable, "//&
              "and add it to the accepted cases if necessary.")
          call register_axis(fileObj, trim(dim_names(i)), dimSizes(i))
      end select
    endif
  enddo

  deallocate(dimSizes)
  deallocate(dim_names)
end subroutine MOM_register_variable_axes_full

!> Get the horizontal grid, vertical grid, and/or time

!> Get the horizontal grid, vertical grid, and/or time dimension names and lengths
!! for a single variable from the hor_grid, t_grid, and z_grid values returned by a prior call to query_vardesc
subroutine get_var_dimension_features(hor_grid, z_grid, t_grid_in, &
                                      dim_names, dim_lengths, num_dims, G, dG, GV)

  character(len=*), intent(in) :: hor_grid !< horizontal grid
  character(len=*), intent(in) :: z_grid !< vertical grid
  character(len=*), intent(in) :: t_grid_in !< time grid
  character(len=*), dimension(:), intent(inout) :: dim_names !< array of dimension names
  integer, dimension(:), intent(inout) :: dim_lengths !< array of dimension sizes
  integer, intent(inout) :: num_dims !< number of axes to register in the restart file
  type(ocean_grid_type),  optional, intent(in) :: G !< The ocean's grid structure
  type(dyn_horgrid_type), optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                      !! is required if the new file uses any
                                                      !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure

  ! local
  logical :: use_lath
  logical :: use_lonh
  logical :: use_latq
  logical :: use_lonq
  character(len=8) :: t_grid
  character(len=8) :: t_grid_read
  integer :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  !integer :: npes
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()
  type(MOM_domain_type), pointer :: domain => NULL() ! Domain used to get the pe count

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.

  ! set the ocean grid coordinates

  if (present(G)) then
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB

    call get_horizontal_grid_logic(hor_grid, use_lath, use_lonh, use_latq, use_lonq)
  elseif (present(dG)) then
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB

    call get_horizontal_grid_logic(hor_grid, use_lath, use_lonh, use_latq, use_lonq)
  endif

  ! add longitude name to dimension name array
  if (use_lonh) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("lonh")) = "lonh"
    dim_lengths(num_dims) = size(gridLonT(isg:ieg))
  elseif (use_lonq) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("lonq")) = "lonq"
    dim_lengths(num_dims) = size(gridLonB(IsgB:IegB))
  endif
  ! add latitude name to dimension name array
  if (use_lath) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("lath")) = "lath"
    dim_lengths(num_dims) = size(gridLatT(jsg:jeg))
  elseif (use_latq) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("latq")) = "latq"
    dim_lengths(num_dims) = size(gridLatB(JsgB:JegB))
  endif

  if (present(GV)) then
    ! vertical grid
    select case (trim(z_grid))
      case ('L')
        num_dims = num_dims+1
        dim_names(num_dims) = ""
        dim_names(num_dims)(1:len_trim("Layer")) = "Layer"
        dim_lengths(num_dims) = GV%ke
      case ('i')
        num_dims = num_dims+1
        dim_names(num_dims) = ""
        dim_names(num_dims)(1:len_trim("Interface")) = "Interface"
        dim_lengths(num_dims) = GV%ke+1
      case ('1') ! Do nothing.
      case default
      call MOM_error(FATAL, "MOM_io: get_var_dimension_features: "//&
                     " has an unrecognized z_grid argument"//trim(z_grid))
    end select
  endif
  ! time
  t_grid = adjustl(t_grid_in)
  select case (t_grid(1:1))
    case ('s', 'a', 'm')
      num_dims = num_dims+1
      dim_names(num_dims) = ""
      dim_names(num_dims)(1:len_trim("Time")) = "Time"
      dim_lengths(num_dims) = unlimited
    case ('p')
      if (len_trim(t_grid(2:8)) <= 0) then
          call MOM_error(FATAL,"MOM_io:get_var_dimension_features: "//&
                           "No periodic axis length was specified in "//trim(t_grid))
      endif
      num_dims = num_dims+1
      dim_names(num_dims) = ""
      dim_names(num_dims)(1:len_trim("Period")) = "Period"
      dim_lengths(num_dims) = unlimited
    case ('1') ! Do nothing.
    case default
      call MOM_error(WARNING, "MOM_io: get_var_dimension_features: "//&
                     "Unrecognized t_grid "//trim(t_grid))
  end select
end subroutine get_var_dimension_features

!> Populate the axis_data structure with axis data and attributes for diagnostic and restart files
subroutine MOM_get_diagnostic_axis_data(axis_data_CS, axis_name, axis_number, G, dG, GV, time_val, time_units)

  type(axis_data_type), intent(inout) :: axis_data_CS !< structure containing the axis data and metadata
  character(len=*), intent(in) :: axis_name !< name of the axis
  integer, intent(in) :: axis_number !< positional value (wrt to file) of the axis to register
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                      !! is required if the file uses any
                                                      !! horizontal grid axes.
  type(verticalGrid_type), target, optional, intent(in) :: GV !< ocean vertical grid structure
  real,dimension(:), target, optional, intent(in) :: time_val !< time value
  character(len=*), optional,intent(in) :: time_units!< units for non-periodic time axis
  ! local
  character(len=40) :: x_axis_units='', y_axis_units=''
  integer :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()

  ! initialize axis_data_CS elements
  axis_data_CS%axis(axis_number)%name = ''
  axis_data_CS%axis(axis_number)%longname = ''
  axis_data_CS%axis(axis_number)%units = ''
  axis_data_CS%axis(axis_number)%horgrid_position = 0
  axis_data_CS%axis(axis_number)%is_domain_decomposed = .false.
  axis_data_CS%axis(axis_number)%positive = ''
  axis_data_CS%data(axis_number)%p => NULL()

  ! set the ocean grid coordinates and metadata
  if (present(G)) then
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif

  select case(trim(axis_name))
    case('lath')
      if (associated(gridLatT)) &
        axis_data_CS%data(axis_number)%p=>gridLatT(jsg:jeg)

      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Latitude'
      axis_data_CS%axis(axis_number)%units = y_axis_units
      axis_data_CS%axis(axis_number)%horgrid_position = CENTER
      axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
    case('lonh')
      if (associated(gridLonT)) &
        axis_data_CS%data(axis_number)%p=>gridLonT(isg:ieg)

      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%horgrid_position  = CENTER
      axis_data_CS%axis(axis_number)%longname = 'Longitude'
      axis_data_CS%axis(axis_number)%units = x_axis_units
      axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
    case('latq')
      if (associated(gridLatB)) &
        axis_data_CS%data(axis_number)%p=>gridLatB(JsgB:JegB)

      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Latitude'
      axis_data_CS%axis(axis_number)%units = y_axis_units
      axis_data_CS%axis(axis_number)%horgrid_position = NORTH_FACE
      axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
    case('lonq')
      if (associated(gridLonB)) &
        axis_data_CS%data(axis_number)%p=>gridLonB(IsgB:IegB)

        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = EAST_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
    case('Layer')
      if (present(GV)) then
        axis_data_CS%data(axis_number)%p=>GV%sLayer(1:GV%ke)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Layer pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
      endif
    case('Interface')
      if (present(GV)) then
        axis_data_CS%data(axis_number)%p=>GV%sInterface(1:GV%ke+1)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Interface pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
      endif
    case('Time')
      if (.not.(present(time_val))) &
           call MOM_error(FATAL, "MOM_io::get_diagnostic_axis_data: requires time_val"//&
                          " and time_units arguments for "//trim(axis_name))

      axis_data_CS%data(axis_number)%p=>time_val
      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Time'

      if (present(time_units)) then
        axis_data_CS%axis(axis_number)%units = time_units
      else
        axis_data_CS%axis(axis_number)%units = 'days'
      endif
    case('Period')
      if (.not.(present(time_val))) &
        call MOM_error(FATAL, "MOM_io::get_diagnostic_axis_data: requires a time_val argument "// &
                       "for "//trim(axis_name))
      axis_data_CS%data(axis_number)%p=>time_val
      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Periods for cyclical variables'
    case default
      call MOM_error(WARNING, "MOM_io::get_diagnostic_axis_data:"//trim(axis_name)//"is an unrecognized axis")
  end select

end subroutine MOM_get_diagnostic_axis_data

!> check a netcdf file to see whether it contains the specified field
function field_exists(file_name, variable_name) result (var_exists)
  character(len=*),  intent(in) :: file_name   !< name of the file to check
  character(len=*),  intent(in) :: variable_name  !< name of the variable to check for

  ! local
  logical :: var_exists ! .true. if variable is found in file
  logical :: file_open_success ! .true. if call to open_file is successful
  type(fmsNetcdfFile_t) :: fileobj ! netcdf file object returned by open_file

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = fms2_open_file(fileobj, file_name, "read", is_restart=.false.)
  var_exists = variable_exists(fileobj, trim(variable_name))
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)
end function field_exists

!> check that latitude or longitude units are valid CF-compliant values
!! return true or false and x_or_y character value corresponding to the axis direction
subroutine validate_lat_lon_units(unit_string, x_or_y, units_are_valid)
character(len=*), intent(in) :: unit_string !< string of units 
character(len=1), intent(out) :: x_or_y !< "x" for longitude or "y" latitude
logical, intent(out) :: units_are_valid !< .true. if units match acceptable values; default is .false.

select case (lowercase(trim(unit_string)))
  case ("degrees_north"); units_are_valid = .true.; x_or_y = "y"
  case ("degree_north"); units_are_valid = .true.; x_or_y = "y"
  case ("degrees_n"); units_are_valid = .true.; x_or_y = "y"
  case ("degree_n"); units_are_valid = .true.; x_or_y = "y"
  case ("degreen"); units_are_valid = .true.; x_or_y = "y"
  case ("degreesn"); units_are_valid = .true.; x_or_y = "y"
  case ("degrees_east"); units_are_valid = .true.; x_or_y = "x"
  case ("degree_east"); units_are_valid = .true.;x_or_y = "x"
  case ("degreese"); units_are_valid = .true.; x_or_y =  "x"
  case ("degreee"); units_are_valid = .true.; x_or_y =  "x"
  case ("degree_e"); units_are_valid = .true.; x_or_y =  "x"
  case ("degrees_e"); units_are_valid = .true.; x_or_y = "x"
  case default; units_are_valid = .false.; x_or_y = ""
end select

end subroutine validate_lat_lon_units

!> get the dimesion sizes of a variable
subroutine field_size(file_name, variable_name, dim_sizes)
  character(len=*),  intent(in) :: file_name   !< name of the file to check
  character(len=*),  intent(in) :: variable_name  !< name of the variable to check for
  integer, dimension(:), intent(inout) :: dim_sizes !< variable dimension sizes
  ! local
  logical :: file_open_success ! .true. if call to open_file is successful
  type(fmsNetcdfFile_t) :: fileobj ! netcdf file object returned by open_file
  integer :: i, ndims

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = fms2_open_file(fileobj, file_name, "read", is_restart=.false.)
  if (.not. (variable_exists(fileobj, trim(variable_name)))) &
    call MOM_error(FATAL, "MOM_io::field_size: variable "//trim(variable_name)// &
                   " not found in file "//trim(file_name))

  ndims = get_variable_num_dimensions(fileObj, trim(variable_name))
  call get_variable_size(fileObj, trim(variable_name), dim_sizes(1:ndims))

  if (check_if_open(fileobj)) call fms2_close_file(fileobj)
end subroutine field_size

!> set the logical variables that determine which diagnositic axes to use
subroutine get_horizontal_grid_logic(grid_string_id, use_lath, use_lonh, use_latq, use_lonq)
  character(len=*), intent(in) :: grid_string_id !< horizontal grid string
  logical, intent(out) :: use_lath !< if .true., y-axis is oriented in CENTER position
  logical, intent(out) :: use_lonh !< if .true., x-axis is oriented in CENTER position
  logical, intent(out) :: use_latq !< if .true., y-axis is oriented in NORTH_FACE position
  logical, intent(out) :: use_lonq !< if .true., x-axis is oriented in EAST_FACE position

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.
  select case (trim(grid_string_id))
     case ('h') ; use_lath = .true. ; use_lonh = .true. ! x=CENTER, y=CENTER
     case ('q') ; use_latq = .true. ; use_lonq = .true. ! x=EAST_FACE, y=NORTH_FACE
     case ('u') ; use_lath = .true. ; use_lonq = .true. ! x=EAST_FACE, y=CENTER
     case ('v') ; use_latq = .true. ; use_lonh = .true. ! x=CENTER, y=NORTH_FACE
     case ('T')  ; use_lath = .true. ; use_lonh = .true. ! x=CENTER, y=CENTER
     case ('Bu') ; use_latq = .true. ; use_lonq = .true. ! x=EAST_FACE, y=NORTH_FACE
     case ('Cu') ; use_lath = .true. ; use_lonq = .true. ! x=EAST_FACE, y=CENTER
     case ('Cv') ; use_latq = .true. ; use_lonh = .true. ! x=CENTER, y=NORTH_FACE
     case ('1') ; ! x=0, y=0
     case default
        call MOM_error(FATAL, "MOM_io:get_var_dimension_features "//&
                        "Unrecognized hor_grid argument "//trim(grid_string_id))
  end select
end subroutine get_horizontal_grid_logic

!> get the size of a variable in bytes
function get_variable_byte_size(hor_grid, z_grid, t_grid, G, num_zlevels) result(var_sz)
  character(len=*), intent(in) :: hor_grid !< horizontal grid string
  character(len=*), intent(in) :: z_grid !< vertical grid string
  character(len=*), intent(in) :: t_grid !< time string
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure;
  integer, intent(in) :: num_zlevels     !< number of vertical levels
  ! local
  integer(kind=8) :: var_sz !< The size in bytes of each variable
  integer :: var_periods
  character(len=8) :: t_grid_read=''

  var_periods = 0

  if (trim(hor_grid) == '1') then
     var_sz = 8
  else
     var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
  endif

  select case (trim(z_grid))
     case ('L') ; var_sz = var_sz * num_zlevels
     case ('i') ; var_sz = var_sz * (num_zlevels+1)
  end select

  if (adjustl(t_grid(1:1)) == 'p') then
     if (len_trim(t_grid(2:8)) > 0) then
        var_periods = -1
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods > 1) var_sz = var_sz * var_periods
     endif
  endif

end function get_variable_byte_size

!> Read the data associated with a named axis in a file
!>\todo: this routine is redundant; replace calls to this routine with MOM_read_data
subroutine read_axis_data(filename, axis_name, var)
  character(len=*),   intent(in)  :: filename  !< Name of the file to read
  character(len=*),   intent(in)  :: axis_name !< Name of the axis to read
  real, dimension(:), intent(out) :: var       !< The axis location data

  ! local
  logical :: file_open_success ! .true. if call to open_file is successful
  type(fmsNetcdfFile_t) :: fileobj ! netcdf file object returned by open_file
  integer :: i, ndims
  integer, dimension(1) :: dim_size ! variable dimension size

  var(:) = -1.0
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = fms2_open_file(fileobj, filename, "read", is_restart=.false.)
  if (.not. (variable_exists(fileobj, trim(axis_name)))) &
    call MOM_error(FATAL, "MOM_io::read_axis_data: variable "//trim(axis_name)// &
                   " not found in file "//trim(filename))

  ndims = get_variable_num_dimensions(fileObj, trim(axis_name))
  if (ndims .ne. 1) &
    call MOM_error(FATAL, "MOM_io::read_axis_data: Variable "//trim(axis_name)// " must be a 1-D array.")

  call get_variable_size(fileobj, trim(axis_name), dim_size, broadcast=.true.)
  call read_data(fileobj, trim(axis_name),var(1:dim_size(1)))

  if (check_if_open(fileobj)) call fms2_close_file(fileobj)

end subroutine read_axis_data

!> Define the time units for the input time value
function get_time_units(time_value) result(time_units_out)
   real, intent(in) :: time_value !< numerical time value in seconds
                                  !! i.e., before dividing by 86400.
   ! local
   character(len=10) :: time_units ! time units
   character(len=10) :: time_units_out ! time units trimmed
   time_units = ''
   time_units_out = ''
   if (time_value < 0.0) then
      time_units = "days" ! The default value.
   elseif (mod(time_value,86400.0)==0.0) then
      time_units = "days"
   elseif ((time_value >= 0.99) .and. (time_value < 1.01)) then
              time_units = "seconds"
   elseif ((time_value >= 3599.0) .and. (time_value < 3601.0)) then
              time_units = "hours"
   elseif ((time_value >= 86399.0) .and. (time_value < 86401.0)) then
              time_units = "days"
   elseif ((time_value >= 3.0e7) .and. (time_value < 3.2e7)) then
              time_units = "years"
   else
       write(time_units,'(es8.2," s")') time_value
   endif
   time_units_out = trim(time_units)
end function get_time_units

!> This function determines how many time levels a variable has.
function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename   !< name of the file to read
  character(len=*),  intent(in) :: varname    !< variable whose number of time levels
                                              !! are to be returned
  integer, optional, intent(in) :: min_dims   !< The minimum number of dimensions a variable must have
                                              !! if it has a time dimension.  If the variable has 1 less
                                              !! dimension than this, then 0 is returned.
  integer :: n_time                           !< number of time levels varname has in filename

  ! local
  logical :: file_open_success ! .true. if call to open_file is successful
  logical :: variableExists ! .true. if variable is found in file
  character(len=200) :: msg
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i, ndims
  type(fmsNetcdfFile_t) :: fileobj !netcdf file object returned by open_file

  n_time = -1

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = fms2_open_file(fileobj, filename, "read", is_restart=.false.)

  ! check that variable is in the file
  if (.not.(variable_exists(fileobj, trim(varname)))) call MOM_error(FATAL, "num_time_levels: variable"//&
    trim(varname)//" not found in "//trim(filename))
  ! get the number of variable dimensions
  ndims = get_variable_num_dimensions(fileobj, trim(varname))

  if (present(min_dims)) then
    if (ndims .LT. min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_time_levels: variable "//trim(varname)//&
        " in file "//trim(filename)//" has fewer than min_dims = "//trim(msg)//&
        " dimensions.")
    elseif (ndims .EQ. min_dims - 1) then
      n_time = 0 ; return
    endif
  endif
  ! check for the unlimited dimension and set n_time to the length of the unlimited dimension
  allocate(dim_names(ndims))

  call get_variable_dimension_names(fileobj, trim(varname), dim_names)
  do i=1,ndims
    if (is_dimension_unlimited(fileobj, trim(dim_names(i)))) &
      call get_dimension_size(fileobj, trim(dim_names(i)), n_time)
  enddo

  deallocate(dim_names)

end function num_timelevels

!> Returns a vardesc type whose elements have been filled with the provided
!! fields.  The argument name is required, while the others are optional and
!! have default values that are empty strings or are appropriate for a 3-d
!! tracer field at the tracer cell centers.
function var_desc(name, units, longname, hor_grid, z_grid, t_grid, &
                  cmor_field_name, cmor_units, cmor_longname, conversion, caller) result(vd)
  character(len=*),           intent(in) :: name               !< variable name
  character(len=*), optional, intent(in) :: units              !< variable units
  character(len=*), optional, intent(in) :: longname           !< variable long name
  character(len=*), optional, intent(in) :: hor_grid           !< variable horizonal staggering
  character(len=*), optional, intent(in) :: z_grid             !< variable vertical staggering
  character(len=*), optional, intent(in) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(in) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(in) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in) :: cmor_longname      !< CMOR long name
  real            , optional, intent(in) :: conversion         !< for unit conversions, such as needed to
                                                               !! convert from intensive to extensive
  character(len=*), optional, intent(in) :: caller             !< calling routine?
  type(vardesc)                          :: vd                 !< vardesc type that is created

  character(len=120) :: cllr
  cllr = "var_desc"
  if (present(caller)) cllr = trim(caller)

  call safe_string_copy(name, vd%name, "vd%name", cllr)

  vd%longname = "" ; vd%units = ""
  vd%hor_grid = 'h' ; vd%z_grid = 'L' ; vd%t_grid = 's'

  vd%cmor_field_name  =  ""
  vd%cmor_units       =  ""
  vd%cmor_longname    =  ""
  vd%conversion       =  1.0

  call modify_vardesc(vd, units=units, longname=longname, hor_grid=hor_grid, &
                      z_grid=z_grid, t_grid=t_grid,                          &
                      cmor_field_name=cmor_field_name,cmor_units=cmor_units, &
                      cmor_longname=cmor_longname, conversion=conversion, caller=cllr)

end function var_desc

!> This routine modifies the named elements of a vardesc type.
!! All arguments are optional, except the vardesc type to be modified.
subroutine modify_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                 cmor_field_name, cmor_units, cmor_longname, conversion, caller)
  type(vardesc),              intent(inout) :: vd              !< vardesc type that is modified
  character(len=*), optional, intent(in)    :: name            !< name of variable
  character(len=*), optional, intent(in)    :: units           !< units of variable
  character(len=*), optional, intent(in)    :: longname        !< long name of variable
  character(len=*), optional, intent(in)    :: hor_grid        !< horizonal staggering of variable
  character(len=*), optional, intent(in)    :: z_grid          !< vertical staggering of variable
  character(len=*), optional, intent(in)    :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in)    :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in)    :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in)    :: cmor_longname   !< CMOR long name
  real            , optional, intent(in)    :: conversion      !< for unit conversions, such as needed
                                                               !! to convert from intensive to extensive
  character(len=*), optional, intent(in)    :: caller          !< calling routine?

  character(len=120) :: cllr
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(name, vd%name, "vd%name", cllr)

  if (present(longname))  call safe_string_copy(longname, vd%longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(units, vd%units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(hor_grid, vd%hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(z_grid, vd%z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(t_grid, vd%t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(cmor_field_name, vd%cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(cmor_units, vd%cmor_units, &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(cmor_longname, vd%cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

end subroutine modify_vardesc

!> This function returns the CMOR standard name given a CMOR longname, based on
!! the standard pattern of character conversions.
function cmor_long_std(longname) result(std_name)
  character(len=*), intent(in) :: longname  !< The CMOR longname being converted
  character(len=len(longname)) :: std_name  !< The CMOR standard name generated from longname

  integer :: k

  std_name = lowercase(longname)

  do k=1, len_trim(std_name)
    if (std_name(k:k) == ' ') std_name(k:k) = '_'
  enddo

end function cmor_long_std

!> This routine queries vardesc
subroutine query_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                         cmor_field_name, cmor_units, cmor_longname, conversion, caller)
  type(vardesc),              intent(in)  :: vd                 !< vardesc type that is queried
  character(len=*), optional, intent(out) :: name               !< name of variable
  character(len=*), optional, intent(out) :: units              !< units of variable
  character(len=*), optional, intent(out) :: longname           !< long name of variable
  character(len=*), optional, intent(out) :: hor_grid           !< horiz staggering of variable
  character(len=*), optional, intent(out) :: z_grid             !< vert staggering of variable
  character(len=*), optional, intent(out) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(out) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(out) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(out) :: cmor_longname      !< CMOR long name
  real            , optional, intent(out) :: conversion         !< for unit conversions, such as needed to
                                                                !! convert from intensive to extensive
  character(len=*), optional, intent(in)  :: caller             !< calling routine?


  character(len=120) :: cllr
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(vd%name, name,         &
                               "vd%name of "//trim(vd%name), cllr)
  if (present(longname))  call safe_string_copy(vd%longname, longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(vd%units, units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(vd%hor_grid, hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(vd%z_grid, z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(vd%t_grid, t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(vd%cmor_field_name, cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(vd%cmor_units, cmor_units,          &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(vd%cmor_longname, cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

end subroutine query_vardesc


!> Copies a string
subroutine safe_string_copy(str1, str2, fieldnm, caller)
  character(len=*),           intent(in)  :: str1    !< The string being copied
  character(len=*),           intent(out) :: str2    !< The string being copied into
  character(len=*), optional, intent(in)  :: fieldnm !< The name of the field for error messages
  character(len=*), optional, intent(in)  :: caller  !< The calling routine for error messages

  if (len(trim(str1)) > len(str2)) then
    if (present(fieldnm) .and. present(caller)) then
      call MOM_error(FATAL, trim(caller)//" attempted to copy the overly long"//&
        " string "//trim(str1)//" into "//trim(fieldnm))
    else
      call MOM_error(FATAL, "safe_string_copy: The string "//trim(str1)//&
                     " is longer than its intended target.")
    endif
  endif
  str2 = trim(str1)
end subroutine safe_string_copy


!> Returns a name with "%#E" or "%E" replaced with the ensemble member number.
function ensembler(name, ens_no_in) result(en_nm)
  character(len=*),  intent(in) :: name       !< The name to be modified
  integer, optional, intent(in) :: ens_no_in  !< The number of the current ensemble member
  character(len=len(name)) :: en_nm  !< The name encoded with the ensemble number

  ! This function replaces "%#E" or "%E" with the ensemble number anywhere it
  ! occurs in name, with %E using 4 or 6 digits (depending on the ensemble size)
  ! and %#E using # digits, where # is a number from 1 to 9.

  character(len=len(name)) :: tmp
  character(10) :: ens_num_char
  character(3)  :: code_str
  integer :: ens_no
  integer :: n, is, ie

  en_nm = trim(name)
  if (index(name,"%") == 0) return

  if (present(ens_no_in)) then
    ens_no = ens_no_in
  else
    ens_no = get_ensemble_id()
  endif

  write(ens_num_char, '(I10)') ens_no ; ens_num_char = adjustl(ens_num_char)
  do
    is = index(en_nm,"%E")
    if (is == 0) exit
    if (len(en_nm) < len(trim(en_nm)) + len(trim(ens_num_char)) - 2) &
      call MOM_error(FATAL, "MOM_io ensembler: name "//trim(name)// &
      " is not long enough for %E expansion for ens_no "//trim(ens_num_char))
    tmp = en_nm(1:is-1)//trim(ens_num_char)//trim(en_nm(is+2:))
    en_nm = tmp
  enddo

  if (index(name,"%") == 0) return

  write(ens_num_char, '(I10.10)') ens_no
  do n=1,9 ; do
    write(code_str, '("%",I1,"E")') n

    is = index(en_nm,code_str)
    if (is == 0) exit
    if (ens_no < 10**n) then
      if (len(en_nm) < len(trim(en_nm)) + n-3) call MOM_error(FATAL, &
        "MOM_io ensembler: name "//trim(name)//" is not long enough for %E expansion.")
      tmp = en_nm(1:is-1)//trim(ens_num_char(11-n:10))//trim(en_nm(is+3:))
    else
      call MOM_error(FATAL, "MOM_io ensembler: Ensemble number is too large "//&
          "to be encoded with "//code_str//" in "//trim(name))
    endif
    en_nm = tmp
  enddo ; enddo

end function ensembler

subroutine scale_data_1d(data, scale_factor)
  real, dimension(:), intent(inout) :: data !< The 1-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor

  if (scale_factor /= 1.0) then
    data(:) = scale_factor*data(:)
  endif
end subroutine scale_data_1d

subroutine scale_data_2d(data, scale_factor, MOM_domain)
  real, dimension(:,:), intent(inout) :: data !< The 2-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    if (present(MOM_domain)) then
      call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
      call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
      data(is:ie,js:je) = scale_factor*data(is:ie,js:je)
    else
      data(:,:) = scale_factor*data(:,:)
    endif
  endif
end subroutine scale_data_2d

subroutine scale_data_3d(data, scale_factor, MOM_domain)
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    if (present(MOM_domain)) then
      call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
      call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
      data(is:ie,js:je,:) = scale_factor*data(is:ie,js:je,:)
    else
      data(:,:,:) = scale_factor*data(:,:,:)
    endif
  endif
end subroutine scale_data_3d

subroutine scale_data_4d(data, scale_factor, MOM_domain)
  real, dimension(:,:,:,:), intent(inout) :: data !< The 4-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    if (present(MOM_domain)) then
      call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
      call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
      data(is:ie,js:je,:,:) = scale_factor*data(is:ie,js:je,:,:)
    else
      data(:,:,:,:) = scale_factor*data(:,:,:,:)
    endif
  endif
end subroutine scale_data_4d

!> convert the variable checksum integer(s) to a single string
!! If there is more than 1 checksum, commas are inserted between
!! each checksum value in the output string
function convert_checksum_to_string(checksum_int) result (checksum_string)
  integer(kind=8), intent(in) :: checksum_int !< checksum integer values
! local
  character(len=64) :: checksum_string
  integer :: i

  checksum_string = ''

  write (checksum_string,'(Z16)') checksum_int ! Z16 is the hexadecimal format code

end function convert_checksum_to_string

!> convert the variable checksum string read in from a restart file
!! to integer value(s)
function convert_checksum_string_to_int(checksum_char) result(checksum_file)
  character(len=*), intent(in) :: checksum_char !< checksum character array
  ! local
  integer(kind=8),dimension(3) :: checksum_file !< checksum string corresponds to
                                                   !< values from up to 3 times
  integer :: last
  integer :: start
  integer(kind=8) :: checksumh
  integer :: num_checksumh
  integer :: k

  start =0
  last = 0
  checksumh = 0
  num_checksumh = 1
  last = len_trim(checksum_char)
  start = index(trim(checksum_char),",") ! A start value of 0 implies only 1 checksum value
  ! Scan checksum character array for the ',' delimiter, which indicates that the corresponding variable
  ! has multiple time levels.
  do while ((start > 0) .and. (start < (last-15)))
     start = start + scan(checksum_char(start:last), "," ) ! move starting pointer after ","
     num_checksumh = num_checksumh + 1
  enddo

  start = 1

  do k = 1, num_checksumh
     read(checksum_char(start:start+15),'(Z16)') checksumh ! Z=hexadecimal integer: Z16 is for 64-bit data types
     checksum_file(k) = checksumh
     start = start+ 17 ! Move start index past the ',' in checksum_char
  enddo

end function convert_checksum_string_to_int

!> function to get the index of a time_value from a netCDF file
function get_time_index(filename, time_to_find) result (time_index)
  character(len=*) :: filename ! name of the file to read in
  real, intent(in) :: time_to_find ! time value to search for in file
  ! local
  type(fmsNetcdfFile_t) :: fileobj ! netCDF file object returned by open_file
  real, allocatable, dimension(:) :: file_times ! array of time values read from file
  integer :: dim_unlim_size, i, time_index
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  logical :: file_open_success

  time_index = 1
  dim_unlim_size = 0
  dim_unlim_name = ""
  file_open_success = .false.

  if (.not. check_if_open(fileobj)) &
    !call MOM_error(FATAL, "get_time_index_nodd: netcdf file object must be open.")
    file_open_success=fms2_open_file(fileobj, trim(filename), "read", is_restart=.false.)

  call get_unlimited_dimension_name(fileobj, dim_unlim_name)
  call get_dimension_size(fileObj, trim(dim_unlim_name), dim_unlim_size)
  ! time index will be one more than the unlimited dimension size if the time_to_find is not in the file
  if (dim_unlim_size .gt. 0) then
    time_index = dim_unlim_size+1
    allocate(file_times(dim_unlim_size))
    call read_data(fileobj,trim(dim_unlim_name), file_times)

    do i=1,dim_unlim_size
      if (ABS(file_times(i)-time_to_find) .gt. TINY(time_to_find)) then
        continue
      else
        time_index = i
        exit
      endif
    enddo
    deallocate(file_times)
  endif
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)
end function get_time_index

!> Initialize the MOM_io module
subroutine MOM_io_init(param_file)
  type(param_file_type), intent(in) :: param_file  !< structure indicating the open file to
                                                   !! parse for model parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_io" ! This module's name.

  call log_version(param_file, mdl, version)

end subroutine MOM_io_init


!> \namespace mom_io
!!
!!  This file contains a number of subroutines that manipulate
!!  NetCDF files and handle input and output of fields.  These
!!  subroutines, along with their purpose, are:
!!
!!   * create_file: create a netCDF file for and register the global axes and variables to the file
!!   * write_field: write a field to an open file generated by a previous call to create_file
!!   * MOM_read_data: read a field from a netcdf file and apply a scaling factor if specified.
!!   * MOM_read_vector : read in the components (u,v) of a vector field and apply a scaling factor to the data
!!    if specified
!!   *

end module MOM_io
