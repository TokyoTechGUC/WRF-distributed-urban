!**********************************************************************************
!     netCDF handling program
!     Originally written by M. Nakayoshi 07.01.01
!     This code modifies the geo_em_d0?.nc file as to include 2nd dominant land use category
!     and land use fraction ratio between 1st and 2nd dominant land use category.
!  (1)   If dominant category is urban, then simply search 2nd dominant category.
!  (2)   Even if urban category is neither nor 2nd dominant category,
!        urban is set to be 2nd dominant category.
!**********************************************************************************
! Compile using gfortran

program case_complete
   use netcdf
   implicit none

   integer, parameter :: water = 17, urban = 13 ! MODIS

   call main()

contains

   subroutine main()
      integer :: ncid, status, domain
      character(len=2) :: num
      character(len=33) :: file_name

      domain = 0
      do while (.true.)
         domain = domain + 1

         write (num, '(i2.2)') domain
         file_name = "geo_em.d"//num//".nc"
         status = nf90_open(file_name, nf90_write, ncid) ! get input file ID
         if (status /= nf90_noerr) exit
         print *, "Processing ", file_name

         call process_dataset(ncid)

         status = nf90_close(ncid)
         call check_status(status)
      end do

      print *, "Finish!"
   end subroutine

   subroutine process_dataset(ncid)
      integer, intent(in) :: ncid

      real, dimension(:, :, :, :), allocatable :: landusef
      real, dimension(:, :, :), allocatable :: &
         lu_index, second, second_fraction, &
         rural_urb_ratio, landmask, dominant_category, &
         mod_frc_urb2d, &
         disp ! displacement height
      integer :: lx, ly, lz, lt, i, j
      integer :: dimid(4)

      ! Get dimensional information
      print *, "  Get dimensional information"
      call dim_info(ncid, 'west_east', dimid(1), lx)
      call dim_info(ncid, 'south_north', dimid(2), ly)
      call dim_info(ncid, 'land_cat', dimid(3), lz)
      call dim_info(ncid, 'Time', dimid(4), lt)

      print *, "  Read dim_info"
      !print *, (dimid(k), k=1, 7)
      !print *, "  lx = ", lx, "ly = ", ly, "lz = ", lz
      !print *, "  lxs =", lxs, "lys =", lys, "lzs =", lzs, "lt =", lt

      ! Allocate
      allocate ( &
         lu_index(lx, ly, lt), &
         landusef(lx, ly, lz, lt), &
         second_fraction(lx, ly, lt), &
         second(lx, ly, lt), &
         rural_urb_ratio(lx, ly, lt), &
         landmask(lx, ly, lt), &
         dominant_category(lx, ly, lt), &
         mod_frc_urb2d(lx, ly, lt), &
         disp(lx, ly, lt) &
         )

      ! Get 3-dimensional real values
      call fill_3d_var_real(ncid, 'LU_INDEX', lu_index)
      call fill_3d_var_real(ncid, 'LANDMASK', landmask)
      call fill_3d_var_real(ncid, 'DISP', disp)
      print *, "  Read 3D variables"

      ! Get 4-dimensional real values
      call fill_4d_var_real(ncid, 'LANDUSEF', landusef)
      print *, "  Read 4D variables"

      ! LU_INDEX is defined as the *non-water* dominant category for land masked cell
      ! Hence, LU_INDEX and dominant category might differ
      do i = 1, lx
         do j = 1, ly
            call calculate( &
               landusef(i, j, :, lt), &
               nint(lu_index(i, j, lt)), &
               nint(landmask(i, j, lt)), &
               dominant_category(i, j, lt), &
               second(i, j, lt), &
               second_fraction(i, j, lt), &
               rural_urb_ratio(i, j, lt), &
               mod_frc_urb2d(i, j, lt), &
               disp(i, j, lt) &
               )
         end do
      end do

      print *, "  Finish of dominant categories searching"

      print *, "  Writing (new) variables"
      call output_3d_real(ncid, dimid((/1, 2, 4/)), &
                          'MOD_FRC_URB2D', mod_frc_urb2d)
      call output_3d_real(ncid, dimid((/1, 2, 4/)), &
                          'DOMINANT_CATEGORY', dominant_category)
      call output_3d_real(ncid, dimid((/1, 2, 4/)), &
                          'SECOND', second)
      call output_3d_real(ncid, dimid((/1, 2, 4/)), &
                          'SECOND_FRACTION', second_fraction)
      call output_3d_real(ncid, dimid((/1, 2, 4/)), &
                          'RURAL_URB_RATIO', rural_urb_ratio)
      call output_3d_real(ncid, dimid((/1, 2, 4/)), &
                          'DISP', disp)
   end subroutine

   subroutine calculate(landusef, lu_index, landmask, dominant_category, &
                        second, second_fraction, rural_urb_ratio, mod_frc_urb2d, &
                        disp)
      real, dimension(:), intent(in) :: landusef
      integer, intent(in) :: lu_index, landmask
      real, intent(out) :: dominant_category, second, second_fraction, &
                           rural_urb_ratio, mod_frc_urb2d
      real, intent(inout) :: disp

      ! First and second dominant category
      integer :: dominant_cat, second_dominant_cat
      ! First and second dominant category fraction
      real :: dominant_frac, second_dominant_frac
      integer :: k

      dominant_frac = -1
      second_dominant_frac = -1
      dominant_cat = -1
      second_dominant_cat = -1

      do k = lbound(landusef, 1), ubound(landusef, 1)
         if (landusef(k) > dominant_frac) then
            second_dominant_frac = dominant_frac
            second_dominant_cat = dominant_cat

            dominant_frac = landusef(k)
            dominant_cat = k
         else if (landusef(k) > second_dominant_frac) then
            second_dominant_frac = landusef(k)
            second_dominant_cat = k
         end if
      end do

      if (dominant_cat == -1 .or. second_dominant_cat == -1) then
         stop "Something was wrong with the computation"
      end if

      dominant_category = dominant_cat

      if (lu_index == urban) then
         if (dominant_cat /= urban) then
            ! If the dominant category is not urban,
            ! it must be water in a land-masked cell
            if (dominant_cat /= water .or. landmask /= 1) then
               stop "Invalid assumption"
            end if
            call swap_i(dominant_cat, second_dominant_cat)
            call swap_r(dominant_frac, second_dominant_frac)
         end if

         if (dominant_cat /= urban) stop "Invalid assumption"

         second = second_dominant_cat
         second_fraction = second_dominant_frac
         rural_urb_ratio = second_dominant_frac / dominant_frac
      else if (landusef(urban) >= 0.01) then
         ! If there is at least 1% urban, set urban as the 2nd dominant
         second = urban
         second_fraction = landusef(urban)
         rural_urb_ratio = landusef(lu_index) / second_fraction
      else
         second = 20
         second_fraction = -999 !error flag
         rural_urb_ratio = -999 !error flag
      end if

      if (lu_index == second) then
         print *, "1st and 2nd dominant are same", &
            second, rural_urb_ratio
      end if

      if (rural_urb_ratio >= 0.0) then
         mod_frc_urb2d = chomp(1.0 / (1 + rural_urb_ratio), 0.0, 0.9)
      else
         mod_frc_urb2d = 0.0
      end if

      if (mod_frc_urb2d <= 0.0) then
         disp = 0.0
      end if

   end subroutine

   subroutine output_3d_real(ncid, dimids, var_name, var_3d, varid_out)
      integer, intent(in) :: ncid
      integer, dimension(:) :: dimids
      character(*), intent(in) :: var_name
      real, dimension(:, :, :) :: var_3d
      integer, optional, intent(out) :: varid_out

      integer :: varid
      integer :: status

      status = nf90_inq_varid(ncid, var_name, varid)
      if (status /= nf90_noerr) then
         ! Define if not exist
         print *, "  Creating a new variable ", var_name

         status = nf90_redef(ncid)
         call check_status(status)

         status = nf90_def_var(ncid, var_name, nf90_float, dimids, varid)
         call check_status(status)

         call check_status(nf90_put_att(ncid, varid, 'FieldType', 104))
         call check_status(nf90_put_att(ncid, varid, 'MemoryOrder', 'XY'))
         call check_status(nf90_put_att(ncid, varid, 'stagger', 'M'))
         call check_status(nf90_put_att(ncid, varid, 'units', '-'))
         call check_status(nf90_put_att(ncid, varid, 'description', '-'))
         call check_status(nf90_put_att(ncid, varid, 'sr_x', 1))
         call check_status(nf90_put_att(ncid, varid, 'sr_y', 1))

         status = nf90_enddef(ncid)
         call check_status(status)
      end if

      ! Write or overwrite
      !print *, "ncid =", ncid, "varid =", varid
      status = nf90_put_var(ncid, varid, var_3d)
      call check_status(status)

      if (present(varid_out)) varid_out = varid
   end subroutine output_3d_real

!--------------------------------------------------------------
!         call dim_info(ncid,4,    'time'    ,dimid,lt  )
!--------------------------------------------------------------
   subroutine dim_info(ncid, dim_name, dimid, nsize)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: dim_name
      integer, intent(out) :: dimid, nsize

      integer :: status

      status = nf90_inq_dimid(ncid, dim_name, dimid)
      call check_status(status)

      status = nf90_inquire_dimension(ncid, dimid, len=nsize)
      call check_status(status)
   end subroutine dim_info

!-------------------------------------------------------------
!         call fill_3d_var_real(ncid, 'U10', U10, varid)
!----------------------------------------------------------------------------
   subroutine fill_3d_var_real(ncid, var_name, values, varid_out)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:, :, :) :: values
      integer, optional, intent(out) :: varid_out

      integer :: varid, status

      status = nf90_inq_varid(ncid, var_name, varid)
      call check_status(status)
      !write (*, *) 'varid for ', var_name, '=', varid
      status = nf90_get_var(ncid, varid, values)
      call check_status(status)

      if (present(varid_out)) varid_out = varid
   end subroutine fill_3d_var_real

!-------------------------------------------------------------
!         call fill_4d_var_real(ncid, 'QCLOUD', U, varid)
!----------------------------------------------------------------------------
   subroutine fill_4d_var_real(ncid, var_name, values, varid_out)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      real, dimension(:, :, :, :), intent(out) :: values
      integer, optional, intent(out) :: varid_out

      integer :: varid, status

      status = nf90_inq_varid(ncid, var_name, varid)
      call check_status(status)
      !write (*, *) 'varid for ', var_name, '=', varid
      status = nf90_get_var(ncid, varid, values)
      call check_status(status)

      if (present(varid_out)) varid_out = varid
   end subroutine fill_4d_var_real

   subroutine check_status(status)
      integer, intent(in) :: status
      if (status /= nf90_noerr) then
         print *, 'error_number =', status
         print *, 'error message =', nf90_strerror(status)
         stop "Non-zero status"
      end if
   end subroutine check_status

   subroutine swap_i(a, b)
      integer, intent(inout) :: a, b
      integer :: tmp
      tmp = a
      a = b
      b = tmp
   end subroutine

   subroutine swap_r(a, b)
      real, intent(inout) :: a, b
      real :: tmp
      tmp = a
      a = b
      b = tmp
   end subroutine

   function chomp(x, lower_bound, upper_bound)
      real, intent(in) :: x, lower_bound, upper_bound
      real :: chomp

      chomp = x
      if (chomp < lower_bound) chomp = lower_bound
      if (chomp > upper_bound) chomp = upper_bound
   end function
end program
