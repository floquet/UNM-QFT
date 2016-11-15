!     open namelist
      open  ( unit = io_unit_nml, file = file_name_param_base, delim = 'apostrophe', iostat = io_status )      
      if ( io_status /= 0 ) then                               ! can't open file
        write ( *,  * )
        write ( *,  * ) 'unable to open file ', file_name_param_base
        write ( *,  * ) 'trying to write namelist run_param_base'
        write ( *,  * ) 'iostat  = ',  io_status
        write ( *,  * ) 'io unit = ',  io_unit_nml
        stop  'stop on error during file open'                    ! stop program
      end if

!     write namelist
      write ( unit = io_unit_nml, nml = run_params, iostat = io_status )
      if ( io_status /= 0 ) then                              ! can't write file
        write ( *,  * )
        write ( *,  * ) 'unable to write file ', file_name_param_base
        write ( *,  * ) 'trying to write namelist run_param_base'
        write ( *,  * ) 'iostat  = ',  io_status
        write ( *,  * ) 'io unit = ',  io_unit_nml
        stop  'stop on error during file write'                   ! stop program
      end if

      write ( *, * ) 'successful write to ', file_name_param_base
      
!     close file
      close ( unit = io_unit_nml, iostat = io_status )
      if ( io_status /= 0 ) then                              ! can't write file
        write ( *,  * )
        write ( *,  * ) 'unable to close file ', file_name_param_base
        write ( *,  * ) 'trying to write namelist run_param_base'
        write ( *,  * ) 'iostat  = ',  io_status
        write ( *,  * ) 'io unit = ',  io_unit_nml
        stop  'stop on error during file close'                   ! stop program
      end if