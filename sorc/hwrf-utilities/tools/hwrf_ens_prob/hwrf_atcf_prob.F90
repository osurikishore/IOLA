program hwrf_atcf_prob
  use tcf_module, only: atcf, read_atcf
  use fileop_module, only: stat

  implicit none


contains

  subroutine read_namelist
    namelist /input/ atcfpath, ensmin, ensmax
    
  end subroutine read_namelist

end program hwrf_atcf_prob
