!*******************************************************************************
!>
!  Numeric kinds.

    module interp_kinds

        use iso_fortran_env, only: real64,int32
    
        private
    
        integer,parameter,public :: dp = real64  !! default real kind
        integer,parameter,public :: ip = int32   !! default integer kind
    
    end module
    !****************************************************************************