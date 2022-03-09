program test_interp

    use interp_kinds
    use interp_module
    
    implicit none

    integer(ip)            :: data_num
    integer(ip)            :: interp_num
    real(dp),dimension(100)  :: t_data, p_interp,t_interp,p_data

    data_num=100
    call sub_grid(0,data_num,0.0_dp,1.0_dp,t_data)
    call sub_grid(0,data_num,0.1_dp,1.1_dp,t_interp)
    p_data     = exp(t_data)
    interp_num = data_num

    call interp_linear(1,data_num,t_data, p_data, interp_num,t_interp, p_interp)

    write(*,'(4(F8.3))'), t_data, p_data, t_interp, p_interp


   contains

    SUBROUTINE sub_grid(n_case,len_grid,min,max,grid_out)
        ! Returns various types of grid
        ! Inputs:
        !           n_case: selects grid type
        !           len_grid : grid length
        !           min, max : smallest and largest grid node
        ! Output:
        !           grid_out
        !USE mod_types

        IMPLICIT NONE

        INTEGER, INTENT(in)     :: n_case, len_grid
        REAL (dp), INTENT(in)   :: min, max
        REAL (dp), INTENT(out), DIMENSION(len_grid) :: grid_out

        INTEGER                 :: i1
        real (dp)               :: den



        IF (n_case<0 .OR. n_case>3) THEN
        PRINT *, "Error, the first argument of grid has to be between 0 and 3"
        STOP
        ENDIF
        den = real(len_grid,dp) - 1.d0

        !PRINT *, 'len_grid', len_grid

        SELECT CASE(n_case)

        CASE(0)
            ! Uniform grid
            grid_out(len_grid) = max
            grid_out(1) = min
            DO i1 = 2, len_grid
                grid_out(i1) = grid_out(1) + REAL(i1-1,dp)/den*&
                    (grid_out(len_grid)-grid_out(1))
            END DO

        CASE(1)
            ! Exponential grid
            grid_out(len_grid) = LOG(max+1)
            grid_out(1) = LOG(min+1)
            DO i1 = 2, len_grid
                grid_out(i1) = grid_out(1) + REAL(i1-1,dp)/den*&
                    (grid_out(len_grid)-grid_out(1))
            END DO
            grid_out = EXP(grid_out)-1

        CASE(2)
            ! Double exponential grid
            grid_out(len_grid) = LOG(LOG(max+1)+1)
            grid_out(1) = LOG(LOG(min+1)+1)
            DO i1 = 2, len_grid
                grid_out(i1) = grid_out(1) + REAL(i1-1,dp)/den*&
                    (grid_out(len_grid)-grid_out(1))
            END DO
            grid_out = EXP(EXP(grid_out)-1)-1

        CASE(3)
            ! Triple exponential grid
            grid_out(len_grid) = LOG(LOG(LOG(max+1)+1)+1)
            grid_out(1) = LOG(LOG(LOG(min+1)+1)+1)
            DO i1 = 2, len_grid
                grid_out(i1) = grid_out(1) + REAL(i1-1,dp)/den*&
                    (grid_out(len_grid)-grid_out(1))
            END DO
            grid_out = EXP(EXP(EXP(grid_out)-1)-1)-1

        END SELECT

    END SUBROUTINE sub_grid   

!*****************************************************************************80
!
!! INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
!
!  Discussion:
!
!    From a space of M dimensions, we are given a sequence of
!    DATA_NUM points, which are presumed to be successive samples
!    from a curve of points P.
!
!    We are also given a parameterization of this data, that is,
!    an associated sequence of DATA_NUM values of a variable T.
!    The values of T are assumed to be strictly increasing.
!
!    Thus, we have a sequence of values P(T), where T is a scalar,
!    and each value of P is of dimension M.
!
!    We are then given INTERP_NUM values of T, for which values P
!    are to be produced, by linear interpolation of the data we are given.
!
!    Note that the user may request extrapolation.  This occurs whenever
!    a T_INTERP value is less than the minimum T_DATA or greater than the
!    maximum T_DATA.  In that case, linear extrapolation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
!    independent variable at the sample points.  The values of T_DATA
!    must be strictly increasing.
!
!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
!    dependent variables at the sample points.
!
!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
!    at which interpolation is to be done.
!
!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
!    independent variable at the interpolation points.
!
!    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
!    values of the dependent variables at the interpolation points.




end program test_interp    