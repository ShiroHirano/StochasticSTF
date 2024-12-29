program main
    use m_stochasticSTF, only : kd, StochasticSTF, output_svg_path
    implicit none
    real(kd), allocatable :: STF(:)   ! Source Time Function to be obtained.
    integer,  parameter   :: n = 1000 ! Length of STF (> 30 is recommended).
    real(kd), parameter   :: r = 1e0  ! (Optional. Ratio of two corner frequencies;
                                      ! default value is 1.0, 
                                      ! resulting in omega-square model).
!    integer,  parameter   :: d = 2    ! (Optional. Dimension of Bessel bridege;
                                      ! default value is 2.
                                      ! Larger it is, the smoother STF is).
    STF = StochasticSTF(n,r)
!    STF = StochasticSTF(n,r,d)
    call output_svg_path(STF) ! The result is output as STF.svg.
                              ! Ignore the first and last line
                              ! of the file to get numerics.
end program main
