module m_stochasticSTF

! Usage
!    see lines 2-7 in main.f90
!
! More information
!    Callable function StochasticSTF(n,r) or StochasticSTF(n,r,d)
!    returns a stochastic Source Time Function (STF) of length n.
!    The floating-point number r (> 0) is the ratio of
!    two corner frequencies and determines roughness of the STF,
!    where r=1.0 results in STF with the omega-square-type spectrum,
!    and the larger abs(log(r)) yields the rougher STF.
!    Integer d (>=2), dimension of the Bessel bridge,
!    is optional; the default value is d=2.
!    The result is normalized so that sum(STF) = 1.0 holds.
!    To specity precision, see the "use" statement
!    and modify as "kd=>real32" or "kd=>real64".
!    The algorithm has been modified after Hirano(2022; 2023).
!    This code generates STFs with arbitrary lengths by using Bessel bridges,
!    while the original model has probabilistic lengths.
!
! Reference
!    Hirano, S. (2022), "Source time functions of earthquakes based on a stochastic differential equation", Scientific Reports, 12:3936, https://doi.org/10.1038/s41598-022-07873-2
!    Hirano, S. (2023), "Stochastic source time functions with double corner frequencies", AGU23 Fall Meeting, S13F-0407, https://agu.confex.com/agu/fm23/meetingapp.cgi/Paper/1299761

    use iso_fortran_env, only : kd=>real64 ! real32 for single, or real64 for double
    implicit none
    private
    public :: kd
    public :: StochasticSTF

    contains

    function StochasticSTF(n,r,d)
        implicit none
        integer,intent(in) :: n ! length of the array (>= 30 is recommended)
        real(kd),intent(in) :: r ! ratio of two corner frequencies (>= 1.0)
        integer,intent(in),optional :: d ! dimension of Brownian bridge (>= 2)
        integer :: d_in
        real(kd),allocatable :: x(:), y(:), StochasticSTF(:)
        integer :: l(2), i
        if (present(d)) then
            d_in = d
        else
            d_in = 2
        end if
        l(1) = nint(n*r/(1d0+r)) ! l(1) / l(2) = r, and
        l(2) = n - l(1) + 1      ! l(1) + l(2) = n+1
        allocate(x(n),y(n),StochasticSTF(n),source=real(0d0,kd))
        x(:l(1)) = BesselBridge(l(1),d_in)
        y(n-l(2)+1:n) = BesselBridge(l(2),d_in)
        do concurrent (i=1:n)
            StochasticSTF(i) = dot_product(x(:i),y(n+1-i:)) ! inner product (x,y)
        end do
        StochasticSTF = StochasticSTF / sum(StochasticSTF) ! normalization
    end function StochasticSTF

    function BesselBridge(n,d)
        implicit none
        integer,intent(in) :: n, d
        real(kd),allocatable :: BesselBridge(:), BB(:,:)
        BB = BrownianBridges(n+1,d)
        BesselBridge = norm2(BB(2:,:),2) ! Bessel bridge is RMS envelope of Brown bridges
    end function BesselBridge

    function BrownianBridges(n,d)
        implicit none
        integer,intent(in) :: n, d
        real(kd),allocatable :: BrownianBridges(:,:), & ! Brownian Bridges
                                              B(:,:), & ! Brownian Noises
                                             dB(:,:)    ! Gaussian Noises
        integer :: i
        allocate(BrownianBridges(n,d),B(n,d),source=real(0d0,kd))
        dB = GaussianNoises(n,d)
        do i = 2,n
            B(i,:) = B(i-1,:) + dB(i,:) ! Generating Brownian noises by numerical integration
        end do
        do concurrent (i=2:n)
            BrownianBridges(i,:) = B(i,:) - dble(i-1)/dble(n-1)*B(n,:) ! Detrending Brownian noises, which results in Brownian bridges
        end do
    end function BrownianBridges

    function GaussianNoises(n,d)
        implicit none
        real(kd),parameter :: pi = acos(-1d0)
        integer,intent(in) :: n, d
        real(kd),allocatable :: GaussianNoises(:,:), rand(:,:,:)
        allocate(rand(n,d,2))
        call random_init(repeatable=.false., image_distinct=.true.)
        call random_number(rand)
        rand = 1d0 - rand ! 0 <= rand < 1, but 0 < rand for the next procedure 
        GaussianNoises = sqrt(-2d0*log(rand(:,:,1)))*cos(2d0*pi*rand(:,:,2)) ! Box-Muller's method
    end function GaussianNoises

end module m_stochasticSTF