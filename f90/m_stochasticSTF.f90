module m_stochasticSTF

! Usage
!    see main.f90
!    To specity precision, see the "use" statement of this module
!    and modify as "kd=>real32" or "kd=>real64".

    use iso_fortran_env, only : kd=>real64 ! real32 for single, or real64 for double
    implicit none
    private
    public :: kd
    public :: StochasticSTF
    public :: output_svg_path

    contains

    function StochasticSTF(n,r_in,d_in)
        implicit none
        integer,intent(in) :: n ! length of the array (>= 30 is recommended)
        real(kd),intent(in),optional :: r_in ! ratio of two corner frequencies (>= 1.0)
        integer,intent(in),optional :: d_in ! dimension of Brownian bridge (>= 2)
        real(kd) :: r
        integer :: d
        real(kd),allocatable :: x(:), y(:), StochasticSTF(:)
        integer :: l(2), i
        r = merge(r_in, real(1d0,kd), present(r_in)) ! Ternary operation
        d = merge(d_in, 2, present(d_in))            ! Ternary operation
        l(1) = nint(n*r/(1d0+r)) ! l(1) / l(2) = r, and
        l(2) = n - l(1) + 1      ! l(1) + l(2) = n+1
        allocate(x(n),y(n),StochasticSTF(n),source=real(0d0,kd))
        x(:l(1)) = BesselBridge(l(1),d)
        y(n-l(2)+1:n) = BesselBridge(l(2),d)
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

    subroutine output_svg_path(y)
        real(kd), intent(in) :: y(:)
        integer :: i, n
        n = size(y)
        open(00,file="STF.svg")
        write(00,'(1a)',advance="no") '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 '
        write(00,'(1G0," ",1G0)',advance="no") n, n/2
        write(00,'(1a)') '"><path d="M'
        do i = 1, n
            write(00,'(1G0," ",1es13.7)') i-1, y(i)
        end do
        write(00, '(1a)',advance="no") '" transform="translate(0,'
        write(00,'(1G0)',advance="no") floor(n*0.5)
        write(00,'(1a)',advance="no") ') scale(1, '
        write(00,'(1G0)',advance="no") -n/2*n/3
        write(00,'(1a)',advance="no") ')" style="fill:#0000ff;fill-opacity:0.25;" /></svg>'
        close(00)
    end subroutine output_svg_path

end module m_stochasticSTF
