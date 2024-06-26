program main
    use m_stochasticSTF, only : kd, StochasticSTF
    implicit none
    real(kd),allocatable :: STF(:) ! Source Time Function to be obtained
    integer,parameter :: n = 1000  ! length of STF (> 30 is recommended)
    real(kd),parameter :: r = 1e0  ! roughness of STF (1.0 for omega-seuqre)
    STF = StochasticSTF(n,r)
    open(00,file="STF.txt",action="write",status="replace")
        write(00,'(1G0)') STF
    close(00)
end program main
