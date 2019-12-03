! npagane | risca lab | dec 2019 | sphere-sphere intersection calculation

! define sphere-sphere intersection module that contains calculation and test
MODULE SphereSphereIntersection
    implicit none
    contains

    ! calculation
    FUNCTION SphereSphereIntersectionCalculation(A1,ra,B1,rb)
        implicit none
        ! A1 is the center of the sphere (x,y,z) of radius ra
        ! B1 is the center of the sphere (x,y,z) of radius rb
        real, intent(in), dimension(3) :: A1
         real, intent(in) :: ra
        real, intent(in), dimension(3) :: B1
        real, intent(in) :: rb
        real, parameter :: dist = 1.0e-5 ! tolerance for collision (should be small)
        logical SphereSphereIntersectionCalculation

        ! default value (i.e. no intersection)
        SphereSphereIntersectionCalculation = .FALSE.

        ! see if radii overalp
        if ( sqrt(dot_product(A1-B1, A1-B1)) - (ra+rb) < dist ) then
            SphereSphereIntersectionCalculation = .TRUE.
            !print*, "sphere collision"
            return
        endif
        return 

    END FUNCTION SphereSphereIntersectionCalculation

    ! testing routines !
    SUBROUTINE SphereSphereIntersectionTestAinB
        implicit none
        real, dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real :: ra = 2
        real, dimension(3) :: B1 = (/0.0,0.0,0.0/)
        real :: rb = 20.0
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestAinB"
            stop
        endif
   
    END SUBROUTINE SphereSphereIntersectionTestAinB

    SUBROUTINE SphereSphereIntersectionTestBinA
        implicit none
        real, dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real :: ra = 20
        real, dimension(3) :: B1 = (/0.0,0.0,0.0/)
        real :: rb = 2
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestBinA"
            stop
        endif
   
    END SUBROUTINE SphereSphereIntersectionTestBinA

    SUBROUTINE SphereSphereIntersectionTestTangent
        implicit none
        real, dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real :: ra = 1
        real, dimension(3) :: B1 = (/3.0,0.0,0.0/)
        real :: rb = 1
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestTangent"
            stop
        endif

    END SUBROUTINE SphereSphereIntersectionTestTangent

    SUBROUTINE SphereSphereIntersectionTestOverlap
        implicit none
        real, dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real :: ra = 1
        real, dimension(3) :: B1 = (/3.0,0.0,0.0/)
        real :: rb = 1
        logical val
    
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestOverlap"
            stop 
        endif

    END SUBROUTINE SphereSphereIntersectionTestOverlap

    SUBROUTINE SphereSphereIntersectionTestNoOverlap
        implicit none
        real, dimension(3) :: A1 = (/1.0,0.0,0.0/)
        real :: ra = 1
        real, dimension(3) :: B1 = (/3.0,0.0,0.0/)
        real :: rb = 0.99
        logical val
   
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestNoOverlap"
            stop
        endif

    END SUBROUTINE SphereSphereIntersectionTestNoOverlap

END MODULE SphereSphereIntersection

! test sphere-line intersection module
PROGRAM SphereSphereTest 
    use SphereSphereIntersection, only: SphereSphereIntersectionTestAinB, SphereSphereIntersectionTestBinA, &
      SphereSphereIntersectionTestTangent, SphereSphereIntersectionTestOverlap, SphereSphereIntersectionTestNoOverlap
    implicit none

    call SphereSphereIntersectionTestAinB()
    call SphereSphereIntersectionTestBinA()
    call SphereSphereIntersectionTestTangent()
    call SphereSphereIntersectionTestOverlap()
    call SphereSphereIntersectionTestNoOverlap()

    print*, "SUCCESS: successful completion of all 5 sphere-sphere collision unit tests"

END PROGRAM
