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
            print*, "sphere collision"
            return
        endif
        return 

    END FUNCTION SphereSphereIntersectionCalculation

    ! test
    SUBROUTINE SphereSphereIntersectionTest
        implicit none
        real, dimension(3) :: A1 = (/5.0,0.0,0.0/)
        real :: ra = 2
        real, dimension(3) :: B1 = (/0.0,0.0,0.0/)
        real :: rb = 20.0
        logical val
        !real, dimension(6,100) :: lines
        !integer :: numLines = 1
        !integer sumLines
        !integer i,j
     
        val = SphereSphereIntersectionCalculation(A1, ra, B1, rb)
        print*, val 
   
        !open(1, file="test/spherelines.dat", status='new')
        ! set random line for initialization
        !lines(1:3,1) = (/multFactor*RAND(), multFactor*RAND(), multFactor*RAND()/)
        !lines(4:6,1) = lines(1:3,1) + (/lineFactor*RAND(), lineFactor*RAND(), lineFactor*RAND()/)
        !write(1,*) lines(1,1), lines(2,1), lines(3,1), &
        !  lines(4,1), lines(5,1), lines(6,1)
        !do i = 1,100 ! MAKE SURE THIS MATCHES LINES DIMENSION
        !    A1(1) = RAND()*multFactor
        !    A1(2) = RAND()*multFactor
        !    A1(3) = RAND()*multFactor
        !    A2 = A1 + (/lineFactor*RAND(), lineFactor*RAND(), lineFactor*RAND()/)
        !    sumLines = 0
        !    do j = 1,numLines
        !        B1 = lines(1:3, j)
        !        B2 = lines(4:6, j)
        !        val = SphereLineIntersectionCalculation(A1, A2, B1, r) 
        !        if (val .NEQV. .TRUE.) then
        !            sumLines = sumLines + 1
        !        endif
        !    enddo
        !    if (sumLines == numLines) then
        !        numLines = numLines + 1
        !        lines(1:3,numLines) = A1
        !        lines(4:6,numLines) = A2
        !        write(1,*) lines(1,numLines), lines(2,numLines), lines(3,numLines), &
        !          lines(4,numLines), lines(5,numLines), lines(6,numLines)
        !    endif
        !enddo
        !close(1)

    END SUBROUTINE SphereSphereIntersectionTest

END MODULE SphereSphereIntersection

! test sphere-line intersection module
PROGRAM SphereSphereTest 
    use SphereSphereIntersection, only: SphereSphereIntersectionTest
    implicit none

    call SphereSphereIntersectionTest()

END PROGRAM
