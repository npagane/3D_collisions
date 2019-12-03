! npagane | risca lab | dec 2019 | sphere-line intersection calculation
! inspired from http://paulbourke.net/geometry/circlesphere/

! define sphere-line intersection module that contains calculation and test
MODULE SphereLineIntersection
    implicit none
    contains

    ! calculation
    FUNCTION SphereLineIntersectionCalculation(A1,A2,B1,r)
        implicit none
        ! A1 is the point: (ax1, ay1, az1)
        ! A2 is the point: (ax2, ay2, az2)
        ! B1 is the center of the sphere (x,y,z) of radius r
        real, intent(in), dimension(3) :: A1
        real, intent(in), dimension(3) :: A2
        real, intent(in), dimension(3) :: B1
        real, intent(in) :: r
        real :: a, b, c, discr, um, up ! quadratic variables
        logical SphereLineIntersectionCalculation

        ! default value (i.e. no intersection)
        SphereLineIntersectionCalculation = .FALSE.

        ! calculate discriminant
        a = dot_product(A2-A1, A2-A1)
        b = 2*dot_product(A1-B1, A1-B1)
        c = dot_product(B1, B1) + dot_product(A1, A1) - 2*dot_product(B1, A1) - r*r
        discr = b*b - 4*a*c
        if ( discr < 0 ) then 
            return
        endif
        up = (-b + sqrt(discr))/(2*a)
        um = (-b - sqrt(discr))/(2*a)
        if ((um > 1 .AND. up < 0) .OR. (um < 0 .AND. up > 1)) then 
            SphereLineIntersectionCalculation = .TRUE.
            print*, "line in sphere"
            return
        else if ((um >= 0 .AND. um <= 1) .OR. (up >= 0 .AND. up <= 1)) then 
            SphereLineIntersectionCalculation = .TRUE.
            print*, "one or more intersections"
            return
        endif
        return 

    END FUNCTION SphereLineIntersectionCalculation

    ! test
    SUBROUTINE SphereLineIntersectionTest
        implicit none
        real, dimension(3) :: A1 = (/-1.1,0.0,0.0/)
        real, dimension(3) :: A2 = (/-2.0,0.0,0.0/)
        real, dimension(3) :: B1 = (/0.0,0.0,0.0/)
        real :: r = 20.0
        logical val
        !real, dimension(6,100) :: lines
        !integer :: numLines = 1
        !integer sumLines
        !integer i,j
     
        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
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

    END SUBROUTINE SphereLineIntersectionTest

END MODULE SphereLineIntersection

! test sphere-line intersection module
PROGRAM SphereLineTest 
    use SphereLineIntersection, only: SphereLineIntersectionTest
    implicit none

    call SphereLineIntersectionTest()

END PROGRAM


