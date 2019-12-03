! npagane | risca lab | dec 2019 | line-line intersection calculation
! inspired from http://paulbourke.net/geometry/pointlineplane/

! define line-line intersection module that contains calculation and test
MODULE LineLineIntersection
    implicit none
    contains

    ! calculation
    FUNCTION LineLineIntersectionCalculation(A1,A2,B1,B2)
        implicit none
        ! A1 is the point: (ax1, ay1, az1)
        ! A2 is the point: (ax2, ay2, az2)
        ! A1->A2 makes the line segment A (likewise for B)
        real, intent(in), dimension(3) :: A1
        real, intent(in), dimension(3) :: A2
        real, intent(in), dimension(3) :: B1
        real, intent(in), dimension(3) :: B2
        real, parameter ::  tol = 1.0e-5 ! tolerance for cooccupancy (should be small to disallow overlap)
        real, parameter :: dist = 1 ! tolerance for collision (pseudo thickness of line)
        real, dimension(3) :: pA ! closest point on A to B
        real, dimension(3) :: pB ! closest point on B to A
        real dotA1B1B2B1, dotB2B1A2A1, dotA1B1A2A1, dotB2B1B2B1, dotA2A1A2A1
        real muA, muB
        real, dimension(3) :: vecA, vecB, tA2, tB1, tB2
        logical LineLineIntersectionCalculation 

        ! default value (i.e. no intersection)
        LineLineIntersectionCalculation = .FALSE.

        ! check for overlap of points
        if (ALL(ABS(A1-B1) <= tol) .OR. ALL(ABS(A1-B2) <= tol) .OR. ALL(ABS(A2-B1) <= tol) .OR. ALL(ABS(A2-B2) <= tol)) then
            print*, "collision, point overlap"
            LineLineIntersectionCalculation = .TRUE.
            return
        endif

        ! check if lines are parallel
        vecA = (A2-A1)/NORM2(A2-A1)
        vecB = (B2-B1)/NORM2(B2-B1)
        if ( ALL(ABS(vecA) - ABS(vecB) <= tol) ) then
            ! try to find overlap
            tA2 = (A2-A1)/vecA
            tB1 = (B1-A1)/vecA
            tB2 = (B2-A1)/vecA
            ! take just the first component
            if (( (tA2(1)*tB1(1) > 0) .OR. (tA2(1)*tB2(1) > 0) ) &
              .AND. ( (abs(tA2(1)) >= abs(tB1(1))) .OR. (abs(tA2(1)) >= abs(tB2(1))) )) then
                print*, "collision, parallel overlap"
                LineLineIntersectionCalculation = .TRUE.
                return
            else
                return ! quit early
            endif
        endif

        ! find shortest line between A and B
        dotA1B1B2B1 = dot_product(A1-B1, B2-B1)
        dotB2B1A2A1 = dot_product(B2-B1, A2-A1)
        dotA1B1A2A1 = dot_product(A1-B1, A2-A1)
        dotB2B1B2B1 = dot_product(B2-B1, B2-B1)
        dotA2A1A2A1 = dot_product(A2-A1, A2-A1)
        muA = (dotA1B1B2B1*dotB2B1A2A1 - dotA1B1A2A1*dotB2B1B2B1) / (dotA2A1A2A1*dotB2B1B2B1 - dotB2B1A2A1*dotB2B1A2A1)
        muB = (dotA1B1B2B1 + mua*dotB2B1A2A1) / dotB2B1B2B1
        pA = A1 + muA * (A2-A1)
        pB = B1 + muB * (B2-B1)
        ! check if dist <= tol
        if ((sqrt(dot_product(pA-pB, pA-pB)) <= dist)  .AND. (abs(muA) <= 1) .AND. (abs(muB) <= 1)) then 
            print*, "collision, intersect"
            LineLineIntersectionCalculation = .TRUE.
            return 
        endif
        return 

    END FUNCTION LineLineIntersectionCalculation

    ! test
    SUBROUTINE LineLineIntersectionTest
        implicit none
        real, dimension(3) :: A1
        real, dimension(3) :: A2 
        real, dimension(3) :: B1 
        real, dimension(3) :: B2
        logical val
        integer, parameter :: multFactor = 10
        integer, parameter :: lineFactor = 5
        real, dimension(6,100) :: lines
        integer :: numLines = 1
        integer sumLines
        integer i,j
        
        open(1, file="test/lines.dat", status='new')
        ! set random line for initialization
        lines(1:3,1) = (/multFactor*RAND(), multFactor*RAND(), multFactor*RAND()/)
        lines(4:6,1) = lines(1:3,1) + (/lineFactor*RAND(), lineFactor*RAND(), lineFactor*RAND()/)
        write(1,*) lines(1,1), lines(2,1), lines(3,1), &
          lines(4,1), lines(5,1), lines(6,1)
        do i = 1,100 ! MAKE SURE THIS MATCHES LINES DIMENSION
            A1(1) = RAND()*multFactor
            A1(2) = RAND()*multFactor
            A1(3) = RAND()*multFactor
            A2 = A1 + (/lineFactor*RAND(), lineFactor*RAND(), lineFactor*RAND()/)
            sumLines = 0
            do j = 1,numLines
                B1 = lines(1:3, j)
                B2 = lines(4:6, j)
                val = LineLineIntersectionCalculation(A1, A2, B1, B2) 
                if (val .NEQV. .TRUE.) then
                    sumLines = sumLines + 1
                endif
            enddo
            if (sumLines == numLines) then
                numLines = numLines + 1
                lines(1:3,numLines) = A1
                lines(4:6,numLines) = A2
                write(1,*) lines(1,numLines), lines(2,numLines), lines(3,numLines), &
                  lines(4,numLines), lines(5,numLines), lines(6,numLines)
            endif
        enddo
        close(1)

    END SUBROUTINE LineLineIntersectionTest

END MODULE LineLineIntersection

! test line-line intersection module
PROGRAM LineLineTest 
    use LineLineIntersection, only: LineLineIntersectionTest
    implicit none

    call LineLineIntersectionTest()

END PROGRAM
