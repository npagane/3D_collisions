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
        real, parameter :: dist = 1.0e-3 ! tolerance for collision (pseudo thickness of line)
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
            !print*, "collision, point overlap"
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
            ! ensure all the components are the same
            if ( tA2(1)==tA2(2) .AND. tA2(2)==tA2(3) .AND. tB1(1)==tB1(2) .AND. tB1(2)==tB1(3) &
              .AND. tB2(1)==tB2(2) .AND. tB2(2)==tB2(3) ) then
                ! take just the first component
                if (( (tA2(1)*tB1(1) > 0) .OR. (tA2(1)*tB2(1) > 0) ) &
                  .AND. ( (abs(tA2(1)) >= abs(tB1(1))) .OR. (abs(tA2(1)) >= abs(tB2(1))) )) then
                    !print*, "collision, parallel overlap"
                    LineLineIntersectionCalculation = .TRUE.
                    return
                else
                    return ! quit early
                endif
            else
                !print*, "no collision, parallel but not overlap"
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
            !print*, "collision, intersect"
            LineLineIntersectionCalculation = .TRUE.
            return 
        endif
        return 

    END FUNCTION LineLineIntersectionCalculation

    ! testing routines !
    SUBROUTINE LineLineIntersectionTestOverlapA1B1
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/) 
        real, dimension(3) :: A2 = (/1,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real, dimension(3) :: B2 = (/-1,-1,-1/)
        logical :: val 
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then 
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA1B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA1B1
  
    SUBROUTINE LineLineIntersectionTestOverlapA1B2
        implicit none
        real, dimension(3) :: A1 = (/-11,-11,-11/)
        real, dimension(3) :: A2 = (/-10,-10,-10/)
        real, dimension(3) :: B1 = (/1,1,1/)
        real, dimension(3) :: B2 = (/-10,-10,-10/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then 
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA1B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA1B2

    SUBROUTINE LineLineIntersectionTestOverlapA2B1
        implicit none
        real, dimension(3) :: A1 = (/1,0,0/)
        real, dimension(3) :: A2 = (/0,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real, dimension(3) :: B2 = (/0,1,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA2B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA2B1

    SUBROUTINE LineLineIntersectionTestOverlapA2B2
        implicit none
        real, dimension(3) :: A1 = (/1,0,0/)
        real, dimension(3) :: A2 = (/0.0001,0.00002,0.00004/)
        real, dimension(3) :: B1 = (/0,0,1/)
        real, dimension(3) :: B2 = (/0.0001,0.00002,0.00004/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestOverlapA2B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestOverlapA2B2

    SUBROUTINE LineLineIntersectionTestSameLine
        implicit none
        real, dimension(3) :: A1 = (/1,1,1/)
        real, dimension(3) :: A2 = (/0,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real, dimension(3) :: B2 = (/1,1,1/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestSameLine"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestSameLine

    SUBROUTINE LineLineIntersectionTestParallelOverlapA1B1
        implicit none
        real, dimension(3) :: A1 = (/2,2,2/)
        real, dimension(3) :: A2 = (/0,0,0/)
        real, dimension(3) :: B1 = (/1,1,1/)
        real, dimension(3) :: B2 = (/3,3,3/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA1B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA1B1

    SUBROUTINE LineLineIntersectionTestParallelOverlapA1B2
        implicit none
        real, dimension(3) :: A1 = (/2,2,2/)
        real, dimension(3) :: A2 = (/-100,-100,-100/)
        real, dimension(3) :: B1 = (/100,100,100/)
        real, dimension(3) :: B2 = (/1.5,1.5,1.5/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA1B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA1B2

    SUBROUTINE LineLineIntersectionTestParallelOverlapA2B1
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/-100,-100,-100/)
        real, dimension(3) :: B1 = (/-99,-99,-99/)
        real, dimension(3) :: B2 = (/-1001,-1001,-1001/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA2B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA2B1

    SUBROUTINE LineLineIntersectionTestParallelOverlapA2B2
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/1,1,1/)
        real, dimension(3) :: B1 = (/2,2,2/)
        real, dimension(3) :: B2 = (/0.999, 0.999, 0.999/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelOverlapA2B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelOverlapA2B2

    SUBROUTINE LineLineIntersectionTestParallelA1B1
        implicit none
        real, dimension(3) :: A1 = (/2,2,2/)
        real, dimension(3) :: A2 = (/0,0,0/)
        real, dimension(3) :: B1 = (/2.1,2.1,2.1/)
        real, dimension(3) :: B2 = (/3,3,3/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA1B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA1B1

    SUBROUTINE LineLineIntersectionTestParallelA1B2
        implicit none
        real, dimension(3) :: A1 = (/1,1,1/)
        real, dimension(3) :: A2 = (/-100,-100,-100/)
        real, dimension(3) :: B1 = (/100,100,100/)
        real, dimension(3) :: B2 = (/1.5,1.5,1.5/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA1B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA1B2

    SUBROUTINE LineLineIntersectionTestParallelA2B1
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/-100,-100,-100/)
        real, dimension(3) :: B1 = (/-100.0001,-100.0001,-100.0001/)
        real, dimension(3) :: B2 = (/-1001,-1001,-1001/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA2B1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA2B1

    SUBROUTINE LineLineIntersectionTestParallelA2B2
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/0.999, 0.999, 0.999/)
        real, dimension(3) :: B1 = (/2,2,2/)
        real, dimension(3) :: B2 = (/1.001,1.001,1.001/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestParallelA2B2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestParallelA2B2

    SUBROUTINE LineLineIntersectionTestIntersectA1
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/1,0,0/)
        real, dimension(3) :: B1 = (/0,-1,0/)
        real, dimension(3) :: B2 = (/0,1,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectA1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectA1

    SUBROUTINE LineLineIntersectionTestIntersectA2
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/0,0,100/)
        real, dimension(3) :: B1 = (/-200,0,100/)
        real, dimension(3) :: B2 = (/200,0,100/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectA2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectA2

    SUBROUTINE LineLineIntersectionTestIntersectB1
        implicit none
        real, dimension(3) :: A1 = (/-1,0,1/)
        real, dimension(3) :: A2 = (/1,0,1/)
        real, dimension(3) :: B1 = (/0,0,1/)
        real, dimension(3) :: B2 = (/0,0,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectB1"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectB1

    SUBROUTINE LineLineIntersectionTestIntersectB2
        implicit none
        real, dimension(3) :: A1 = (/0,-1,0/)
        real, dimension(3) :: A2 = (/0,1,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real, dimension(3) :: B2 = (/0,1,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectB2"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectB2

    SUBROUTINE LineLineIntersectionTestIntersectMiddle
        implicit none
        real, dimension(3) :: A1 = (/0,-2,0/)
        real, dimension(3) :: A2 = (/0,2,0/)
        real, dimension(3) :: B1 = (/-2,0,0/)
        real, dimension(3) :: B2 = (/2,0,0/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectMiddle"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectMiddle

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollide
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/0,0,5/)
        real, dimension(3) :: B1 = (/-1,0,2/)
        real, dimension(3) :: B2 = (/1,0,2/)
        logical :: val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollide"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollide


    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollide
        implicit none
        real, dimension(3) :: A1 = (/-1,0,5/)
        real, dimension(3) :: A2 = (/1,0,5/)
        real, dimension(3) :: B1 = (/-1,0,2/)
        real, dimension(3) :: B2 = (/1,0,2/)
        logical :: val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollide"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollide


END MODULE LineLineIntersection

! test line-line intersection module
PROGRAM LineLineTest 
    use LineLineIntersection, only: LineLineIntersectionTestOverlapA1B1, LineLineIntersectionTestOverlapA1B2, &
      LineLineIntersectionTestOverlapA2B1, LineLineIntersectionTestOverlapA2B2, LineLineIntersectionTestSameLine, &
      LineLineIntersectionTestParallelOverlapA1B1, LineLineIntersectionTestParallelOverlapA1B2, &
      LineLineIntersectionTestParallelOverlapA2B1, LineLineIntersectionTestParallelOverlapA2B2, &
      LineLineIntersectionTestParallelA1B1, LineLineIntersectionTestParallelA1B2, LineLineIntersectionTestParallelA2B1, &
      LineLineIntersectionTestParallelA2B2, LineLineIntersectionTestIntersectA1, LineLineIntersectionTestIntersectA2, &
      LineLineIntersectionTestIntersectB1, LineLineIntersectionTestIntersectB2, LineLineIntersectionTestIntersectMiddle, &
      LineLineIntersectionTestIntersectProjectionCollide, LineLineIntersectionTestIntersectProjectionNoCollide
    implicit none

    call LineLineIntersectionTestOverlapA1B1()
    call LineLineIntersectionTestOverlapA1B2()
    call LineLineIntersectionTestOverlapA2B1()
    call LineLineIntersectionTestOverlapA2B2()
    call LineLineIntersectionTestSameLine()
    call LineLineIntersectionTestParallelOverlapA1B1()
    call LineLineIntersectionTestParallelOverlapA1B2()
    call LineLineIntersectionTestParallelOverlapA2B1()
    call LineLineIntersectionTestParallelOverlapA2B2()
    call LineLineIntersectionTestParallelA1B1()
    call LineLineIntersectionTestParallelA1B2()
    call LineLineIntersectionTestParallelA2B1()
    call LineLineIntersectionTestParallelA2B2()
    call LineLineIntersectionTestIntersectA1()
    call LineLineIntersectionTestIntersectA2()
    call LineLineIntersectionTestIntersectB1()
    call LineLineIntersectionTestIntersectB2()
    call LineLineIntersectionTestIntersectMiddle()
    call LineLineIntersectionTestIntersectProjectionCollide()
    call LineLineIntersectionTestIntersectProjectionNoCollide()

    print*, "SUCCESS: successful completion of all 20 line-line collision unit tests"

END PROGRAM
