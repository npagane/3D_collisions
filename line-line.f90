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
        real, parameter ::  tol = 1.0e-8 ! tolerance for cooccupancy (should be small to disallow overlap)
        real, parameter :: dist = 1.0e-8 ! tolerance for collision (pseudo thickness of line)
        real, dimension(3) :: pA ! closest point on A to B
        real, dimension(3) :: pB ! closest point on B to A
        real dotA1B1B2B1, dotB2B1A2A1, dotA1B1A2A1, dotB2B1B2B1, dotA2A1A2A1
        real muA, muB
        real, dimension(3) :: vecA, vecB, tA2, tB1, tB2
        integer LineLineIntersectionCalculation

        ! default value set to 0
        LineLineIntersectionCalculation = 0

        ! check for overlap of points
        if (ALL(ABS(A1-B1) <= tol) .OR. ALL(ABS(A1-B2) <= tol) .OR. ALL(ABS(A2-B1) <= tol) .OR. ALL(ABS(A2-B2) <= tol)) then
            !print*, "collision, point overlap"
            LineLineIntersectionCalculation = 1
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
                    LineLineIntersectionCalculation = 1
                    return ! quit early, coincident lines
                else
                    return ! quit early, non-coincident lines
                endif
            else
                !print*, "no collision, parallel but not overlap"
                return ! quit early, parallel and separated lines
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
        if ((sqrt(dot_product(pA-pB, pA-pB)) <= dist)  .AND. (abs(muA) <= 1) .AND. (abs(muB) <= 1) &
          .AND. .NOT. (abs(-1 - muA*muB) <= dist)) then ! check for -1 and 1 pairs  
            !print*, "collision, intersect"
            LineLineIntersectionCalculation = 1
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
        integer val 
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then 
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then 
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
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
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectMiddle"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectMiddle

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideZ
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/0,0,5/)
        real, dimension(3) :: B1 = (/-1,0,2/)
        real, dimension(3) :: B2 = (/1,0,2/)
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollideZ"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideZ

    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideZ
        implicit none
        real, dimension(3) :: A1 = (/-1,0,5/)
        real, dimension(3) :: A2 = (/1,0,5/)
        real, dimension(3) :: B1 = (/-1,0,2/)
        real, dimension(3) :: B2 = (/1,0,2/)
        integer val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideZ"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideZ

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideY
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/0,5,0/)
        real, dimension(3) :: B1 = (/-1,2,0/)
        real, dimension(3) :: B2 = (/1,2,0/)
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollideY"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideY

    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideY
        implicit none
        real, dimension(3) :: A1 = (/-1,5,0/)
        real, dimension(3) :: A2 = (/1,5,0/)
        real, dimension(3) :: B1 = (/-1,2,0/)
        real, dimension(3) :: B2 = (/1,2,0/)
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideY"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideY

    SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideX
        implicit none
        real, dimension(3) :: A1 = (/0,0,0/)
        real, dimension(3) :: A2 = (/5,0,0/)
        real, dimension(3) :: B1 = (/2,0,-1/)
        real, dimension(3) :: B2 = (/2,0,1/)
        integer val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 0) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionCollideX"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionCollideX
        
    SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideX
        implicit none
        real, dimension(3) :: A1 = (/5,0,-1/)
        real, dimension(3) :: A2 = (/5,0,1/)
        real, dimension(3) :: B1 = (/2,0,-1/)
        real, dimension(3) :: B2 = (/2,0,1/)
        integer val
        
        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjectionNoCollideX"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjectionNoCollideX

    SUBROUTINE LineLineIntersectionTestIntersectProjection
        implicit none
        real, dimension(3) :: A1 = (/-8.8910375546638831E-002 ,  0.45626156257315231      ,   6.6668324441179081/)!(/0,0,0/)
        real, dimension(3) :: A2 = (/5.5511151231257827E-017 ,  -4.4408920985006262E-016 ,   8.2500003278255480/)!(/5,5,0))
        real, dimension(3) :: B1 = (/8.5784016210477236E-018 ,  -6.3001451707972334E-017 ,   9.9000003933906555/)!(/0,2,2/)
        real, dimension(3) :: B2 = (/9.0729117231340004E-019 ,  -7.8930676690238801E-017 ,   11.550000458955765 /)!(/2,0,1/)
        integer val

        val = LineLineIntersectionCalculation(A1,A2,B1,B2)
        if (val == 1) then
            print*, "FAILURE: failed LineLineIntersectionTestIntersectProjection"
            stop
        endif

    END SUBROUTINE LineLineIntersectionTestIntersectProjection


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
      LineLineIntersectionTestIntersectProjectionCollideZ, LineLineIntersectionTestIntersectProjectionNoCollideZ, &
      LineLineIntersectionTestIntersectProjectionCollideY, LineLineIntersectionTestIntersectProjectionNoCollideY, &
      LineLineIntersectionTestIntersectProjectionCollideX, LineLineIntersectionTestIntersectProjectionNoCollideX, &
      LineLineIntersectionTestIntersectProjection
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
    call LineLineIntersectionTestIntersectProjectionCollideZ()
    call LineLineIntersectionTestIntersectProjectionNoCollideZ()
    call LineLineIntersectionTestIntersectProjectionCollideY()
    call LineLineIntersectionTestIntersectProjectionNoCollideY()
    call LineLineIntersectionTestIntersectProjectionCollideX()
    call LineLineIntersectionTestIntersectProjectionNoCollideX()
    call LineLineIntersectionTestIntersectProjection()

    print*, "SUCCESS: successful completion of all 25 line-line collision unit tests"

END PROGRAM
