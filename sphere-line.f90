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
        b = 2*dot_product(A2-A1, A1-B1)
        c = dot_product(B1, B1) + dot_product(A1, A1) - 2*dot_product(B1, A1) - r*r
        discr = b*b - 4*a*c
        if ( discr < 0 ) then 
            return
        endif
        up = (-b + sqrt(discr))/(2*a)
        um = (-b - sqrt(discr))/(2*a)
        if ((um > 1 .AND. up < 0) .OR. (um < 0 .AND. up > 1)) then 
            SphereLineIntersectionCalculation = .TRUE.
            !print*, "line in sphere"
            return
        else if ((um >= 0 .AND. um <= 1) .OR. (up >= 0 .AND. up <= 1)) then 
            SphereLineIntersectionCalculation = .TRUE.
            !print*, "one or more intersections"
            return
        endif
        return 

    END FUNCTION SphereLineIntersectionCalculation

    ! testing routines !
    SUBROUTINE SphereLineIntersectionTestLineInside
        implicit none
        real, dimension(3) :: A1 = (/-1,0,0/)
        real, dimension(3) :: A2 = (/-2,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real :: r = 20.0
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineInside"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineInside

    SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA1
        implicit none
        real, dimension(3) :: A1 = (/-2,0,0/)
        real, dimension(3) :: A2 = (/-1,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real :: r = 2
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineInsideEdgeA1"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA1

    SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA2
        implicit none
        real, dimension(3) :: A1 = (/-1,0,0/)
        real, dimension(3) :: A2 = (/-2,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real :: r = 2
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineInsideEdgeA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineInsideEdgeA2

    SUBROUTINE SphereLineIntersectionTestLineTangent
        implicit none
        real, dimension(3) :: A1 = (/-1,0,0/)
        real, dimension(3) :: A2 = (/1,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineTangent"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineTangent

    SUBROUTINE SphereLineIntersectionTestLineOutsideA1
        implicit none
        real, dimension(3) :: A1 = (/1,0,0/)
        real, dimension(3) :: A2 = (/2,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA1"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineOutsideA1

    SUBROUTINE SphereLineIntersectionTestLineOutsideA2
        implicit none
        real, dimension(3) :: A1 = (/4,0,0/)
        real, dimension(3) :: A2 = (/3,0,0/)
        real, dimension(3) :: B1 = (/1,0,0/)
        real :: r = 2
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineOutsideA2

    SUBROUTINE SphereLineIntersectionTestLineOutsideBoth
        implicit none
        real, dimension(3) :: A1 = (/5,0,0/)
        real, dimension(3) :: A2 = (/-5,0,0/)
        real, dimension(3) :: B1 = (/0,0,0/)
        real :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .TRUE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideBoth"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineOutsideBoth


    SUBROUTINE SphereLineIntersectionTestLineCloseA1
        implicit none
        real, dimension(3) :: A1 = (/12,0,0/)
        real, dimension(3) :: A2 = (/11.001,0.0,0.0/)
        real, dimension(3) :: B1 = (/10,10,10/)
        real :: r = 1
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineCloseA1"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineCloseA1

    SUBROUTINE SphereLineIntersectionTestLineCloseA2
        implicit none
        real, dimension(3) :: A1 = (/-55,0,0/)
        real, dimension(3) :: A2 = (/-50.003,0.0,0.0/)
        real, dimension(3) :: B1 = (/-10,-10,-10/)
        real :: r = 40
        logical val

        val = SphereLineIntersectionCalculation(A1, A2, B1, r)
        if (val .NEQV. .FALSE.) then
            print*, "FAILURE: failed LineLineIntersectionTestLineOutsideA2"
            stop
        endif

    END SUBROUTINE SphereLineIntersectionTestLineCloseA2

END MODULE SphereLineIntersection

! test sphere-line intersection module
PROGRAM SphereLineTest 
    use SphereLineIntersection, only: SphereLineIntersectionTestLineInside, SphereLineIntersectionTestLineInsideEdgeA1, &
      SphereLineIntersectionTestLineInsideEdgeA2, SphereLineIntersectionTestLineTangent, &
      SphereLineIntersectionTestLineOutsideA1, SphereLineIntersectionTestLineOutsideA2, &
      SphereLineIntersectionTestLineOutsideBoth, SphereLineIntersectionTestLineCloseA1, &
      SphereLineIntersectionTestLineCloseA2
    implicit none

    call SphereLineIntersectionTestLineInside()
    call SphereLineIntersectionTestLineInsideEdgeA1()
    call SphereLineIntersectionTestLineInsideEdgeA2()
    call SphereLineIntersectionTestLineTangent()
    call SphereLineIntersectionTestLineOutsideA1()
    call SphereLineIntersectionTestLineOutsideA2()
    call SphereLineIntersectionTestLineOutsideBoth()
    call SphereLineIntersectionTestLineCloseA1()
    call SphereLineIntersectionTestLineCloseA2()

    print*, "SUCCESS: successful completion of all 9 sphere-line collision unit tests"

END PROGRAM


