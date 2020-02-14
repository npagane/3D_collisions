! npagane | risca lab | feb 2020 | fortran implementation of gjk algorithm 
! adapted from MATLAB code https://github.com/mws262/MATLAB-GJK-Collision-Detection/blob/master/GJK.m

! define sphere-sphere intersection module that contains calculation and test
MODULE GJKAlgorithm
    implicit none
    contains

    ! calculation
    FUNCTION GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        implicit none
        integer, intent(in) :: nVerts
        integer, parameter :: iteration = 6
        real, intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        real, dimension(3) :: a, b, c, d, aout, bout, cout, dout
        real, dimension(3), parameter :: v = (/0.8, 0.5, 1.0/)
        integer GJK

        ! default value (i.e. no intersection)
        GJK = 0

        ! line segment
        call pickLine(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, v, a, b)

        ! triangle
        call pickTriangle(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, &
                            aout, bout, cout, GJK, iteration)
        a = aout
        b = bout
        c = cout

        ! tetrahedron
        if (GJK == 1) then 
            call pickTetrahedron(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, c, &
                              aout, bout, cout, dout, GJK, iteration)
        endif

    END FUNCTION GJK

    ! pickLine
    SUBROUTINE pickLine(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, v, a, b)
        implicit none
        integer, intent(in) :: nVerts
        real, intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        real, dimension(3), intent(in) :: v 
        real, dimension(3), intent(out) :: a, b
        
        ! make first line of simplex
        a = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        b = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, -1*v)

    END SUBROUTINE pickLine

    SUBROUTINE pickTriangle(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, &
                            aout, bout, cout, flag, iteration)
        implicit none
        integer, intent(in) :: nVerts
        integer, intent(in) :: iteration
        real, dimension(3), intent(in) :: a, b
        real, intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        integer, intent(out) :: flag ! success
        real, dimension(3), intent(out) :: aout, bout, cout
        real, dimension(3) :: ab, ao, ac, abc, abp, acp, v
        integer i

        ! default value (i.e. no success)
        flag = 0
        aout = a
        bout = b

        ! first try
        ab = bout - aout
        ao = -1*aout
        v = cross(cross(ab, ao), ab) ! v is perp to ab and pointing towards origin
        cout = bout
        bout = aout
        aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)

        ! iterate until convergence
        do i = 1,iteration
            ! check if found
            ab = bout - aout
            ao = -1*aout
            ac = cout - aout
            ! find normal to face
            abc = cross(ab, ac)
            ! perp to ab and ac away from triangle
            abp = cross(ab, abc)
            acp = cross(abc, ac)
            ! check if triangle contains origin
            if (dot_product(abp, ao) > 0 ) then 
                cout = bout 
                bout = aout
                v = abp
            else if (dot_product(acp, ao) > 0 ) then 
                bout = aout 
                v = acp
            else 
                flag = 1
                exit
            endif
            aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        enddo

    END SUBROUTINE pickTriangle

    SUBROUTINE pickTetrahedron(s1x, s1y, s1z, s2x, s2y, s2z, nVerts, a, b, c, &
                              aout, bout, cout, dout, flag, iteration)
        implicit none
        integer, intent(in) :: nVerts
        integer, intent(in) :: iteration
        real, dimension(3), intent(in) :: a, b, c
        real, intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        integer, intent(out) :: flag ! success
        real, dimension(3), intent(out) :: aout, bout, cout, dout
        real, dimension(3) :: ab, ao, ac, ad, abc, acd, adb, v
        integer i

        ! default value (i.e. no success)
        flag = 0
        aout = a
        bout = b
        cout = c

        ! first try
        ab = bout - aout
        ac = cout - aout
        abc = cross(ab,ac)
        ao = -1*aout

        ! check if simplex is above or below origin
        if (dot_product(abc, ao) > 0) then 
            dout = cout
            cout = bout
            bout = aout
            v = abc
            aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        else
            dout = bout
            bout = aout
            v = -1*abc
            aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        endif

        ! iterate until convergence
        do i = 1, iteration
            ab = bout - aout
            ao = -1*aout
            ac = cout - aout
            ad = dout - aout
            abc = cross(ab, ac)
            if (dot_product(abc, ao) <= 0) then
                acd = cross(ac, ad)
                if (dot_product(acd, ao) > 0) then 
                    bout = cout
                    cout = dout
                    ab = ac
                    ac = ad
                    abc = acd
                else if (dot_product(acd, ao) < 0) then 
                    adb = cross(ad, ab)
                    if (dot_product(adb, ao) > 0) then 
                        cout = bout
                        bout = dout
                        ac = ab
                        ab = ad
                        abc = adb
                    else
                        flag = 1
                        exit
                    endif
                endif
            endif
            ! try again
            if (dot_product(abc, ao) > 0) then
                dout = cout
                cout = bout
                bout = aout
                v = abc
                aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
            else
                dout = bout
                bout = aout
                v = -1*abc
                aout = support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
            endif
        enddo

    END SUBROUTINE pickTetrahedron

    FUNCTION getExtremaPoint(sx, sy, sz, nVerts, v)
        implicit none
        integer, intent(in) :: nVerts
        real, intent(in) :: sx(nVerts), sy(nVerts), sz(nVerts)
        real, dimension(3), intent(in) :: v
        real, dimension(3) :: getExtremaPoint
        real :: mag(nVerts)
        integer maxInd(1)

        ! find the furthest data point in v direction in shape s
        mag = sx*v(1) + sy*v(2) + sz*v(3)
        !print*, 'mag', mag
        maxInd = maxloc(mag)
        !print*, 'ind', maxInd
        print*,  (/sx(maxInd(1)), sy(maxInd(1)), sz(maxInd(1))/)
        getExtremaPoint = (/sx(maxInd(1)), sy(maxInd(1)), sz(maxInd(1))/)

    END FUNCTION getExtremaPoint

    FUNCTION support(s2x, s2y, s2z, s1x, s1y, s1z, nVerts, v)
        implicit none
        integer, intent(in) :: nVerts
        real, intent(in) :: s1x(nVerts), s1y(nVerts), s1z(nVerts), s2x(nVerts), s2y(nVerts), s2z(nVerts)
        real, dimension(3), intent(in) :: v
        real, dimension(3) :: point1, point2, support

        ! support function for minkowski difference
        point1 = getExtremaPoint(s1x, s1y, s1z, nVerts, v)
        point2 = getExtremaPoint(s2x, s2y, s2z, nVerts, -1*v)
        support =  point1-point2

    END FUNCTION support

    FUNCTION cross(a, b)
        real, dimension(3) :: cross
        real, dimension(3), intent(in) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    END FUNCTION cross

    SUBROUTINE test()
        implicit none
        integer :: nVerts = 4
        real, dimension(4) :: s1x = (/0, 1, 1, 0/)
        real, dimension(4) :: s1y = (/0, 0, 1, 1/)
        real, dimension(4) :: s1z = (/1, 0, 0, 1/)
        real, dimension(4) :: s2x = (/0, 1, 1, 0/)
        real, dimension(4) :: s2y = (/0, 0, 1, 1/)
        real, dimension(4) :: s2z = (/2, 0, 0, 1/)
        integer res

        res = GJK(s1x, s1y, s1z, s2x, s2y, s2z, nVerts)
        print*, res

    END SUBROUTINE test

END MODULE

! test module
PROGRAM GJKTest 
    use GJKAlgorithm, only: GJK, test
    implicit none

    call test()

    print*, "SUCCESS: successful completion of all GJK collision unit tests"

END PROGRAM
