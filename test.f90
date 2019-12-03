! npagane | risca lab | dec 2019 | test for sphere and line collisisons

! test line-line intersection module
PROGRAM Test
    use LineLineIntersection, only: LineLineIntersectionCalculation
    use SphereLineIntersection, only: SphereLineIntersectionCalculation
    use SphereSphereIntersection, only: SphereSphereIntersectionCalculation
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

END PROGRAM

