PROGRAM flameSpeed 
    IMPLICIT NONE
    REAL :: deltaT, deltax, loc, locnot=2.85e-3, eta, D = 0.3, cp = 0.2524, p=101000
    REAL, DIMENSION(3) :: mw, h
    REAL, PARAMETER :: pi=3.142, domainLength=5.7e-2, endtime=1e-2, IGC=1.987204118, IGCJ=8.3144621
    INTEGER :: l, j, time 
    INTEGER, PARAMETER :: mwdata=43, hdata=44, nSpecies=3, Nstepst=100, Nstepsx=100
    REAL, DIMENSION(3,Nstepsx) :: y, omega
    REAL, DIMENSION(Nstepsx) :: temp, energyEvolution, rho, at, bt, ct, cwrk, r
    REAl, DIMENSION(Nstepsx,Nstepsx) :: A
    deltaT = endtime/Nstepst
    deltax = domainLength/Nstepsx
    deltaT = 0.00001

    OPEN(mwdata, FILE='../generationCalculation/mwdata', STATUS='OLD')
    DO j=1,nSpecies
    READ(mwdata, 100) mw(j)
    ENDDO
    OPEN(hdata, FILE='../generationCalculation/hdata', STATUS='OLD')
    DO j=1,nSpecies
    READ(hdata, 100) h(j)
    ENDDO
    100 FORMAT (F10.2)
    CLOSE(mwdata)
    CLOSE(hdata)

!initialize
    eta = 0
    j = 0
    DO 
    j = j+1
    IF ( eta < 1.20 ) THEN
    temp(j) = 300.0*(1.0+(19.0*COS(pi/2.0*(eta/1.2)**7)**5)/6.0)
    y(1,j) = 5e-3*COS(pi/2*(eta/1.2)**7)**5
    y(2,j) = (2+COS(pi/2*(eta/1.2)**7)**5)/3 - y(1,j)
    y(3,j) = 1 - y(2,j) - y(3,j)
    ELSE IF ( eta < ((Nstepsx*deltax)/locnot) ) THEN
    temp(j) = 300.0
    y(1,j) = 0.0 
    y(2,j) = 2.0/3.0
    y(3,j) = 1.0/3.0
    ELSE
    EXIT
    ENDIF
    eta = (eta*locnot + deltax)/locnot
    WRITE(*,*) temp(j)
    ENDDO
    DO j=1,Nstepsx
    WRITE(*,*) temp(j), y(2,j), y(3,j)
    ENDDO

!start the loop
    DO time=1,2
    DO j=1,Nstepsx

!omega calculation
    CALL generation(y(:,j),temp(j),omega(:,j))

!rho calculation
    rho(j) = 0
    DO l=1,nSpecies
    rho(j) = rho(j) + y(l,j)/mw(j)
    ENDDO
    rho(j) = (1e-6)*p/(IGCJ*temp(j))/rho(j)
    ENDDO

!mass fraction calculation
    DO l=1,nSpecies
!setting up the mass fraction matrix
    j=1
    r(j) = (1-(2*D*deltat/deltax**2)/2)*y(l,j) + omega(l,j)*deltat/rho(j)
    A(j,j) = 1+(2*D*deltat/deltax**2)/2
    A(j,j+1) = -(2*D*deltat/deltax**2)/2
    at(j) = 0
    ct(j) = A(j,j+1)
    bt(j) = A(j,j)
    DO j=2,(Nstepsx-1)
    r(j) = (1-(2*D*deltat/deltax**2)/2)*y(l,j) + omega(l,j)*deltat/rho(j)
    A(j,j) = 1+2*D*deltat/deltax**2/2
    A(j,j+1) = -D*deltat/deltax**2/2
    A(j,j-1) = -D*deltat/deltax**2/2
    at(j) = 0
    ct(j) = A(1,2)
    bt(j) = A(1,1)
    ENDDO
    j=Nstepsx
    r(j) = (1-2*D*deltat/deltax**2/2)*y(l,j) + omega(l,j)*deltat/rho(j)
    A(j,j) = 1+2*D*deltat/deltax**2/2
    A(j,j-1) = -2*D*deltat/deltax**2/2
    ct(j) = 0
    at(j) = A(j,j-1)
    bt(j) = A(j,j)
!solving mass fraction    
    CALL thomas(at,bt,ct,cwrk,r,Nstepsx) 
!writing solution
    DO j=1,Nstepsx
    y(l,j) = r(j)
    ENDDO
    ENDDO		!end of mass fraction calculation

!temperature/energy calculation
!solving for the energy of evolution
    DO j=1,Nstepsx
    energyEvolution(j) = 0
    DO l=1,nSpecies
    energyEvolution(j) = energyEvolution(j) + omega(l,j)*h(l)/(rho(j)*cp)
    ENDDO
    ENDDO
!setting up temperature/energy matrix 
    j=1
    r(j) = (1-(2*D*deltat/deltax**2)/2)*temp(j) + energyEvolution(j)*deltat
    A(j,j) = 1+(2*D*deltat/deltax**2)/2
    A(j,j+1) = -(2*D*deltat/deltax**2)/2
    at(j) = 0
    ct(j) = A(j,j+1)
    bt(j) = A(j,j)
    DO j=2,(Nstepsx-1)
    r(j) = (1-(2*D*deltat/deltax**2)/2)*temp(j) + energyEvolution(j)*deltat
    A(j,j) = 1+2*D*deltat/deltax**2/2
    A(j,j+1) = -D*deltat/deltax**2/2
    A(j,j-1) = -D*deltat/deltax**2/2
    at(j) = 0
    ct(j) = A(1,2)
    bt(j) = A(1,1)
    ENDDO
    j=Nstepsx
    r(j) = (1-2*D*deltat/deltax**2/2)*temp(j) + energyEvolution(j)*deltat
    A(j,j) = 1+2*D*deltat/deltax**2/2
    A(j,j-1) = -2*D*deltat/deltax**2/2
    ct(j) = 0
    at(j) = A(j,j-1)
    bt(j) = A(j,j)
!solving temperature    
    CALL thomas(at,bt,ct,cwrk,r,Nstepsx) 
!writing solution
    DO j=1,Nstepsx
    temp(j) = r(j)
    ENDDO
    ENDDO		!end time loop
    

CONTAINS

ENDPROGRAM flameSpeed 
