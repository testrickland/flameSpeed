PROGRAM flameSpeed
    IMPLICIT NONE
    REAL, DIMENSION(3) :: x
    REAL :: temp
    temp = 300
    x(1) = 0.3
    x(2) = 0.69
    x(3) = 0.01
CALL reaction(x(1),x(2),x(3),1200)

CONTAINS

SUBROUTINE reaction(x(1),x(2),x(3),temp)
    IMPLICIT NONE
    REAL, DIMENSION(3) :: omega, mw, y, massFrac, reactionRateF, reactionRateB, eF, eB, bF, bB
    REAL, INTENT(INOUT), DIMENSION(3) :: x
    REAL, INTENT(IN) :: T
    REAL, DIMENSION(3,7) :: nuF, nuB
    REAL, DIMENSION(7) :: rF, rB
    REAL :: rho, temp, p
    REAL, PARAMETER :: sF=2.5, sB=3.5, R=8.3144621
    INTEGER :: m, k, u, n, i, nEquations=7, nSpecies=3
!   files
    INTEGER :: bFdata=20, bBdata=21, nuFdata=22, nuBdata=23, eFdata=24, eBdata=25
    p = 100000
    mw(1) = 16
    mw(2) = 32
    mw(3) = 48 
    rho = (x(1)*MW(1) + x(2)*MW(2) + x(3)*MW(3))/1
    y(1) = (x(1)*MW(1))/rho
    y(2) = (x(2)*MW(2))/rho
    y(3) = (x(3)*MW(3))/rho
    omega(:) = 0
        
    OPEN(bFdata, FILE='bFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(bFdata, 100) bF(1)
    ENDDO
   
    OPEN(bBdata, FILE='bBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(bBdata, 100) bB(1)
    ENDDO
   
    OPEN(nuFdata, FILE='nuFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(nuFdata, *) nuF(1,i), nuF(2,i), nuF(3,i)
    ENDDO
    
    OPEN(nuBdata, FILE='nuBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(nuBdata, *) nuB(1,i), nuB(2,i), nuB(3,i)
    ENDDO
    
    OPEN(eFdata, FILE='eFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(eFdata, 100) eF(1)
    ENDDO
    
    OPEN(eBdata, FILE='eBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(eBdata, 100) eB(1)
    ENDDO

    100 FORMAT (F10.2)
    CLOSE(bFdata)
    CLOSE(bBdata)
    CLOSE(nuFdata)
    CLOSE(nuBdata)
    CLOSE(eFdata)
    CLOSE(eBdata)

    omega(:) = 0
    reactionRateF(1) = bF(1)*(temp**sF)*EXP(-eF(1)/(R*temp))
    reactionRateB(1) = bB(1)*(temp**sB)*EXP(-eB(1)/(R*temp))
    rF(1) = reactionRateF(1)*((y(1)*rho/mw(1))**nuF(1,1))*((y(2)*rho/mw(2))**nuF(2,1))*((y(3)*rho/mw(3))**nuF(3,1))
    rB(1) = reactionRateB(1)*(y(1)*rho/mw(1))**nuB(1,1)*(y(2)*rho/mw(2))**nuB(2,1)*(y(3)*rho/mw(3))**nuB(3,1)

    k = 1
    m = 1
    DO k=1,nSpecies
    omega(k) = omega(k) + mw(k)*(nuB(k,m) - nuF(k,m))*(rF(m) - rB(m))
    WRITE(*,*)  omega(k)
    ENDDO
END SUBROUTINE reaction
END PROGRAM flameSpeed


