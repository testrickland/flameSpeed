SUBROUTINE generation(y,temp,omega)
    IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: y
    REAL, DIMENSION(3), INTENT(OUT) :: omega
    REAL, INTENT(IN) :: temp
    REAL, DIMENSION(3) :: mw, x, mol, concentration, massFrac
    REAL, DIMENSION(3,7) :: nuF, nuB
    REAL, DIMENSION(7) :: eF, eB, bF, bB, r, sF, sB, reactionRateF, reactionRateB
    REAL :: psi, mtot, ntot, p, concentrationTotal, prodConcen, reactConcen
    REAL, PARAMETER :: IGC=1.987204118, IGCJ=8.3144621
    INTEGER :: n, k, u, i, l, nEquations=7, nSpecies=3, firstReaction=1, lastReaction=7
!   files
    INTEGER :: sFdata=26, sBdata=27, bFdata=20, bBdata=21, nuFdata=22, nuBdata=23, eFdata=24, eBdata=25
    p = 101000
    concentrationTotal = (1e-6)*p/(IGCJ*temp)
    mw(1) = 16
    mw(2) = 32
    mw(3) = 48

    OPEN(sFdata, FILE='../generationCalculation/sFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(sFdata, 100) sF(i)
    ENDDO
    OPEN(sBdata, FILE='../generationCalculation/sBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(sBdata, 100) sB(i)
    ENDDO
    OPEN(bFdata, FILE='../generationCalculation/bFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(bFdata, 100) bF(i)
    ENDDO
    OPEN(bBdata, FILE='../generationCalculation/bBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(bBdata, 100) bB(i)
    ENDDO
    OPEN(nuFdata, FILE='../generationCalculation/nuFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(nuFdata, *) nuF(1,i), nuF(2,i), nuF(3,i)
    ENDDO
    OPEN(nuBdata, FILE='../generationCalculation/nuBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(nuBdata, *) nuB(1,i), nuB(2,i), nuB(3,i)
    ENDDO
    OPEN(eFdata, FILE='../generationCalculation/eFdata', STATUS='OLD')
    DO i=1,nEquations
    READ(eFdata, 100) eF(i)
    ENDDO
    OPEN(eBdata, FILE='../generationCalculation/eBdata', STATUS='OLD')
    DO i=1,nEquations
    READ(eBdata, 100) eB(i)
    ENDDO

    100 FORMAT (F10.2)
    CLOSE(bFdata)
    CLOSE(bBdata)
    CLOSE(nuFdata)
    CLOSE(nuBdata)
    CLOSE(eFdata)
    CLOSE(eBdata)
    CLOSE(sFdata)
    CLOSE(sBdata)

    omega(:) = 0
    DO n=firstReaction,lastReaction
    reactionRateF(n) = bF(n)*(temp**sF(n))*EXP(-eF(n)/(IGC*temp))
    reactionRateB(n) = bB(n)*(temp**sB(n))*EXP(-eB(n)/(IGC*temp))
    ENDDO
    DO n=firstReaction,lastReaction
    prodConcen = 1 
    reactConcen = 1 
    psi = 0
    mtot = 0
    ntot = 0
    DO l=1,nSpecies
    psi = psi + y(l)/mw(l)
    ENDDO
    DO l=1,nSpecies
    x(l) = y(l)/mw(l)/psi
    concentration(l) = x(l)*concentrationTotal
    prodConcen = prodConcen*(concentration(l)**nuB(l,n))
    reactConcen = reactConcen*(concentration(l)**nuF(l,n))
    ENDDO
    r(n) = reactionRateF(n)*(reactConcen) - reactionRateB(n)*(prodConcen)
    ENDDO
    DO l=1,nSpecies
    omega(l) = 0
    DO n=firstReaction,lastReaction
    omega(l) = omega(l) + mw(l)*(nuB(l,n) - nuF(l,n))*r(n)  
    ENDDO
!    WRITE(*,*) omega(l)
    ENDDO
END SUBROUTINE generation
