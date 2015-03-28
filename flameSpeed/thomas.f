c=======================================================================
      subroutine thomas (a,b,c,cwrk,r,n)
c-----------------------------------------------------------------------
c solve tridiagonal system with thomas algorithm.  tri(a,b,c)*x = r.
c solution returned in r array.
c input
c   a,b,c   arrays containing the tridiag matrix (these are not changed)
c   r       array containing the right hand side
c   n       size of the system
c output
c   r       solution returned in r array
c work
c   cwrk    work array so that original matrix is not destroyed
c-----------------------------------------------------------------------
      implicit none
      integer i,n
      real*8 a(n),b(n),c(n),cwrk(n),r(n),fac
c-----------------------------------------------------------------------
c     forward elimination

      fac     = 1.d0/b(1)
      cwrk(1) = c(1) * fac
      r(1)    = r(1) * fac

      do i = 2,n
        fac     = 1.d0/( b(i) - a(i)*cwrk(i-1) )

        cwrk(i) =                   c(i)   * fac
        r(i)    = ( r(i) - r(i-1) * a(i) ) * fac
      end do
c-----------------------------------------------------------------------
c     back substitution
c     put solution into r array, r(n) already has solution
      do i = n-1,1,-1
        r(i) = r(i) - cwrk(i)*r(i+1)
      end do
c-----------------------------------------------------------------------
      return
      end
c=======================================================================
