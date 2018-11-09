
 module lsq
!  LSQ - Interface module for LAPACK QR routines, and a special-purpose
!  subroutine for simple Y versus X linear regressions.
!
!  SYNOPSIS
!    use lsq
!
!  DESCRIPTION
!   linfit subroutine performs least squares regression for linear system or
!   power law and returns the required coefficients..
!
!  PUBLIC ROUTINES DEFINED
!    linfit - Y versus X linear or power law regression.
!
!  DEPENDENCI
!    precision - defines KINDs for single- and double-precision floating point.
!
!  REVISION HISTORY
!    10/15/18 - first implementation. 10/25/18 - final implementation
!
!  PROGRAMMER
!    Akshay Singhvi, asinghv3@asu.edu
!
   use precision
   implicit none

!  Generic interface (for single- or double-precision) for LAPACK
!  QR factorizations SGELS and DGELS.
   interface xgels
      subroutine sgels(tr,n,k,nrhs,x,ldx,y,ldy,work,lwork,info)
         import
         character, intent(in):: tr
         integer, intent(in):: n, k, ldx, ldy, nrhs
         integer, intent(inout):: lwork
         real(SP), intent(inout):: x(ldx,k), y(ldy,nrhs), work(lwork)
         integer, intent(out):: info
      end subroutine sgels

      subroutine dgels(tr, n, k, nrhs, x, ldx, y, ldy, work, lwork, info)
         import
         character, intent(in):: tr
         integer, intent(in):: n, k, ldx, ldy, nrhs
         integer, intent(inout):: lwork
         real(DP), intent(inout):: x(ldx,k), y(ldy,nrhs), work(lwork)
         integer, intent(out):: info
      end subroutine dgels
   end interface
!
   contains
!----------
   subroutine linfit(model, n, x, y, param, ierr)
!  LINFIT: uses the xgels interface to call dgels or sgels for solving the least
!  square proble from the intel mkl library 
!  Inputs:
!  model=character array which provides the option for power law or linear regression.
!  n: total number of data points
!  x: State matrix
!  y: measurments
!  Outputs:
!  param: Result coefficients b1 and b2 
!  ierr: error status, 0 if no error non zero represents a particular error 
   character(*), intent(in):: model
   integer, intent(in):: n
   real(DP):: x(n,2), y(n,1), work(100) 
   real(DP), intent(out):: param(2)
   integer, intent(out):: ierr
   integer:: k=2, nrhs=1,ldx, ldy, lwork=100, info
    character:: tr
    ldx=n
    ldy=n
   ! info=ierr     
 select case(model(1:1))
  case('P','p')
    x(:,2)= log10(x(:,2))   
    y(:,1)=log10(y(:,1))
    call xgels('N',n,2,1,x,ldx,y,ldy,work,lwork,info)
    param(1)=10**(y(1,1))
    param(2)=y(2,1)
    ierr=info   
  case('L','l')
    tr='N'
   call xgels(tr,n,k,nrhs,x,ldx,y,ldy,work,lwork,info)
    param(1)=y(1,1)
    param(2)=y(2,1)
    ierr=info 
 case default
    ierr=-1
  end select

   end subroutine linfit
!-----------------------
   end module lsq


