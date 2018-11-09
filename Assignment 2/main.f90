!  DESCRIPTION
!    The basin dimension has been estimated using the following steps:
!   1) A 4096 by 4096 grid wit range -3 to 3 in both dimensions is
!   generated.
!   2) Referance henon map is created for the grid containing 1 and 0 
!   3) The grid is shifted to the righ and left by epsilon and henon map
!   is
!   generated respectively.
!   4) Reference basin is compared to both right and left henon maps and
!   number
!   of different elemnts is calculated 
!   5) The process is repeated for all the epsilons 
!   6) All the epsilons and corresponding number of differnt elemnts is
!   sent to
!   the linfit subroutine to compute the power law coefficients
!   7) The Basin dimension is calculated as 1+b2 and printed.
!
!   Subroutines:
!   1) basin_henon()- generates the henon map for a given N*N grid
!   2) algorithmD()- Computes the number of different elements between
!   the base
!   henon map and the epsilon shifted henon maps
!   3) combine()- combines all the function to implement the do loop
!   parallely 
!
!  REQUIRED DEPENDENCIES
!    precision - to define KINDs for single and double floating-point.
!    lsq- to call linfit for least squares regression
!
!  REVISION HISTORY
!   10/15/18 - First implementation.
!   10/25/18 - Final implementation
!
!  PROGRAMMER
!   Akshay Singhvi, asinghv3@asu.edu

!Generates a henon map with elements either 0 or 1 for a given N*N grid 
 pure subroutine basin_henon(basin,xmax,xmin,ymax,ymin,NGRID)
!   Input
!   xmax,xmin,ymax,ymin: grid boundaries in x and y direction
!   NGRID: Grid dimension
!   a,b: algoritm parameters
!   l: lockout parameter
!   k: iteration parameter 
!   outputs
!   basin: Henon basin map

    implicit none
    integer, intent(in):: NGRID
    integer::i,j,ii,k
    integer, intent(inout):: basin(NGRID,NGRID)
    real, intent(in):: xmax,xmin,ymax,ymin
    real:: dy, dx, a, b, l
    real, dimension(1:NGRID)::x0
    real, dimension(1:NGRID)::y0
    real, dimension(1:NGRID)::ytemp
    real:: X(NGRID,2), Y(NGRID,2), basin1(NGRID,NGRID)
    basin(:,:)=0
    dx=(xmax-xmin)/(NGRID-1)
    x0(1)=xmin
    k=50
    a=2.12
    b=-0.3
    l=100.0
!generates the x axis of the grid     
    do i=1,NGRID
        x0(i)=xmin+(i-1)*dx
    enddo

    dy=(ymax-ymin)/(NGRID-1)
!generates the y axis of the grid 
    do i=1,NGRID
        y0(i)=ymin+(i-1)*dx
    enddo

!loop that will form the B matrix by computing 1 column at a time 
    do i=1,NGRID
        ytemp(:)=y0(i)
        X(:,1)=x0
        X(:,2)=ytemp
        do j=1,k
            Y(:,2)=X(:,1)
            Y(:,1)=a-X(:,1)**2+b*X(:,2)
            X=Y
            basin1(:,i)=SQRT(Y(:,1)**2+Y(:,2)**2)
            do ii=1,NGRID
                if (basin1(ii,i)>100.0) then
                    basin(ii,i)=1
                end if
            enddo
        enddo
    enddo
end subroutine basin_henon

!Computes the number of different elements between the base henon map
!and the epsilon shifted henon maps
pure subroutine algorithmD(Bref,Bplus,Bminus,NGRID,N)
!   inputs
!   Bref: Refereence basin
!   Bplus: Right shifted basin 
!   Bminus: left shifted basin
!   NGRID: Grid dimension
!   Outputs
!   N: Total number of different elements    
    implicit none
    integer, intent(in):: NGRID
    integer, intent(out):: N    
    integer,intent(in):: Bref(NGRID,NGRID), Bplus(NGRID,NGRID), Bminus(NGRID,NGRID)
    integer:: i,j, N1,N2
    N1=0
    N2=0
    do i=1,NGRID
        do j=1,NGRID
            if (Bref(i,j)/= Bminus(i,j)) then
                N1=N1+1
            end if
            if (Bref(i,j)/=Bplus(i,j)) then
                N2=N2+1
            end if
        enddo
    enddo
    N=N1+N2
end subroutine algorithmD

!combines all the function to implement the do loop parallely
subroutine combine(Bref,xmax,xmin,ymax,ymin,NGRID,epsi,N)
     
     integer, intent(in) :: NGRID,Bref(NGRID,NGRID)
     integer::Bplus(NGRID,NGRID),Bminus(NGRID,NGRID)
     real, intent(in)::xmax,xmin,ymax,ymin,epsi
     integer, intent(out):: N
!Gnenrates the Bplus matrix which corresponds to the right shifted grid     
     call basin_henon(Bplus,xmax+epsi,xmin+epsi,ymax,ymin,NGRID)
!generates the Bminus matrix which correspondes to the left shifted grid 
     call basin_henon(Bminus,xmax-epsi,xmin-epsi,ymax,ymin,NGRID)
!gives the number of different elemnts between the B matrix and Bplus,Bminus
     call algorithmD(Bref,Bplus,Bminus,NGRID,N)
    
end subroutine combine

!This program implements the entire algorithm and prints the estimated
!basin dimension
program main
    use lsq
    implicit none
    integer, parameter:: NGRID=4096
    real:: xmax=3.0, xmin=-3.0, ymax=3.0,ymin=-3.0,epsi(10)
    real(DP)::x(10,2), y(10,1), param(2), b_dim
    integer::Bref(NGRID,NGRID),i,N_final(10),ierr,tnr
    epsi=(/2.0**(-12),2.0**(-13),2.0**(-14),2.0**(-15),2.0**(-16),2.0**(-17),2.0**(-18),2.0**(-19),2.0**(-20),2.0**(-21)/)
!genrates the B matrix     
    call basin_henon(Bref,xmax,xmin,ymax,ymin,NGRID)
! this loop runs for all the 10 epsilon values and the fucntion produces
! N_final
! as number of different elemnts between Bref and Bplus,Bminus
    call omp_set_num_threads(4)
    !$OMP PARALLEL DO 
    do i=1,10    
        call combine(Bref,xmax,xmin,ymax,ymin,NGRID,epsi(i),N_final(i))
    enddo
    x(:,1)=1
    x(:,2)=epsi
    y(:,1)=N_final    
!The function performos least squares regression for power law and
!returns the
!coefficients C and P (Power Law: Y=CX**P)
    call linfit('Power',10,x,y,param,ierr)
    b_dim=1.0+param(2)
    write(*,*) "The estimated dimension of the basin is: ",  b_dim
end program main
