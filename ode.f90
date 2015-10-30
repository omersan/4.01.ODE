!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Basic ODE solvers
!     Euler forward, Euler backward, Trapezoidal (Crank-Nicolson)
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 8, 2015
!-----------------------------------------------------------------------------!


program ode
implicit none
real*8 ::y0,h,tmax,t,f,lamda
integer::j,n,np
real*8,allocatable ::y1(:),y2(:),y3(:)

!Solve y'=lamda*y
lamda = -0.5d0

h = 1.00d0
tmax = 22.0d0
y0 = 1.0d0

n = nint(tmax/h)

allocate(y1(0:n))
allocate(y2(0:n))
allocate(y3(0:n))

!Initial condition
y1(0) = y0
y2(0) = y0
y3(0) = y0

!Solvers
do j=1,n
y1(j) = y1(j-1)*(1.0d0+h*lamda) !Euler forward
y2(j) = y2(j-1)/(1.0d0-h*lamda) !Euler backward
y3(j) = y3(j-1)*(1.0d0+h*lamda/2.0d0)/(1.0d0-h*lamda/2.0d0) !Crank-Nicolson
end do

!Plot
open(12, file="numerical.plt")
write(12,*)'variables ="t","EulerF","EulerB","CN"'
do j=0,n
	t = dfloat(j)*h
	write(12,*) t,y1(j),y2(j),y3(j)
end do
close(12)
   
! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="t","f"'
	do j=0,np
		t = dfloat(j)*(tmax)/(dfloat(np))
		f = y0*dexp(-0.5d0*t)
		write(12,*) t,f
	end do
close(12)



end






