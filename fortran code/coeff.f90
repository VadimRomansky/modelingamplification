REAL*8 function u(x)
implicit none
real*8 x,dx
real*8 u1, u2

u1=1.0D0
u2=0.25D0*u1

if  (x<0.0D0) then
	u=u1
end if
if  (x>0.0D0) then
    u=u2
end if
end function



REAL*8 function kappa(x,y)
implicit none
real*8 x, y
real*8 kappa0

kappa0=0.1D0
kappa=kappa0*exp(y)
!kappa=kappa0

end function



REAL*8 function Qinj(x,dx,y,dy)
implicit none
real*8 x,dx, y, dy
real*8 u1, u2

Qinj=0.0D0

if (x<0.1D0*dx.AND.x>-0.1D0*dx.AND.y<0.5D0*dy.AND.y>-0.5D0*dy) then
	Qinj=1.0D0
end if

end function



REAL*8 function lim(fp,f0,fm)
implicit none
real*8 fp,f0,fm,r


if  ((fp-f0)==0.0D0) then
	lim=2.0
else
    r=(f0-fm)/(fp-f0)
    lim=(r+abs(r))/(1+abs(r))
end if

end function