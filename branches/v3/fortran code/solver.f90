subroutine solver(a,ymin,x,dy,dt,Nx,Np,gn,g)
implicit none
integer Nx, Np, i, k
real*8 a, ymin, x(Nx), dy, dt
real*8 gn(Nx,Np-1), g(Nx,Np-1) 
real*8 lim, u, Qinj, kappa
real*8 a1(1:(Nx-1)),c1(1:Nx),b1(1:(Nx-1)),f(1:Nx), Xg(1:Nx)
real*8 y,dx, dxp, dxm, xp, xm, gip, gim, gkp, gkm


do k=1,Np-1
	y=ymin+k*dy	
    gkp=gn(1,k)
	if (k==1) then
		gkm=0.0D0
	else
		gkm=gn(1,k-1)
	end if
	dx=(x(2)+a)/2.0D0
	dxp=x(2)-x(1)
	dxm=x(1)+a
	xp=(x(2)+x(1))/2.0D0
	xm=(x(1)-a)/2.0D0
	gip=(gn(2,k)+gn(1,k))/2.0D0
	gim=gn(1,k)/2.0D0
	c1(1)=1.0D0+(dt/(2.0D0*dx))*(kappa(xp,y)/dxp+kappa(xm,y)/dxm)
	b1(1)=-(dt/(2.0D0*dx))*(kappa(xp,y)/dxp)
	f(1)=gn(1,k)+(dt/(2.0D0*dx))*(kappa(xp,y)*gn(2,k)/dxp)-(dt/(2.0D0*dx))*(kappa(xp,y)/dxp+kappa(xm,y)/dxm)*gn(1,k)-(dt/dx)*(u(xp)*gip-u(xm)*gim)+(dt/3.0D0)*((u(xp)-u(xm))/dx)*((gkp-gkm)/dy)+dt*Qinj(x(1),dx,y,dy)
	do i=2,Nx-1
		gkp=gn(i,k)
		if (k==1) then
			gkm=0.0D0
		else
			gkm=gn(i,k-1)
		end if	
		dx=(x(i+1)-x(i-1))/2.0D0
		dxp=x(i+1)-x(i)
		dxm=x(i)-x(i-1)
		xp=(x(i+1)+x(i))/2.0D0
		xm=(x(i)+x(i-1))/2.0D0
		gip=(gn(i+1,k)+gn(i,k))/2.0D0
		gim=(gn(i,k)+gn(i-1,k))/2.0D0
		a1(i-1)=-(dt/(2.0D0*dx))*(kappa(xm,y)/dxm)
		c1(i)=1.0D0+(dt/(2.0D0*dx))*(kappa(xp,y)/dxp+kappa(xm,y)/dxm)
		b1(i)=-(dt/(2.0D0*dx))*(kappa(xp,y)/dxp)
		f(i)=gn(i,k)+(dt/(2.0D0*dx))*(kappa(xp,y)*(gn(i+1,k)-gn(i,k))/dxp-kappa(xm,y)*(gn(i,k)-gn(i-1,k))/dxm)-(dt/dx)*(u(xp)*gip-u(xm)*gim)+(dt/3.0D0)*((u(xp)-u(xm))/dx)*((gkp-gkm)/dy)+dt*Qinj(x(i),dx,y,dy)						
	end do
    gkp=gn(Nx,k)
	if (k==1) then
		gkm=0.0D0
	else
		gkm=gn(Nx,k-1)
	end if
	dx=(x(Nx)-x(Nx-1))/2.0D0
	dxm=x(Nx)-x(Nx-1)
	xp=x(Nx)
	xm=(x(Nx)+x(Nx-1))/2.0D0
	gip=gn(Nx,k)
	gim=(gn(Nx,k)+gn(Nx-1,k))/2.0D0
	a1(Nx-1)=-(dt/(2.0D0*dx))*(kappa(xm,y)/dxp)
	c1(Nx)=1.0D0+(dt/(2.0D0*dx))*(kappa(xm,y)/dxm)
	f(Nx)=gn(Nx,k)-(dt/(2.0D0*dx))*(kappa(xm,y)*(gn(Nx,k)-gn(Nx-1,k))/dxm)-(dt/dx)*(u(xp)*gip-u(xm)*gim)+(dt/3.0D0)*((u(xp)-u(xm))/dx)*((gkp-gkm)/dy)+dt*Qinj(x(Nx),dx,y,dy)
	call Progon(a1,c1,b1,(Nx-1),f,Xg)
	g(:,k)=Xg(:)
end do


end subroutine