program DSAimpl
implicit none
integer, parameter :: Nx=200,Np=100, Nt=300000
real*8, parameter ::  a=10000.0D0, b=10000.0D0, Pmin=0.01D0, Pmax=100000.0D0
real*8  gn(Nx,Np-1), g(Nx,Np-1), x(1:Nx)
real*8 dt, h1, h2, ymin, dy, y
integer i,k,j,iter


ymin=dlog(Pmin)
dy=dlog(Pmax/Pmin)/Np

!создание сетки по х
h1=0.5D0*Nx/dlog(1.0D0+a)
h2=0.5D0*Nx/dlog(1.0D0+b)
do i=1,Nx/2
	x(i)=1.0D0-exp(-(1.0D0*i-0.5D0*Nx)/h1)
end do
do i=(Nx/2+1),Nx
	x(i)=exp((1.0D0*i-0.5D0*Nx)/h1)-1.0D0
end do
!создание сетки по х

!do i=1,Nx
!	x(i)=-a+(b+a)*i/Nx
!end do


dt=0.001D0
write(*,*) dt, x(Nx/2)-x(Nx/2-1), x(Nx/2+1)-x(Nx/2)
gn=0.0D0
open(8,file='data/gfull.dat')
	do i=1, Nx
		do k=1, Np-1
			read(8,*) gn(i,k)
		end do
	end do
close(UNIT=8)
!gn=0.0D0
g=0.0D0
iter=0
do j=1,Nt
	iter=iter+1
	call solver(a,ymin,x,dy,dt,Nx,Np,gn,g)
	gn(:,:)=g(:,:) 
	write(*,*) iter	
end do


	open(10,file='data/fp4.dat')
	do k=1, Np-1
		y=ymin+k*dy
		write(10,*) k, exp(y), g(Nx/2,k)*exp(y)
	end do
	close(UNIT=10)

	open(8,file='data/Nx.dat')
		write(8,*) Nx
	close(UNIT=8)

	open(12,file='data/Np.dat')
		write(12,*) Np
	close(UNIT=12)

	open(14,file='data/gfull.dat')
	do i=1, Nx
		do k=1, Np-1
			write(14,*) g(i,k)
		end do
	end do
	close(UNIT=14)

	open(18,file='data/x.dat')
	do i=1, Nx
		write(18,*) x(i)
	end do
	close(UNIT=18)


end