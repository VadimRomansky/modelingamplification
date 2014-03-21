SUBROUTINE Progon(a,c,b,N,f,x)
implicit none
INTEGER j,N 
real*8 a(1:N),c(1:(N+1)),b(1:N),f(1:(N+1)),x(1:(N+1)),alfa(1:N),betta(1:N)
alfa(1:N)=0
betta(1:N)=0
alfa(1)=-b(1)/c(1)
betta(1)=f(1)/c(1)
  do j=2,N
     alfa(j)=(-b(j))/(a(j-1)*alfa(j-1)+c(j))      
	 betta(j)=(f(j)-a(j-1)*betta(j-1))/(a(j-1)*alfa(j-1)+c(j))
  end do
x(N+1)=(f(N+1)-a(N)*betta(N))/(a(N)*alfa(N)+c(N+1))
  do j=1,N
     x(N+1-j)=alfa(N-j+1)*x(N-j+2)+betta(N-j+1)
  end do
END SUBROUTINE
