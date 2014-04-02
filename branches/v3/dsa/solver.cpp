#include "coeff.h"
#include "progon.h"

void solver(double a, double ymin, double* x, double dy, double dt, int Nx, int Np, double** gn, double** g){
	double* a1 = new double[Nx-1];
	double* c1 = new double[Nx];
	double* b1 = new double[Nx-1];
	double* f = new double[Nx];
	double* Xg = new double[Nx];
	double y,dx, dxp, dxm, xp, xm, gip, gim, gkp, gkm;


	for(int k =0; k < Np-1; ++k){
		y = ymin + k*dy;
		gkp = gn[0][k];
		if (k == 0) {
			gkm=0.0;
		} else {
			gkm = gn[1][k-1];
		}
		dx = (x[1] + a)/2;
		dxp = x[1] - x[0];
		dxm = x[0] + a;
		xp = (x[1] + x[0])/2;
		xm = (x[0] - a)/2;
		gip = (gn[1][k] + gn[0][k])/2;
		gim = gn[0][k]/2;
		c1[0] = 1 + (dt/(2*dx))*(kappa(xp,y)/dxp + kappa(xm,y)/dxm);
		b1[0] = -(dt/(2*dx))*(kappa(xp,y)/dxp);
		f[0]=gn[0][k] + (dt/(2*dx))*(kappa(xp,y)*gn[1][k]/dxp)
					  - (dt/(2*dx))*(kappa(xp,y)/dxp+kappa(xm,y)/dxm)*gn[0][k] 
					  - (dt/dx)*(u(xp)*gip - u(xm)*gim)
					  + (dt/3)*((u(xp)-u(xm))/dx)*((gkp-gkm)/dy)+dt*Qinj(x[0],dx,y,dy);
		for(int i = 1; i < Nx-1; ++i){
			gkp = gn[i][k];
			if (k == 0) {
				gkm=0;
			} else {
				gkm=gn[i][k-1];
			}
			dx = (x[i+1] - x[i-1])/2;
			dxp=x[i+1]-x[i];
			dxm=x[i]-x[i-1];
			xp=(x[i+1]+x[i])/2;
			xm=(x[i]+x[i-1])/2;
			gip=(gn[i+1][k] + gn[i][k])/2;
			gim=(gn[i][k] + gn[i-1][k])/2;
			a1[i-1] = -(dt/(2*dx))*(kappa(xm,y)/dxm);
			c1[i] = 1 + (dt/(2*dx))*(kappa(xp,y)/dxp+kappa(xm,y)/dxm);
			b1[i] = -(dt/(2*dx))*(kappa(xp,y)/dxp);
			f[i] = gn[i][k] + (dt/(2*dx))*(kappa(xp,y)*(gn[i+1][k] - gn[i][k])/dxp
							- kappa(xm,y)*(gn[i][k] - gn[i-1][k])/dxm)
							- (dt/dx)*(u(xp)*gip - u(xm)*gim)
							+ (dt/3)*((u(xp) - u(xm))/dx)*((gkp - gkm)/dy)
							+ dt*Qinj(x[i],dx,y,dy);
		}
		gkp = gn[Nx-1][k];
		if (k==0) {
			gkm=0;
		} else {
			gkm = gn[Nx-1][k-1];
		}
		dx = (x[Nx-1] - x[Nx-2])/2;
		dxm = x[Nx-1] - x[Nx-2];
		xp = x[Nx-1];
		xm = (x[Nx-1] + x[Nx-2])/2;
		gip = gn[Nx-1][k];
		gim = (gn[Nx-1][k] + gn[Nx-2][k])/2;
		a1[Nx-2] = -(dt/(2*dx))*(kappa(xm,y)/dxp);
		c1[Nx-1] = 1 + (dt/(2*dx))*(kappa(xm,y)/dxm);
		f[Nx-1] = gn[Nx-1][k] - (dt/(2*dx))*(kappa(xm,y)*(gn[Nx-1][k] - gn[Nx-2][k])/dxm)
							  - (dt/dx)*(u(xp)*gip - u(xm)*gim)
							  + (dt/3)*((u(xp) - u(xm))/dx)*((gkp - gkm)/dy)
							  + dt*Qinj(x[Nx-1],dx,y,dy);
		progon(a1,c1,b1,(Nx-1),f,Xg);
		for(int i = 0; i < Nx; ++i){
			g[i][k]= Xg[i];
		}
	}
	delete[] a1;
	delete[] c1;
	delete[] b1;
	delete[] f;
	delete[] Xg;

}