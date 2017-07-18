	subroutine bsstep(y,dydx,nv,x,htry,eps,i0,yscal,hdid,hnext)
	implicit real*8 (a-h,o-z)
	parameter (nmax=10,imax=14,nuse=7,one=1.d0,shrink=0.95d0,grow=1.2d0)
	dimension y(nv),dydx(nv),yscal(nv)
	dimension nseq(imax)
c	dimension yerr(nmax),ysav(nmax),dysav(nmax),yseq(nmax)
c	The above line is that specified in Numerical Recipes.  However, it
c	is inconsistent with the dimension statements in subroutines
c	MMID and RZEXTR.  The next line has dimension statements
c	consistent with those subroutines.
	dimension yerr(5),ysav(5),dysav(5),yseq(5)
	data nseq/2,4,6,8,12,16,24,32,48,64,96,128,192,256/
	h=htry
	xsav=x
	do 11 i=1,nv
	  ysav(i)=y(i)
	  dysav(i)=dydx(i)
11	continue
1	do 10 i=1,imax
c		write(6,*)'BSSTEP: entering mmid'
	  call mmid(ysav,dysav,nv,i0,xsav,h,nseq(i),yseq)
c		write(6,*)'bsstep: out of mmid'
c		write(6,*)'yseq=',yseq
	  xest=(h/(dble(nseq(i))))**2
	  call rzextr(i,xest,yseq,y,yerr,nv,nuse)
c		write(6,*)'out of rzextr: yerr=',yerr
c		write(6,*)'yscal=',yscal
	  errmax=0.d0
	do 12 j=1,nv
	  t=dabs(yerr(j)/yscal(j))
	  if(errmax.lt.t) errmax=t
c	  errmax=max(errmax,dabs(yerr(j)/yscal(j)))
c		write(6,*)'j=',j,yerr(j)/yscal(j)
12	continue
c		write(6,*)'** i=',i,' ** errmax=',errmax
c		write(6,*)'eps=',eps
c		write(6,*)'  '
	  errmax=errmax/eps
c		write(6,*)'NEW errmax=',errmax,'one=',one
	  if(errmax.lt.one) then
c		write(6,*)'yerr=',yerr,'yscal=',yscal
	  	x=x+h
		hdid=h
		if(i.eq.nuse) then
		  hnext=h*shrink
		else if(i.eq.nuse-1) then
		  hnext=h*grow
		else
		  hnext=(h*dble(nseq(nuse-1)))/dble(nseq(i))
		endif
		return
	  endif
10	continue
	h=0.25d0*h/dble(2**((imax-nuse)/2))
	if(x+h.eq.x) pause 'step size underflow'
	go to 1
	end	
	subroutine mmid(y,dydx,nvar,i0,xs,htot,nstep,yout)
	implicit real*8 (a-h,o-z)
	parameter (nmax=10)
	dimension y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax)
	h=htot/dble(nstep)
	do 11 i=1,nvar
	  ym(i)=y(i)
	  yn(i)=y(i)+h*dydx(i)
11	continue
	x=xs+h
	if(nvar.eq.4) call derivs4(i0,x,yn,yout)
	if(nvar.eq.5) call derivs5(i0,x,yn,yout)
	if(nvar.eq.6) call derivs6(i0,x,yn,yout)
	h2=2.d0*h
	do 13 n=2,nstep
	do 12 i=1,nvar
	  swap=ym(i)+h2*yout(i)
	  ym(i)=yn(i)
	  yn(i)=swap
12	continue
	  x=x+h
	if(nvar.eq.4) call derivs4(i0,x,yn,yout)
	if(nvar.eq.5) call derivs5(i0,x,yn,yout)
	if(nvar.eq.6) call derivs6(i0,x,yn,yout)
13	continue
	do 14 i=1,nvar
	yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14	continue
	return
	end
	subroutine derivs4(j,r,yn,yout)
	real*8 r,s0
	real*8 yn(4),yout(4),aj(4,4)
	real*4 kappa,mu,mu2,mupr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
	call amat4(aj,j,r)
	call prod(aj,yn,4,yout) 
	return
	end
	subroutine derivs5(j,r,yn,yout)
	real*8 r,s0
	real*8 yn(5),yout(5),aj(5,5)
	real*4 kappa,mu,mu2,mupr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
	call amat5(aj,j,r)
	call prod(aj,yn,5,yout) 
	return
	end
	subroutine derivs6(j,r,yn,yout)
	real*8 r,s0
	real*8 yn(6),yout(6),aj(6,6)
	real*4 kappa,mu,mu2,mupr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
	call amat6(aj,j,r)
	call prod(aj,yn,6,yout) 
	return
	end
    	subroutine amat4(a,i0,r)
c 	Use to determine A(r) (returned in a) in layer at radius r.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu,mu2,mupr 
	real*8 mus,lams,sfac,flp,s0
	real*8 rho,g0
	real*8 tau1,tau2
	real*8 a(4,4) 
	real*8 fl21,r,biga,r2
	common/erad/bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
c
c	write(6,*)'AMAT4: kappa=',kappa
c	write(6,*)'AMAT4: mu=',mu
c	write(6,*)'AMAT4: mu2=',mu2
c	write(6,*)'AMAT4: eta1=',eta1 
c
	tau1=dble(mu(i0)/eta1(i0))
	tau2=dble(mu2(i0)/eta(i0))
	mus=dble(mu(i0))*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+
     &	dble(mu(i0))*s0/dble(eta(i0)))
	lams=dble(kappa(i0)-2.*mus/3.)
	biga=lams+2.d0*mus   
	fl21=dble(l*(l+1))
	flp=dble(l)
	r2=r*r
c		write(6,*)'AMAT4: s0,tau1,tau2=',s0,tau1,tau2,'eta1(',i0,')=',eta1(i0)
c		write(6,*)'AMAT4: mus,lams,r,l=',mus,lams,r,l
c	include gravitation for layers 80-94.  Otherwise
c	use layer matrices for the nongraviational case.
c**
c	Use density contrast of 0.5 g-cm**{-3} between asthenosphere
c	and lithospshere.  This version: assume that earth radius
c	is 2250 km.
cc	if(i0.ge.79) rho=3.3d0
cc	if(i0.lt.79) rho=3.8d0
c	rho=0.5d0
c	g0=980.d0*(2250.d+5)*1.e-11
	rho=dble(dens(i0))
c	g0=980.d0*(6371.d+5)*1.e-11
	g0=980.d0*dble(bigr*1.d+5)*1.e-11
c	  g0=0.d0
	r2=r*r
	a(1,1)=-2.d0*lams/(r*biga)
	a(1,2)=1.d0/biga
	a(1,3)=fl21*lams/(r*biga)
	a(1,4)=0.d0
	a(2,1)=(4.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,1)-4.d0*rho*g0/r2
	a(2,2)=-2.d0/r+(2.d0/r)*lams*a(1,2)
	a(2,3)=-fl21*(2.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,3)
     &	+rho*g0*fl21/r2
	a(2,4)=fl21/r
	a(3,1)=-1.d0/r
	a(3,2)=0.d0
	a(3,3)=1.d0/r
	a(3,4)=1.d0/mus
	a(4,1)=-(lams/r)*a(1,1)+2.d0*(mus-biga)/r2+rho*g0/r2 
	a(4,2)=-(lams/r)*a(1,2)
	a(4,3)=-(lams/r)*a(1,3)+(fl21*biga-2.d0*mus)/r2
	a(4,4)=-3.d0/r
c	  write(6,*)'A: aj=',aj
	return
	end 
    	subroutine amat5(aj,i0,r)
c 	Use to determine A(r) (returned in aj) in layer at radius r.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu,mu2,mupr 
	real*8 mus,lams,sfac,flp,s0
	real*8 rho,g0
	real*8 tau1,tau2
	real*8 aj(5,5),a(5,5) 
	real*8 fac
	real*8 fl21,r,biga,r2
	dimension id(4,4)
	common/erad/bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
	tau1=dble(mu(i0)/eta1(i0))
	tau2=dble(mu2(i0)/eta(i0))
	mus=dble(mu(i0))*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+
     &	dble(mu(i0))*s0/dble(eta(i0)))
	lams=dble(kappa(i0)-2.*mus/3.)
	biga=lams+2.d0*mus   
	fl21=dble(l*(l+1))
	flp=dble(l)
	r2=r*r
c	include gravitation for layers 80-94.  Otherwise
c	use layer matrices for the nongraviational case.
c**
c	Use density contrast of 0.5 g-cm**{-3} between asthenosphere
c	and lithospshere.  This version: assume that earth radius
c	is 2250 km.
cc	if(i0.ge.79) rho=3.3d0
cc	if(i0.lt.79) rho=3.8d0
c	rho=0.5d0
cc	g0=980.d0*(2250.d+5)*1.e-11
	rho=dble(dens(i0))
	g0=980.d0*dble(bigr*1.d+5)*1.e-11
c	  g0=0.d0
	r2=r*r
	do 4 i=1,5
	do 3 j=1,5
3	aj(i,j)=0.d0
4	continue
c	specify id-array, which identifies the index of the second order minor.
	id(1,1)=0
	id(1,2)=1
	id(1,3)=2
	id(1,4)=3
	id(2,1)=id(1,2)
	id(2,2)=0
	id(2,3)=4
	id(2,4)=5
	id(3,1)=id(1,3)
	id(3,2)=id(2,3)
	id(3,3)=0
	id(3,4)=6
	id(4,1)=id(1,4)
	id(4,2)=id(2,4)
	id(4,3)=id(3,4)
	id(4,4)=0
c
	a(1,1)=-2.d0*lams/(r*biga)
	a(1,2)=1.d0/biga
	a(1,3)=fl21*lams/(r*biga)
	a(1,4)=0.d0
	a(2,1)=(4.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,1)-4.d0*rho*g0/r2
	a(2,2)=-2.d0/r+(2.d0/r)*lams*a(1,2)
	a(2,3)=-fl21*(2.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,3)
     &	+rho*g0*fl21/r2
	a(2,4)=fl21/r
	a(3,1)=-1.d0/r
	a(3,2)=0.d0
	a(3,3)=1.d0/r
	a(3,4)=1.d0/mus
	a(4,1)=-(lams/r)*a(1,1)+2.d0*(mus-biga)/r2+rho*g0/r2 
	a(4,2)=-(lams/r)*a(1,2)
	a(4,3)=-(lams/r)*a(1,3)+(fl21*biga-2.d0*mus)/r2
	a(4,4)=-3.d0/r
	do 5 j=1,4
	do 10 k=j,4
	do 15 n=1,4
	m0=id(j,k)
	if(m0.eq.0.or.m0.eq.6) go to 15
	m=id(n,k)
	if(m.eq.0) go to 16
	fac=1.d0
	if(m.eq.6) fac=-1.d0/fl21
	if(m.eq.6) m=1
	if(n.gt.k) fac=-fac
	aj(m0,m)=aj(m0,m)+a(j,n)*fac
c	  write(6,*)'m0,m,=',m0,m,'j,n=',j,n,'fac=',fac
16	m=id(j,n)
	if(m.eq.0) go to 15
	fac=1.d0
	if(m.eq.6) fac=-1.d0/fl21
	if(m.eq.6) m=1
	if(j.gt.n) fac=-fac
	aj(m0,m)=aj(m0,m)+a(k,n)*fac
c	  write(6,*)'m0,m,=',m0,m,'k,n=',k,n,'fac=',fac
15	continue
10	continue
5	continue
c	  write(6,*)'A: aj=',aj
	return
	end 
    	subroutine amat6(a,i0,r)
c 	Use to determine A(r) (returned in a) in layer at radius r.
c	dY(r)/dr = A(r)Y.
	real*4 kappa,mu,mu2,mupr 
	real*8 mus,lams,sfac,flp,s0
	real*8 rho,g0,bigg,pi
	real*8 tau1,tau2
	real*8 a(6,6) 
	real*8 fl21,r,biga,r2
	common/erad/bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
	common/gval/g0
	pi=3.14159265358979d0
	tau1=dble(mu(i0)/eta1(i0))
	tau2=dble(mu2(i0)/eta(i0))
	mus=dble(mu(i0))*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+
     &	dble(mu(i0))*s0/dble(eta(i0)))
	lams=dble(kappa(i0)-2.*mus/3.)
	biga=lams+2.d0*mus   
	fl21=dble(l*(l+1))
	flp=dble(l)
	r2=r*r
c	Include self-gravitation for both g-terms and G-terms.
c**
	rho=dble(dens(i0))
c	g0 at radius r is now given in common.
c	g0=980.d0*dble(bigr*1.d+5)*1.e-11
c	  g0=0.d0
	bigg=(6.672d-11)*(1.d+2)*bigr*bigr
c		write(6,*)'amat6: g0,bigg=',g0,bigg,'bigr=',bigr
	r2=r*r
	a(1,1)=-2.d0*lams/(r*biga)
	a(1,2)=1.d0/biga
	a(1,3)=fl21*lams/(r*biga)
	a(1,4)=0.d0
	a(1,5)=0.
	a(1,6)=0.
	a(2,1)=(4.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,1)-4.d0*rho*g0/r
	a(2,2)=-2.d0/r+(2.d0/r)*lams*a(1,2)
	a(2,3)=-fl21*(2.d0*(lams+mus)/r)/r+(2.d0/r)*lams*a(1,3)
     &	+rho*g0*fl21/r
	a(2,4)=fl21/r
	a(2,5)=(flp+1.d0)*rho/r
	a(2,6)=-rho
	a(3,1)=-1.d0/r
	a(3,2)=0.d0
	a(3,3)=1.d0/r
	a(3,4)=1.d0/mus
	a(3,5)=0.
	a(3,6)=0.
	a(4,1)=-(lams/r)*a(1,1)+2.d0*(mus-biga)/r2+rho*g0/r2 
	a(4,2)=-(lams/r)*a(1,2)
	a(4,3)=-(lams/r)*a(1,3)+(fl21*biga-2.d0*mus)/r2
	a(4,4)=-3.d0/r
	a(4,5)=-rho/r
	a(4,6)=0.
	a(5,1)=4.*pi*bigg*rho
	a(5,2)=0.
	a(5,3)=0.
	a(5,4)=0.
	a(5,5)=-(flp+1.d0)/r
	a(5,6)=1.d0
	a(6,1)=4.*pi*bigg*rho*(flp+1.d0)/r
	a(6,2)=0.
	a(6,3)=-4.*pi*bigg*rho*fl21/r
	a(6,4)=0.
	a(6,5)=0.
	a(6,6)=(flp-1.d0)/r
c		write(6,*)'a-matrix='
c		write(6,*) (a(1,k1), k1=1,6)
c		write(6,*) (a(2,k1), k1=1,6)
c		write(6,*) (a(3,k1), k1=1,6)
c		write(6,*) (a(4,k1), k1=1,6)
c		write(6,*) (a(5,k1), k1=1,6)
c		write(6,*) (a(6,k1), k1=1,6)
c	  write(6,*)'A: aj=',aj
	return
	end 
	subroutine prod(a,c,n,b)
c	Form matrix product b=(A)c
	real*8 a,b,c 
	dimension a(n,n),c(n),b(n)
	do 5 i=1,n
	b(i)=0.d0
	do 4 j=1,n
4	b(i)=b(i)+a(i,j)*c(j)
5	continue
	return
	end 
	subroutine rzextr(iest,xest,yest,yz,dy,nv,nuse)
	implicit real*8 (a-h,o-z)
	parameter (imax=14,nmax=10,ncol=7)
c	The next line fixes an inaccuracy in Numerical Recipes,
c	which has the xval-array (actually named x) and d-array
c	in regular dimension statements.  They should rather be in common
c	because their values are needed in subsequent calls of rzextr.
	common/xvals/xval(imax),d(nmax,ncol)
	dimension yest(nv),yz(nv),dy(nv),fx(ncol)
	xval(iest)=xest
c		write(6,*)'A xval(',iest,')=',xval(iest)
	if(iest.eq.1) then
	do 11 j=1,nv
		  yz(j)=yest(j)
	 	  d(j,1)=yest(j)
		  dy(j)=yest(j)
11	continue
	else
		m1=min(iest,nuse)
c		write(6,*)'m1=',m1
	do 12 k=1,m1-1
c		write(6,*)'B xval(',iest-k,')=',xval(iest-k)
		fx(k+1)=xval(iest-k)/xest
c		write(6,*)'fx(',k+1,')=',fx(k+1)
12	continue
	do 14 j=1,nv
		yy=yest(j)
c		write(6,*)'in do 14, yy=',yy
		v=d(j,1)
c		write(6,*)'in do 14, v=',v
		c=yy
		d(j,1)=yy
		do 13 k=2,m1
		  b1=fx(k)*v
		  b=b1-c
		  if(b.ne.0.d0) then
			b=(c-v)/b
			ddy=c*b
			c=b1*b
		  else
			ddy=v
		  endif
		  if(k.ne.m1) v=d(j,k)
		  d(j,k)=ddy
		  yy=yy+ddy
13	continue
		dy(j)=ddy
c		write(6,*)'RZEXTR: dy(',j,')=',dy(j)
		yz(j)=yy
14	continue
	endif
	return
	end

