	program decay
c	Find singular values s=-s(j) of N-layer earth model
c	for toroidal modes.  Input file
c	'earth.model' has format
c	N  Nd 	(Number of layers; Number of different layers)
c	rb(1)	     |  rt(1)   |dens(1)|  kappa(1)  | mu(1)  | eta(1)
c	...		...	 ...	   ...	       ...      ...
c	rb(N)	     |  rt(N)   |dens(N)|  kappa(N)  | mu(N)  | eta(N)
c	bottom radius|top radius|density|bulk modulus|rigidity|viscosity
c	(km)          (km)      (g/cm**2) (10**11 dyne/cm**2)  10**19
c							       dyne-s/cm**2
c							       10**18 Pa-s
	character*80 aread
	dimension dens(200)
	real*4 kappa(200),mu,mu2,mupr
	dimension si(600)
	real*8 aj(2,2),zcr(100),b(2),c(2),s0,s1(29),sfrac(29)
	dimension nfrac(29),tcr0(29),tcr1(29),iter(29)
	common/mat/rb(200),rt(200),mu(200),mu2(200),eta(200),
     &	eta1(200),l,s0,mupr(200)
	common/ifrdum/ifr
	pi=3.14159265
c	bigr=6371.
c	bigr=earth's radius in km.
c
c	Read in earth model
	open(4,file='earth.model')
	rewind(4)
	read(4,5) n,nd,bigr,depfac 
	write(6,5) n,nd,bigr,depfac
5	format(2i2,2f10.3)
	nm1=2*nd
c ** *
	  etala=1.e+13
	do 10 j=1,n
	i=2*j-1
c
	read(4,14) aread 
	read(aread,18,err=17) rb(i),rt(i),dens(i),kappa(i),mu(i),
     &	mupr(i),eta(i),eta1(i)
	  etaval=min(eta(i),eta1(i))
	if(eta1(i).eq.0.0) go to 17
	go to 12
17	eta1(i)=1.e+13
	read(aread,16,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),
     &	mupr(i),eta(i)
	  etaval=eta(i)
	if(eta(i).eq.0.0) go to 11
	go to 12
11	eta1(i)=1.e+13
	mupr(i)=0.
	read(aread,15) rb(i),rt(i),dens(i),kappa(i),mu(i),eta(i)
	  etaval=eta(i)
12	continue
	mu2(i)=mupr(i)*mu(i)/(mu(i)-mupr(i))
	write(6,18) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i),eta1(i)
c***	Above lines: specify [mu prime] for standard linear solid
c***	in same units as mu.
c***	mu prime = 0 corresponds to Maxwell material.
c***	Reference: Cohen, Time dependent deformation following an earthquake,
c***	JGR, vol. 87, pp. 5414-5421 (1982).
c***	eta1 is an additional viscosity for a Burghers solid.
c***	The following correspondences apply:
c***	mu == mu sub 1, mupr==mu', eta==eta sub 2, eta1==eta sub 1.
c***	mu'=(mu sub 1)*(mu sub 2)/(mu sub 1 + mu sub 2)
c***	where mu sub 1, eta sub 1, mu sub 2, eta sub 2 are as described
c***	by Peltier (1985 Science, v. 318, p. 615, eqn 1).
c--	Determine elastic plate thickness.
	if(etaval.ge.1.e+10.and.etala.lt.1.e+10) helast=rb(i)
	etala=etaval
c--
	dens(i+1)=dens(i)
	mu(i+1)=mu(i)
	mu2(i+1)=mu2(i)
	mupr(i+1)=mupr(i)
	eta(i+1)=eta(i)
	eta1(i+1)=eta1(i)
	kappa(i+1)=kappa(i)
	rb(i+1)=(rb(i)+rt(i))/2.
	rt(i+1)=rt(i)
	rt(i)=rb(i+1)
	rb(i+1)=rb(i+1)/bigr + (1.-6371./bigr)
	rt(i+1)=rt(i+1)/bigr + (1.-6371./bigr)
	rb(i)=rb(i)/bigr + (1.-6371./bigr)
	rt(i)=rt(i)/bigr + (1.-6371./bigr)
10	continue
c		pause
14	format(a80)
15	format(5f9.3,e13.6e2)
16	format(6f9.3,e13.6e2)
18	format(6f9.3,2e13.6e2)
	n=2*n 
	close(4)
c--
	  helast=6371.-helast
	  write(6,*)'elastic plate thickness=',helast
c--
c ** *
c	Since the characteristic relaxation time is equal
c	to eta/(mu + mub) rather than eta/mu, rescale eta accordingly.
c	Note that mub=mupr*mu/(mu-mupr).
c---	do 19 i=1,n
c---	eta(i)=eta(i)*(1.-mupr(i)/mu(i))
c---19	continue
c***
	do 20 i=1,n
	si(3*i-2)=(mu(i)+mu2(i))/eta(i)
	si(3*i-1)=(mu(i)/eta(i))
20	si(3*i)=mu(i)/eta1(i)
	m=3*n
	call order(si,m)
		write(6,*)'si=',si
	open(2,file='decay.out')
	write(2,21) nd,bigr,depfac
21	format(i2,3f10.3)
	write(6,*)'minimum l, maximum l?'
	read(5,*) lmin,lmax 
	do 50 l=lmin,lmax 
	write(6,*)'l=',l
	depmax=1.-helast/bigr-3.0*(2.*pi/(real(l)+0.5))
	depmin=1.-helast/bigr+1.5*(2.*pi/(real(l)+0.5))
		dwrite=bigr-bigr*depmax
		write(6,*)'start integration at depth=',dwrite,' km'
		dwrite=bigr-bigr*depmin
		write(6,*)'top integration at depth=',dwrite,' km'
c	depmax=dimensionless radius 4.0 wavelengths into the earth.
C	Look for zero crossings in surface traction boundary condition.
c	Start with partition of smin -- smax into 100 units.
c	Subsequent zero crossing intervals zcr1-zcr2 are divided
c	into 3 units.
	npole=0
	nzero=0
24	nfrac(1)=1600
c	if(npole.eq.1) nfrac(1)=800
	  sidif=si(npole+2)-si(npole+1)
c	sfrac(1)=dble(sidif)/(dble(nfrac(1))-0.96d0)
c	 s1(1)=-sfrac(1)*0.98d0+dble(si(npole+1))
        sfrac(1)=dble(sidif)/(dble(nfrac(1))-0.996d0)
	 s1(1)=-sfrac(1)*0.998d0+dble(si(npole+1))
	iter(1)=0
	ifr=1
	tcr0(ifr)=0.
25	iter(ifr)=iter(ifr)+1
	if(iter(ifr).gt.nfrac(ifr)) go to 45
	s0=-s1(ifr)-dble(iter(ifr))*sfrac(ifr)
c	Start matrix propagation.
c	b(1)=1.d0
c	b(2)=0.d0
	is=0
	do 30 i=1,n
	if(i.lt.2) go to 30
c	Start integration at most 2.5 wavelengths into the earth.
	if(rb(i).lt.depmax) go to 30
	if(rt(i).gt.depmin) go to 30
c	Give boundary conditions at top of layer (i-1) appropriate for
c	homogeneous solid in 0 < r < rt(i-1).
	if(is.ne.0) go to 28
	is=1
	j=i-1
	call matra(aj,j,1)
	b(1)=aj(1,1)
	b(2)=aj(2,1)
28	continue
c	Form matrix aj in layer i, evaluated at bottom radius.
	call matra(aj,i,0)
c***	Check against mus too close to zero.
	call ainver(c,aj,b,2)
	call matra(aj,i,1)
	call prod(aj,c,2,b)
	bmag=abs(real(b(1)))
	if(bmag.lt.1.e+10) go to 30
	b(1)=b(1)*1.d-10
	b(2)=b(2)*1.d-10
30	continue
	tcr1(ifr)=real(b(2))
c		write(6,*)'s0=',s0
c		write(6,*)'previous/current traction=',tcr0(ifr),tcr1(ifr)
	prodx=tcr0(ifr)*tcr1(ifr)
	if(tcr1(ifr).eq.0.) go to 35
	if(prodx.lt.0.) go to 35
	tcr0(ifr)=tcr1(ifr)
	go to 25
35	ifr=ifr+1
	if(ifr.gt.21) go to 40
	sfrac(ifr)=sfrac(ifr-1)/3.d0
	s1(ifr)=-(s0+sfrac(ifr-1))
	tcr0(ifr)=tcr0(ifr-1)
	iter(ifr)=0
	nfrac(ifr)=3
	go to 25
40	nzero=nzero+1
	zcr(nzero)=s1(21)
c		write(6,*)'line 40: s0=',s0,'s1(21)=',s1(21)
c--
	rat=abs(b(2)/b(1))
	write(6,*)'stress/displacement ratio=',rat
		rat=0.01
	if(rat.lt.0.1) then
	write(2,41) l,zcr(nzero)
	tyr=3.17/real(s1(21))
	endif
	write(6,*)'decay time #',nzero,'=',tyr,'years'
c--
41	format(i4,d17.10)
	if(nzero.eq.nm1) go to 50
	ifr=1
	tcr0(1)=tcr1(1)
	go to 25
45	write(6,*)'no more zero crossings: npole=',npole
	npole=npole+1
	if(npole.lt.(2*nd)) go to 24
50	continue
	end 
	subroutine matra(aj,i,n)
c	Compute matrix elements for toroidal modes (decay times).
c	n=0: evaluate at bottom of layer i
c	n=1: evaluate at top of layer i
	real*4 mu,mu2,mupr
	real*8 mus,s0,tau1,tau2
	common/mat/rb(200),rt(200),mu(200),mu2(200),eta(200),
     &	eta1(200),l,s0,mupr(200)
	common/ifrdum/ifr
	real*8 aj(2,2)
	tau1=dble(mu(i)/eta1(i))
	tau2=dble(mu2(i)/eta(i))
	mus=mu(i)*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
c		if(s0.lt.-466.d0) then
c		write(6,*)'mu2(',i,')=',mu2(i),'eta=',eta(i)
c		write(6,*)'s0=',s0,'tau2=',tau2
c		write(6,*)'mus=',mus,'s0+tau2=',s0+tau2
c		write(6,*)'-----------------'
c		endif
c		if(ifr.eq.21) then
c		write(6,*)'mus(',i,')=',mus
c		write(6,*)'eta1,eta=',eta1(i),eta(i)
c		write(6,*)'tau1,tau2=',tau1,tau2
c		write(6,*)'denominator=',((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
c		write(6,*)'------------'
c		endif
c
	if(n.ne.0) r=rt(i)
	if(n.eq.0) r=rb(i)
	aj(1,1)=dble((r/rb(i))**l)
	aj(1,2)=dble((rb(i)/r)**(l+1))
	aj(2,1)=mus*dble(l-1)*dble((r/rb(i))**(l-1)*(1./rb(i)))
	aj(2,2)=-mus*dble(l+2)*dble((rb(i)/r)**(l+2)*(1./rb(i)))
	return
	end 
	subroutine prod(a,c,n,b)
c	Form matrix product b=(A)c
	real*8 a(n,n),c(n),b(n)
	do 5 i=1,n
	b(i)=0.d0
	do 4 j=1,n
4	b(i)=b(i)+a(i,j)*c(j)
5	continue
	return
	end 
	subroutine ainver(a,b,c,n)
c       Find solution (a) of the matrix equation (B)(a)=(c).
c       Method:  Gaussian elimination.
	integer*2 perm 
	dimension perm(10)
	real*8 a(n),b(n,n),c(n),bsave,fac
        do 5 i=1,n
5       perm(i)=i
	i=0
10      i=i+1
	if (i.gt.n) go to 35
c       Find maximum in row i.
	amax=0.
	imax=i
	do 15 j=i,n
	t=real(b(i,j))
	if (t.lt.0.) t=-t
	if (t.lt.amax) go to 15
	amax=t
	imax=j
15      continue
	j=imax
c       Switch columns i and j.
	do 20 m=1,n
	bsave=b(m,i)
	b(m,i)=b(m,j)
20      b(m,j)=bsave
	iperm=perm(i)
	perm(i)=perm(j)
	perm(j)=iperm
c       Eliminate ith column.
c	  write(6,*)'pivot row',i,'=',b(i,i)
c	  if(i.eq.n) pause  
	do 25 j=1,n
	if (j.eq.i) go to 25
	fac=b(j,i)/b(i,i)
	do 30 k=i,n
30      b(j,k)=b(j,k)-fac*b(i,k)
	c(j)=c(j)-fac*c(i)
25      continue
	go to 10
35      do 40 i=1,n
	k=perm(i)
40      a(k)=c(i)/b(i,i)
	return
	end
	subroutine order(y,n)
c	Rearrange elements of positive-valued array y
c	in order of increasing values.  n.le.200.
c	Only recognize distinct layers -- layers of equal
c	material properties are all considered as one.
	dimension y(600)
	do 10 i=1,n-1
	do 15 j=i+1,n
	ymag=abs(y(j)-y(i))
	if(ymag.lt.1.e-06) go to 11
	if(y(j).gt.y(i)) go to 15
	if(y(j).lt.y(i)) go to 14
11	y(j)=1.e+12
	go to 15
14	continue 
	ysav=y(i)
	y(i)=y(j)
	y(j)=ysav
15	continue
10	continue 
	return 
	end 


