	program decay4m
c**
c	OCT 9, 1996.  Same as DECAY4 but gravitation is included
c	and method of second order minors is used.
c**
c	Find singular values s=-s(j) of N-layer earth model
c	for spheroidal modes.
c       Input file
c	'earth.model' has format
c	N  Nd 	(Number of layers; Number of different layers)
c	rb(1)	     |  rt(1)   |dens(1)|  kappa(1)  | mu(1)  | eta(1)
c	...		...	 ...	   ...	       ...      ...
c	rb(N)	     |  rt(N)   |dens(N)|  kappa(N)  | mu(N)  | eta(N)
c	bottom radius|top radius|density|bulk modulus|rigidity|viscosity
c	(km)          (km)      (g/cm**3) (10**11 dyne/cm**2)  10**19
c							       dyne-s/cm**2
c							       10**18 Pa-s
	character*80 aread
	real*4 kappa,mu,mu2,mupr 
	dimension si(600)
	real*8 zcr(100)
	real*8 b1(5),s0,s1(29),sfrac(29)
		real*8 b1sav(5),htot,eps,rsav,hdid,hnext,yscal(5),dydr(5)
		real*8 fl
	real*8 aj4(4,4)
	real*8 b14(4),b24(4)
	real*8 r,dr
	dimension nfrac(29),tcr0(29),tcr1(29),iter(29)
	common/erad/bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
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
	nm1=7*nd
	n2m1=2*nd-1
	n3m3=3*nd-3
c ** *
cNEW03
	relast=0.0
c
	do 10 j=1,n
	i=2*j-1
c
	read(4,14) aread 
	read(aread,18,err=17) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i),eta1(i)
	if(eta1(i).eq.0.0) go to 17
	go to 12
17	eta1(i)=1.e+13
	read(aread,16,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i)
	if(eta(i).eq.0.0) go to 11
	go to 12
11	eta1(i)=1.e+13
	mupr(i)=0.
	read(aread,15) rb(i),rt(i),dens(i),kappa(i),mu(i),eta(i)
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
c	  write(6,*)'i=',i,'rt(i),rt(i+1)=',rt(i),rt(i+1)
10	continue
14	format(a80)
15	format(5f9.3,e13.6e2)
16	format(6f9.3,e13.6e2)
18	format(6f9.3,2e13.6e2)
	n=2*n
	n3m0=3*n 
	close(4)
c		write(6,*)'A: bigr=',bigr
c ** *
c	Since the characteristic relaxation time is equal
c	to eta/(mu + mub) rather than eta/mu, rescale eta accordingly.
c	Note that mub=mupr*mu/(mu-mupr).
c----	do 19 i=1,n
c----19	eta(i)=eta(i)*(1.-mupr(i)/mu(i))
c***
c	do 20 i=1,n
c	si(3*i-2)=(mupr(i)/eta(i))
c	si(3*i-1)=mu(i)/eta(i)
cc	si(2*i)=(mu(i)/eta(i))*kappa(i)/(kappa(i)+4.*mu(i)/3.)
c	si(3*i)=(mu(i)/eta(i))*(kappa(i)-2.*mupr(i)/3.)/
c&	(kappa(i)-2.*mu(i)/3.)
cc	si(2*i)=(mu(i)/eta(i))*(kappa(i)-2.*mu(i)/3.)/(kappa(i)+4.*mu(i)/3.)
c20	continue 
	do 20 i=1,n
	si(3*i-2)=(mu(i)+mu2(i))/eta(i)
	si(3*i-1)=(mu(i)/eta(i))
20	si(3*i)=mu(i)/eta1(i)
	call order(si,n3m0)
	if(nd.eq.1) si(1)=0.
	  write(6,*)'si=',si
c	Look for zeros outside of maximum mu/eta if necessary.
	open(2,file='decay4.out')
	write(2,21) nd,bigr,depfac
21	format(i2,2f10.3)
	write(6,*)'minimum/maximum l?'
	read(5,*) lmin,lmax
	do 50 l=lmin,lmax 
	write(6,*)'l=',l
	depmax=1.-3.0*(2.*pi/(real(l)+0.5))
cTEST
c		depmax=1.-200./6371.
	depmin=1.
c-------------
cNOTE:	Subroutine GETMIN assumes that depmin=1.0, i.e., Earth's surface.
c	In current applications this is always the case.  However, the current
c	version of VSPHM would not work correctly if depmin were < 1.  
c	Make sure it always equals 1!
c-------------
		dwrite=bigr-bigr*depmax
		write(6,*)'start depth of integration=',dwrite,'km'
c	depmax=dimensionless radius 2 wavelengths into the earth.
c	Look for zero crossings in surface traction boundary condition.
c	Start with partition of smin -- smax into 100 units.
c	Subsequent zero crossing intervals zcr1-zcr2 are divided
c	into 3 units.
	npole=0
	nzero=0
	fl21=real(l*(l+1))
	flp1=real(l+1)
24	nfrac(1)=1600
c	if(npole.eq.1) nfrac(1)=200
c	if(npole.eq.2) nfrac(1)=200
c		if(npole.gt.0) nfrac(1)=200
c		if(nzero.gt.1) nfrac(1)=50
c		if(nzero.gt.2) nfrac(1)=25
c	if(l.gt.600.and.npole.eq.2) nfrac(1)=400
	  sidif=si(npole+2)-si(npole+1)
        sfrac(1)=dble(sidif)/(dble(nfrac(1))-0.98d0)
	 s1(1)=-sfrac(1)*0.99d0+dble(si(npole+1))
	iter(1)=0
	ifr=1
	tcr0(ifr)=0.
25	iter(ifr)=iter(ifr)+1
	if(iter(ifr).gt.nfrac(ifr)) go to 45
	s0=-s1(ifr)-dble(iter(ifr))*sfrac(ifr)
c		if(-3.17d0/s0.gt.39.90d0.or.-3.17d0/s0.lt.17.8d0) go to 25
c	Start matrix propagation.
	is=0
	do 30 i=1,n
c	if(i.lt.2) go to 30
c	Start integration at most 2 wavelengths into the earth.
	if(rb(i).lt.depmax) go to 30
cc	Give boundary conditions at top of layer (i-1) appropriate
cc	for rigid solid in 0 < r < rt(i-1).
c	Give boundary conditions at top of layer (i-1) appropriate for
c	homogeneous solid in 0 < r < rt(i-1).
	if(is.ne.0) go to 28
	is=1
c	b14(1)=0.d0
c	b14(2)=1.d0
c	b14(3)=0.d0
c	b14(4)=0.d0
c	b24(1)=0.d0
c	b24(2)=0.d0
c	b24(3)=0.d0
c	b24(4)=1.d0
	j=i
		if(ifr.eq.25) write(6,*)'calling matra at r=',rb(j)*6371.,'km'
	call matra(aj4,j,0)
	do 27 k=1,4
	b14(k)=aj4(k,1)
	b24(k)=aj4(k,2)
27	continue
c	form second order minors corresponding to chosen b's above.
	b1(1)=b14(1)*b24(2)-b14(2)*b24(1)
	b1(2)=b14(1)*b24(3)-b14(3)*b24(1)
	b1(3)=b14(1)*b24(4)-b14(4)*b24(1)
	b1(4)=b14(2)*b24(3)-b14(3)*b24(2)
	b1(5)=b14(2)*b24(4)-b14(4)*b24(2)
		if(ifr.eq.25) write(6,*)'initial b14=',b14,'b24=',b24
28	continue
cc*
c	nrad=int(sqrt(real(l/10)))+1
	  nrad=2
	dr=dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions up to
c	the surface.
	r=dble(rb(i))
c---
c	do 31 j=1,nrad
c	call amat5(aj,i,r)
cc	New b(r+dr)=exp(aj*(dr)).b(r)
c	call matexp(aj,5,dr,ak)
c	call prod(ak,b1,5,c1) 
c	call equal(c1,b1,5)
c32	r=r+dr    
c31	continue 
c---
	b1sav(1)=b1(1)
	b1sav(2)=b1(2)
	b1sav(3)=b1(3)
	b1sav(4)=b1(4)
	b1sav(5)=b1(5)
	fl=dble(l)
36	continue
	r=dble(rb(i))
c	Do b1 using bsstep.
c		write(6,*)'nrad=',nrad
	iwk=0
c	if(dr.eq.0.0) go to 31
	do 31 j=1,nrad
	rsav=r
c		write(6,*)'entering derivs5'
	call derivs5(i,r,b1,dydr)
	htot=dr
	eps=3.d-7
	yscal(1)=dabs(b1(1))+fl*dabs(b1(3))+dabs(b1(2))*fl**2+dabs(b1(4))*fl
     &		 +dabs(b1(5))
	yscal(2)=yscal(1)/fl**2
	yscal(3)=yscal(1)/fl
	yscal(4)=yscal(3)
	yscal(5)=yscal(1)
	call bsstep(b1,dydr,5,r,htot,eps,i,yscal,hdid,hnext)
	if((hdid-htot)/htot.ne.0.d0) iwk=1
c		write(6,*)'j=',j,'r=',r,'b1(5)=',b1(5)
	r=rsav+dr    
31	continue 
	if(iwk.ne.0) then
		write(6,*)'*** iwk=1 ***: redoing layer ',i
		b1(1)=b1sav(1)
		b1(2)=b1sav(2)
		b1(3)=b1sav(3)
		b1(4)=b1sav(4)
		b1(5)=b1sav(5)
		write(6,*)'old nrad=',nrad
		nrad=2*nrad
		write(6,*)'new nrad=',nrad
		dr=dr/2.d0
		go to 36
	endif
c---
	bmag=abs(real(b1(1)))+real(l)*abs(real(b1(3)))
	if(bmag.lt.1.e+5) go to 39
	do 37 j1=1,5
37	b1(j1)=b1(j1)*1.d-5
39	continue
c		if(n.gt.0) stop
30	continue
	tcr1(ifr)=real(b1(5))
c	  write(6,*)'s0=',s0,tcr0(ifr),tcr1(ifr)
c	  write(6,*)'s0=',s0,'previous/current traction=',
c     &	  tcr0(ifr),tcr1(ifr) 
	  if(tcr1(ifr).eq.0.) s1(21)=-s0
	  if(tcr1(ifr).eq.0.) go to 40
c	if(prodx.lt.0.) go to 35
	if(tcr1(ifr).lt.0.0.and.tcr0(ifr).gt.0.0) go to 35
	if(tcr0(ifr).lt.0.0.and.tcr1(ifr).gt.0.0) go to 35
	tcr0(ifr)=tcr1(ifr)
	go to 25
35	ifr=ifr+1
	if(ifr.gt.25) go to 40
	sfrac(ifr)=sfrac(ifr-1)/3.d0
	s1(ifr)=-(s0+sfrac(ifr-1))
	tcr0(ifr)=tcr0(ifr-1)
	iter(ifr)=0
	nfrac(ifr)=3
	go to 25
40	nzero=nzero+1
	zcr(nzero)=s1(21)
	write(2,53) l,zcr(nzero),depmax,depmin
	write(6,53) l,zcr(nzero),depmax,depmin
41	format(i4,d17.10)
cOLD53	format(i4,d22.15,2f10.6)
53	format(i4,3d22.15)
	tyr=3.17/real(s1(21))
	write(6,*)'decay time #',nzero,'=',tyr,'years'
	write(6,*)'b1=',b1
	if(nzero.eq.nm1) go to 50
93	ifr=1
	tcr0(1)=tcr1(1)
	go to 25
45	write(6,*)'no more zero crossings: npole=',npole
	npole=npole+1
	if(npole.le.n3m3) go to 24
50	continue
	end
	subroutine matra(aj,i,n)
c	Compute matrix elements for spheroidal modes (decay times).
c	n=0: evaluate at bottom of layer i
c	n=1: evaluate at top of layer i
	real*4 mu,mu2,kappa,mupr
	real*8 mus,lams,sfac,biga,r,flp,flm,r2,flp3,flm3,s0
	real*8 tau1,tau2
	real*8 a1p,a1m,a2p,a2m,a3p,a3m,a4p,a4m
	real*8 b1p,b1m,b2p,b2m,b3p,b3m,b4p,b4m,rlp,rlm  
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0
	real*8 aj(4,4)
	tau1=dble(mu(i)/eta1(i))
	tau2=dble(mu2(i)/eta(i))
	mus=mu(i)*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
c		write(6,*)'MATRA: tau1,tau2,mus=',tau1,tau2,mus,'r=',r
	lams=dble(kappa(i)-2.*mus/3.)
	biga=lams+2.d0*mus   
10	if(n.ne.0) r=dble(rt(i))
	if(n.eq.0) r=dble(rb(i))
c	Determine matrix elements (Takeuchi and Saito, 1972, p. 243-244).
	flp=dble(l)
	flm=-flp-1.d0
	flp3=2.d0*(2.d0*flp+3.d0)
	flm3=2.d0*(2.d0*flm+3.d0)
	r2=r*r
	aj(1,1)=flp*(flp+1.d0)
	aj(1,3)=flm*(flm+1.d0)
	aj(2,1)=2.d0*mus*(flp*(flp*flp-1.d0))
	aj(2,3)=2.d0*mus*(flm*(flm*flm-1.d0))
	aj(3,1)=(flp+1.d0)
	aj(3,3)=(flm+1.d0)
	aj(4,1)=2.d0*mus*(flp*flp-1.d0)
	aj(4,3)=2.d0*mus*(flm*flm-1.d0)
	a1p=-(flp+2.d0)*r2/flp3
	a1m=-(flm+2.d0)*r2/flm3
	a2p=-(lams+2.d0*mus*(flp+2.d0)*(flp+1.d0)/flp3)*r2
	a2m=-(lams+2.d0*mus*(flm+2.d0)*(flm+1.d0)/flm3)*r2
	a3p=-r2/flp3
	a3m=-r2/flm3
	a4p=-2.d0*mus*(flp+1.d0)*r2/flp3
	a4m=-2.d0*mus*(flm+1.d0)*r2/flm3
	b1p=flp*(flp+1.d0)*r2/flp3
	b1m=flm*(flm+1.d0)*r2/flm3
	b2p=2.d0*mus*flp*(flp+1.d0)*(flp+1.d0)*r2/flp3
	b2m=2.d0*mus*flm*(flm+1.d0)*(flm+1.d0)*r2/flm3
	b3p=(flp+3.d0)*r2/flp3
	b3m=(flm+3.d0)*r2/flm3
	b4p=mus*2.d0*flp*(flp+2.d0)*r2/flp3
	b4m=mus*2.d0*flm*(flm+2.d0)*r2/flm3
	aj(1,2)=mus*(flp+1.d0)*a1p+biga*b1p
c		write(6,*)'MATRA2: aj(1,2)=',aj(1,2)
	aj(1,4)=mus*(flm+1.d0)*a1m+biga*b1m
	aj(2,2)=mus*(flp+1.d0)*a2p+biga*b2p
c		write(6,*)'lams,mus,flp=',lams,mus,flp
c		fdum=real(2.d0*mus*(flp+1.d0)*(lams*(flp*flp-flp-3.d0)+
c	&		mus*(flp*flp-flp-2.d0)))/real(2.d0*(2.d0*flp+3.d0))
c		write(6,*)'aj(2,2)again=',fdum
	aj(2,4)=mus*(flm+1.d0)*a2m+biga*b2m
	aj(3,2)=mus*(flp+1.d0)*a3p+biga*b3p
	aj(3,4)=mus*(flm+1.d0)*a3m+biga*b3m
	aj(4,2)=mus*(flp+1.d0)*a4p+biga*b4p
	aj(4,4)=mus*(flm+1.d0)*a4m+biga*b4m 
c	Normalize matrix elements.
	rlp=(r/dble(rb(i)))**l
	rlm=(dble(rb(i))/r)**(l+1)
	do 20 m=1,2
	aj(1,m)=aj(1,m)*rlp/r
	aj(1,m+2)=aj(1,m+2)*rlm/r 
	aj(2,m)=aj(2,m)*rlp/r2
	aj(2,m+2)=aj(2,m+2)*rlm/r2
	aj(3,m)=aj(3,m)*rlp/r
	aj(3,m+2)=aj(3,m+2)*rlm/r 
	aj(4,m)=aj(4,m)*rlp/r2
	aj(4,m+2)=aj(4,m+2)*rlm/r2
20	continue
c		write(6,*)'MATRA2A: aj(1,2)=',aj(1,2)
	return
	end 
	subroutine order(y,n)
c	Rearrange elements of positive-valued array y
c	in order of increasing values.  n.le.20.
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
	subroutine equal(a,b,n)
c	Set the n x 1 matrix (b) equal to (a).
	real*8 a(n),b(n)
	do 10 j=1,n
10	b(j)=a(j)
	return
	end 
