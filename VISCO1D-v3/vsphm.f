	program vsphm
c	NOTE: the problem is that go 30 goes all the way
c	down to layer 1.  It should stop at the radius level corresponding
c	to radmin
c*
c*
c	NOTE: 10.29.98 modified do 52 to write out k=1-500
c
c	Determine eigenfunctions using l-values determined in RWAVESM.
c	---------------------------
c	UP TO 125 modes and 250 frequencies.
c	Sequences of (l,n,f,ur,dur,vr,dvr) are written into
c	file 'reigen.out', where l=total degree, n=radial order number,
c	f=frequency in rad/sec, ur,dur,vr, and dvr are normalized eigenfunctions
c	evaluated at top 30 radius points as given by the dimensionless
c	array rw(500).  Eigenfunctions written out in cgs units and
c	are normalized according to
c	INT {r=0 to r=a} w0**2 [rho*r**2*{U(r)**2+l(l+1)V(r)**2}] dr = 1. 
c	Input file   
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
	real*4 kappa,mu,mu2,mupr,mus,lams 
c	real*4 deter 
	real*8 mudum,lamdum,denss
	real*8 s0,zwrite
c	real*8 zcr(125)
c------
cNOTE
	parameter (maxmod=21)
	parameter (lmax0=2500)
c	NOTE: do not expect more than 19 mode branches
c	or more than 2500 l's.
c------
	real*8 htot,eps,rsav,hdid,hnext,yscal(4),dydr(4)
	real*8 valt
	real*8 r,r0,dr
		real*8 fl
	real*8 a1b,a2b,a3b,a4b,c5
	real*8 bmf1,bmf2,bmf3,bmf4,bmf5,bmftop
	real*8 b3(4),bm(5)
	dimension mdec(lmax0),rdep1(lmax0,maxmod),rdep2(lmax0,maxmod)
c	real*4 dens1(200),kapp1(200),mu1(200)
	common/botcof/a1b(2),a2b(2),a3b(2),a4b(2),c5
c---
cOLD	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),l,w0,r0
	common/erad/bigr
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0,r0
c---
	common/jfac0/j01,j02,jm1
	common/jdif/jdif1,jdif3,jdif4,jdif6
	common/jfac/jfac1,jfac2,jfac3,jfac4
c*
	real*8 sdecr,sdec(maxmod,lmax0)
	real*8 ya1,ya2,ya3,ya4
	parameter (ndep=50)
	common/upd1/rw(500),ur(500),dur(500),vr(500),dvr(500),densi(500),
     &	sur(500)
	common/upd/ya1(51),ya2(51),ya3(51),ya4(51)
c*
	common/minfn/bmf1(201,51),bmf2(201,51),bmf3(201,51),
     &	bmf4(201,51),bmf5(201,51),bmftop
c*
	pi=3.14159265
c----------------
cnew
c	Read in earth model
	open(4,file='earth.model')
	rewind(4)
	read(4,5) n,nd,bigr,depfac 
	write(6,5) n,nd,bigr,depfac 
5	format(2i2,2f10.3)
c ** *
	do 10 j=1,n
	i=2*j-1
c
	read(4,14) aread 
	read(aread,18,err=27) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i),eta1(i)
	if(eta1(i).eq.0.0) go to 27
	go to 12
27	eta1(i)=1.e+13
	read(aread,17,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i)
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
		write(6,*)'eta1(',i,')=',eta1(i)
		write(6,*)'eta1(',i+1,')=',eta1(i+1)
	kappa(i+1)=kappa(i)
	rb(i+1)=(rb(i)+rt(i))/2.
	rt(i+1)=rt(i)
	rt(i)=rb(i+1)
	rb(i+1)=rb(i+1)/bigr + (1.-6371./bigr)
	rt(i+1)=rt(i+1)/bigr + (1.-6371./bigr)
	rb(i)=rb(i)/bigr + (1.-6371./bigr)
	rt(i)=rt(i)/bigr + (1.-6371./bigr)
10	continue
14	format(a80)
15	format(5f9.3,e13.6e2)
17	format(6f9.3,e13.6e2)
18	format(6f9.3,2e13.6e2)
	n=2*n 
c ** *
	close(4)
c--------------
c	vsph20.out will have eigenfunctions written out
c	over some depth range depmax to depmin km.  This info is
c	not needed in strains/strainw/strainA, which reads only vsph.out.
c	So just pick a range here.
	depmax=250.
	depmin=0.
	ddep=(depmax-depmin)/real(ndep-1)
	do ir=1,500
	rw(ir)=1.-(depmax-real(ir-1)*ddep)/bigr
	enddo
c
        do idep=1,40
        de=28.+(0.001 - 28.)*real(idep-1)/39.
        rw(idep+ndep)=1.-de*depfac/bigr
        enddo
c*
	write(6,*)'depth level of extra output Greens function (km)?'
	read(5,*) edep
	dep41=1.-edep/bigr
	rw(41+ndep)=dep41
	rw(42+ndep)=1.
c--------------
c**	Read in decay times already determined by DECAY4M.
	write(6,*)'*** read in decay times ***'
	open(4,file='decay4.out')
	rewind(4)
	read(4,21) ndd,bigrr,depfac
	if(ndd.ne.nd) write(6,*)'read incorrect #layers in decay4m.out'
	write(6,*)'ndd=',ndd,'nd=',nd
	if(ndd.ne.nd) stop 
	write(6,21) nd,bigrr,depfac
	idec=0
	is=0 
16	idec=idec+1 
cNEW
	read(4,53,end=20) l,sdecr,depmax,depmin
	if(is.eq.0) llast=l
	lmax=l
	if(l.eq.llast) mdec(l)=idec
	if(l.eq.llast) rdep1(l,mdec(l))=depmax
	if(l.eq.llast) rdep2(l,mdec(l))=depmin
c-----
	if(l.eq.llast) sdec(idec,l)=sdecr
	if(is.eq.0) lmin=l  
	if(is.eq.0) is=1
	if(l.eq.llast) go to 16 
	sdec(1,l)=sdecr
	idec=1
	mdec(l)=1
cNEW
	rdep1(l,mdec(l))=depmax
	rdep2(l,mdec(l))=depmin
c------
	llast=l 
	go to 16
20	close(4)
cOLD 53	format(i4,d22.15,2f10.6)
53	format(i4,3d22.15)
	write(6,*)'read decay times: lmin,lmax=',lmin,lmax
c**
c-------
cnew
	open(4,file='vsph20.out',form='unformatted')
	open(8,file='vsph.out',form='unformatted')
	write(4) nd
	write(8) nd
c	write(8) dep41
c----
21	format(i2,2f10.3)
	do 50 l=lmin,lmax
	fl=dble(l)
	write(6,*)'l=',l
	fl21=real(l*(l+1.))
c	Revise starting radius for numerical integration.
c	Put bottom of radius integration [fnwav] wavelengths down 
c	into earth, where fnwav=1.0*nzero+2.5  Thus fundamental mode
c	integration started at 2.5 wavelengths down, 1st overtone
c	started at 3.5 wavelengths down, etc.
c	radmin=rb(2)
cOLD	radmin=(1.-(depmax)/bigr) - 2.5*(2.*pi/(real(l)+0.5))
c	NOTE: above line used to TEST how well programs do with fixed radmin
c	(combined with homogeneous material beneath).
	write(6,*)'bottom of integration at dimensionless radius',radmin
	write(6,*)'pi=',pi,'l=',l
c	write(6,*)'rb(4),rb(5)=',rb(4),rb(5)
	do 48 nzero=1,mdec(l)
		write(6,*)'doing nzero=',nzero,' out of',mdec(l)
	radmin=rdep1(l,nzero)
	s0=-sdec(nzero,l)
	s0r=real(s0)
	do 22 k=1,500
	ur(k)=0.
	dur(k)=0.
	vr(k)=0.
	dvr(k)=0.
	densi(k)=0.
	sur(k)=0.
22	continue
	wo=0.
c	om2=dble(w0*w0)
c * *	Prescribe the minors for this (l,w0) combination from the free surface
c	down to 500 km depth.
		write(6,*)'calling getmin: l,s0=',l,s0
	call getmin(n,radmin)
		write(6,*)'bmftop=',bmftop
c		if(k.ge.0) stop
c	bmftop is the dimensionless ratio between the minor m5 and the
c	other minors at the free surface.  It should theoretically = zero.
c	Badly determined modes listed in 'DECAY4M.OUT' can be eliminated
c	by using bmftop as a diagnostic, as in the following lines.
c	if(bmftop.gt.1.d-4) then
c	  gpvel=1.e-3
c	  wo=1.
c	  go to 152
c	endif
c		if(iw0.eq.39.and.nzero.gt.1) pause
c * *
c		write(6,*)'out of getmin: bmf3(19,2)=',bmf3(19,2)
c	Start matrix propagation.
c 
47	continue
	do 30 i=n,2,-1
	if(rb(i).lt.radmin) go to 30
c-	ird=0 
	i0=i
	rbot=rb(i)
	rtop=rt(i)
c-	if(rb(i).lt.radmin.and.rt(i).ge.radmin) ird=1
c-49	if(ird.eq.2) then
c-	  i0=n+1
c-	  rbot=rb(i)
c-	  rtop=radmin
c-	  ird=0
c-	endif
c-	if(ird.eq.1) then
c-	  rbot=radmin
c-	  ird=2
c-	endif
c		write(6,*)'rbot,rtop=',rbot*6371.,rtop*6371.,'i0=',i0
c*
c       want to divide the radial integration into steps of one-tenth
c       of a wavelength.
        nrad=50
c		write(6,*)'rbot,rtop=',rbot*6371.,rtop*6371.
c		write(6,*)'nrad=',nrad
        dr=-dble(rtop-rbot)/dble(nrad)
c	Prescribe initial values at free surface.
	if(i.eq.n) then
	  b3(1)=1.d0
	  b3(2)=0.d0
	  b3(3)=(bmf1(n,nrad+1)/bmf3(n,nrad+1))*(-1.d0/(fl*(fl+1.d0)))
	  b3(4)=0.d0
		write(6,*)'initial b3 at surface=',b3
	rf=1.0
	  b1w=real(b3(1))
	  b2w=real(b3(2))
	  b3w=real(b3(3))
	  b4w=real(b3(4))
	do 58 k=500,1,-1
	if(rf.lt.rw(k)) go to 58
c		lk=k-50*(k/50)
c
	r0=dble(rw(k))
	if(r0.lt.dble(radmin)) r0=dble(radmin)	
	tau1=mu(i)/eta1(i)
	tau2=mu2(i)/eta(i)
	mus=mu(i)*s0r*(s0r+tau2)/((s0r+tau2)*(s0r+tau1)+mu(i)*s0r/eta(i))
	lams=kappa(i)-2.*mus/3.
	biga=lams+2.*mus   
c
	densi(k)=real(denss)
	ur(k)=b1w
	dur(k)=(b2w-(lams/rf)*(2.*b1w-fl21*b3w))/biga
	vr(k)=b3w
	dvr(k)=(1./rf)*(b3w-b1w)+b4w/mus 
	sur(k)=b2w
58	continue  
c		write(6,*)'n,nrad+1=',n,nrad+1
c		write(6,*)'initial bmf1(n,nrad+1)=',bmf1(n,nrad+1)
c		write(6,*)'initial bmf2(n,nrad+1)=',bmf2(n,nrad+1)
c		write(6,*)'initial bmf3(n,nrad+1)=',bmf3(n,nrad+1)
c		write(6,*)'initial bmf4(n,nrad+1)=',bmf4(n,nrad+1)
c		write(6,*)'initial bmf5(n,nrad+1)=',bmf5(n,nrad+1)
c		write(6,*)'INITIAL b3=',b3
c		pause
	endif
c**
c	Store b3-array in ya1,ya2,ya3,ya4 in order to perform spline
c	interpolation of required quantities accurately.
c	Here we are at top of layer #i and want to store things with
c	index nrad+1 (index 1 is bottom of layer, index nrad+1
c	is top of layer, and we are starting at top of layer).
	ya1(nrad+1)=b3(1)
	ya2(nrad+1)=b3(2)
	ya3(nrad+1)=b3(3)
	ya4(nrad+1)=b3(4)
c**
c	Integrate equations of motion to propagate solutions down from
c	rtop to rbot
c	mus=mu(i)
c	lams=kappa(i)-2.*mu(i)/3.
c	biga=lams+2.*mus 
	r=dble(rtop)
	do 31 j=1,nrad
	rsav=r
c	Now do b3.
	r=rsav
	htot=dr
	eps=3.d-7
	yscal(1)=dabs(b3(1))+fl*dabs(b3(3))+dabs(b3(2))/fl+dabs(b3(4))
	yscal(2)=dabs(b3(2))+fl*dabs(b3(4))+dabs(b3(1))*fl+dabs(b3(3))*fl**2
	yscal(3)=dabs(b3(1))/fl+dabs(b3(3))+dabs(b3(2))/fl**2+dabs(b3(4))/fl
	yscal(4)=dabs(b3(2))/fl+dabs(b3(4))+dabs(b3(1))+dabs(b3(3))*fl
c		write(6,*)'j=',j,'yscal=',yscal
c		write(6,*)'entering derivs4 #1'
	call derivs4(i,r,b3,dydr)
c		write(6,*)'entering bsstep #3F'
c		write(6,*)'b3=',b3
c		write(6,*)'dydr=',dydr
c		write(6,*)'r,htot=',r,htot
c		write(6,*)'eps=',eps
c		write(6,*)'yscal=',yscal
	call bsstep(b3,dydr,4,r,htot,eps,i,yscal,hdid,hnext)
	if((hdid-htot)/htot.ne.0.d0) pause '(hdid-htot)/htot.ne.0 in vsphm'
c		write(6,*)'out of bsstep: r=',r
	bm(1)=bmf1(i0,nrad+1-j)
	bm(2)=bmf2(i0,nrad+1-j)
	bm(3)=bmf3(i0,nrad+1-j)
	bm(4)=bmf4(i0,nrad+1-j)
	bm(5)=bmf5(i0,nrad+1-j)
c		write(6,*)'bm=',bm
	if(bm(1).eq.0.d0.or.bm(2).eq.0.d0) go to 29
	if(bm(3).eq.0.d0.or.bm(4).eq.0.d0) go to 29
c	The minors will have zero values only at the base of the final layer.
	valt=dlog(dabs(bm(4))+(1.d0/fl**2)*dabs(bm(1)))
     &	-dlog(dabs(bm(2))+(1.d0/fl**2)*dabs(bm(4)))-dlog(fl)
c			if(i.eq.2) write(6,*)'before b3-bm relation: b3(1)=',b3(1)
	if(valt.lt.1.5d0) then
c		write(6,*)'use 1st relation'
c	Update values of b3(2) and b3(4) using the relation between
c	the minors and the displacement-stress vector components.
	b3(2)=(bm(4)/bm(2))*b3(1)+(bm(1)/bm(2))*b3(3)
	b3(4)=(1.d0/(fl*(fl+1.d0)))*(bm(1)/bm(2))*b3(1)+
     &	(bm(3)/bm(2))*b3(3)
c		write(6,*)'used 1st relation'
c			if(i.eq.2) write(6,*)'used 1st relation'
	endif
	if(valt.ge.1.5d0) then
c		write(6,*)'use 2nd relation'
c	Update values of b3(1) and b3(3) using the relation between
c	the minors and the displacement-stress vector components.
c			if(i.eq.2) write(6,*)'bm(3)=',bm(3)
c			if(i.eq.2) write(6,*)'bm(5)=',bm(5)
c			if(i.eq.2) write(6,*)'b3(2)=',b3(2)
c			if(i.eq.2) write(6,*)'bm(1)=',bm(1)
c			if(i.eq.2) write(6,*)'bm(5)=',bm(5)
c			if(i.eq.2) write(6,*)'b3(4)=',b3(4)
	b3(1)=(bm(3)/bm(5))*b3(2)-(bm(1)/bm(5))*b3(4)
	b3(3)=-(1.d0/(fl*(fl+1.d0)))*(bm(1)/bm(5))*b3(2)+
     &	(bm(4)/bm(5))*b3(4)
c		write(6,*)'used 2nd relation'
c			if(i.eq.2) write(6,*)'used 2nd relation'
	endif
c			write(6,*)'after b3-bm relation: b3(1)=',b3(1)
29	continue
c**
c	Store b3-array in ya1,ya2,ya3,ya4 in order to perform spline
c	interpolation of required quantities accurately.
	ya1(nrad+1-j)=b3(1)
	ya2(nrad+1-j)=b3(2)
	ya3(nrad+1-j)=b3(3)
	ya4(nrad+1-j)=b3(4)
c**
	r=rsav+dr
31	continue	
c		write(6,*)'out of do 31'
c**
c	Update values of x1- and x2- integrals and eigenfunctions in layer #i.
c		write(6,*)'entering update'
	call update(radmin,rtop,rbot,nrad,x1)
c		write(6,*)'out of update'
c**
c		write(6,*)'BOTTOM LAYER ',i,'x1=',x1,'x2=',x2,'dr=',dr 
	wo=wo+x1*real(-dr)
c		write(6,*)'BOTTOM LAYER ',i,'NEW wo=',wo
c		write(6,*)'BOTTOM LAYER ',i,'delta-wo=',x1*real(-dr)
c		write(6,*)'b3=',b3
c-	if(ird.eq.2) go to 49
c**
38	bmag=abs(real(b3(1)))+real(l)*abs(real(b3(3)))
	if(bmag.lt.1.e+5) go to 30
		write(6,*)'re-normalizing'
	do 37 j1=1,4
37	b3(j1)=b3(j1)*1.d-10
	do 34 k=1,500
	ur(k)=ur(k)*1.e-10
	dur(k)=dur(k)*1.e-10
	sur(k)=sur(k)*1.e-10
	vr(k)=vr(k)*1.e-10
34	dvr(k)=dvr(k)*1.e-10
	wo=wo*1.e-20
	go to 38
30	continue
c------------
152	continue
c	Normalize the eigenfunctions for output.
	fac=(2.464e-2)*(6371./bigr)**2
	fac=sqrt(fac)/sqrt(wo)
c	fac=normalizing constant involving the units of mu,eta,and earth radius.
c------
c	Get rid of poorly-integrated modes
	if(bmftop.gt.1.e-4) fac=0.
c------
	sj=real(sdec(nzero,l))
	do k=1,ndep
	ur(k)=(ur(k)/rw(k))*fac/sj
	dur(k)=dur(k)*fac/sj
	vr(k)=(vr(k)/rw(k))*fac/sj
	dvr(k)=dvr(k)*fac/sj
	write(4) l,ur(k),dur(k),vr(k),dvr(k)
c		write(6,*)'reigen20','r=',rw(k)*6371,l,dur(k),ur(k),dvr(k),vr(k)
	enddo
c		write(6,*)'BEFORE NORMALIZATION: ur(40+ndep),dur(40+ndep)=',ur(40+ndep),dur(40+ndep)
c		write(6,*)'BEFORE NORMALIZATION: ur(42+ndep),dur(42+ndep)=',ur(42+ndep),dur(42+ndep)
	do k=ndep+1,ndep+40
	ur(k)=(ur(k)/rw(k))*fac/sj
	dur(k)=dur(k)*fac/sj
	vr(k)=(vr(k)/rw(k))*fac/sj
	dvr(k)=dvr(k)*fac/sj
	write(8) l,dur(k),ur(k),dvr(k),vr(k)
c		write(6,*)'VSPH10: rw(',k,')=',rw(k)*bigr
c		write(6,*)'k=',k,'r=',rw(k)*6371,l,dur(k),ur(k),dvr(k),vr(k)
	enddo
	ur(41+ndep)=(ur(41+ndep)/rw(41+ndep))*fac
	dur(41+ndep)=dur(41+ndep)*fac
	vr(41+ndep)=(vr(41+ndep)/rw(41+ndep))*fac
	dvr(41+ndep)=dvr(41+ndep)*fac
	ur(42+ndep)=(ur(42+ndep)/rw(42+ndep))*fac
	dur(42+ndep)=dur(42+ndep)*fac
	vr(42+ndep)=(vr(42+ndep)/rw(42+ndep))*fac
	dvr(42+ndep)=dvr(42+ndep)*fac
	write(8) int(l),ur(42+ndep),vr(42+ndep),dur(42+ndep),dvr(42+ndep)
	write(8) int(l),ur(41+ndep),vr(41+ndep),dur(41+ndep),dvr(41+ndep)
c	write(6,*)'k=',42+ndep,int(l),dur(42+ndep),ur(42+ndep),dvr(42+ndep),vr(42+ndep)
c	write(6,*)'reigen10',int(l),ur(42+ndep),vr(42+ndep),dur(42+ndep),
c     &	dvr(42+ndep)
c	write(6,*)'reigen10',int(l),ur(41+ndep),vr(41+ndep),dur(41+ndep),
c     &	dvr(41+ndep)
c------------------
48	continue 
50	continue
141	format(i5,e13.6e2)
	end 

