	program vsphdep
c*
c	Calculate Greens functions at surface and one other depth
c	specified in input. (differs from VSPH which does Greens
c	functions at surface and one other pre-specified depth within program).
c*
c	Read in singular values s=-s(j) of N-layer earth model
c	for spheroidal modes of degree l=lmin-lmax, read from
c	input file 'decay4.out'.
c	Expect maximum
c	number of [maxmod] decay times. 
c	Find coefficients which control
c	geomeetric dependence of displacement field.  Write
c	out these coefficients in the same order as the decay times were
c	read in.  Residues at s=-s(j) are evaluated analytically.
c	The displacement coefficients are written for
c	the vertical and longitudinal components in that order.  
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
	parameter (maxmod=19)
	character*80 aread
	real*4 kappa,mu,mu2,dens(200),mupr
cN	
	dimension depf(41)
c-
	real*8 sdecr,cof1r,cof2r,cof1l,cof2l 
	real*8 cof1(maxmod,6500),cof2(maxmod,6500)
	real*8 aj(4,4),d1(4),d2(4),s0
	real*8 sdec(maxmod,6500)
	real*8 x1,dmds,dlds,third
	real*8 r,dr,r2,biga,fl21,flp1,sfac
	real*8 bigde,pdds  
	real*8 tau1,tau2
cN
	real*8 us(41),ups(41),vs(41),vps(41)
c-
c	real*8 us1,us2,us3,us4,us5,us6,us7,us8,us9,us10,us11
c	real*8 ups1,ups2,ups3,ups4,ups5,ups6,ups7,ups8,ups9,ups10,ups11
c	real*8 vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8,vs9,vs10,vs11
c	real*8 vps1,vps2,vps3,vps4,vps5,vps6,vps7,vps8,vps9,vps10,vps11
	real*8 wo,ur,vr,dur,dvr,c1cof,c2cof,c3cof,c4cof,dc11,dc12
	real*8 up,vp,y7,r0
	real*8 lams,mus
	real*8 c(4),b1(4),b2(4)
cN
	dimension y1(maxmod,40),y2(maxmod,40),y3(maxmod,40),y4(maxmod,40)
c-
	dimension urf(maxmod,2),vrf(maxmod,2),durf(maxmod,2),dvrf(maxmod,2)
	dimension mdec(6500)
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),
     &	eta(200),eta1(200),l,s0,r0,mupr(200)
	pi=3.14159265
	twopi=2.*pi
	third=1.d0/3.d0
c	bigr=6371.
c	bigr=earth's radius in km.
c	The displacement coefficients are evaluated at depths
c	(1,4,7,...,28)* (depfac) km.
c	Read in earth model
	open(4,file='earth.model')
	rewind(4)
	read(4,5) n,nd,bigr 
	write(6,5) n,nd,bigr  
5	format(2i2,f10.3)
c ** *	
	  etala=1.e+13
	do 10 j=1,n
	i=2*j-1
c
	read(4,14) aread 
	read(aread,18,err=37) rb(i),rt(i),dens(i),kappa(i),mu(i),
     &	mupr(i),eta(i),eta1(i)
	  etaval=min(eta(i),eta1(i))
	if(eta1(i).eq.0.0) go to 37
	go to 12
37	eta1(i)=1.e+13
	read(aread,17,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),
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
	write(6,18) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),
     &	eta(i),eta1(i)
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
14	format(a80)
15	format(5f9.3,e13.6e2)
17	format(6f9.3,e13.6e2)
18	format(6f9.3,2e13.6e2)
	n=2*n 
c--
	  helast=6371.-helast
	  write(6,*)'elastic plate thickness=',helast
c--
c ** *
	close(4)
	  write(6,*)'reading decay times'
	open(4,file='decay4.out')
	rewind(4)
	read(4,21) ndd,bigrr,depfac
	if(ndd.ne.nd) write(6,*)'read incorrect #layers in decay4.out'
	if(ndd.ne.nd) stop 
c	Since the characteristic relaxation time is equal
c	to eta/(mu + mub) rather than eta/mu, rescale eta accordingly.
c	Note that mub=mupr*mu/(mu-mupr).
c----	do 18 i=1,n
c----18	eta(i)=eta(i)*(1.-mupr(i)/mu(i))
	idec=0
	is=0 
16	idec=idec+1 
c	if(l.lt.0) lgra=llast   
	read(4,44,end=20) l,sdecr,cof1r,cof2r
	if(is.eq.0) llast=l
	lmax=l
	if(l.eq.llast) mdec(l)=idec
	if(l.eq.llast) sdec(idec,l)=sdecr
	if(is.eq.0) lmin=l  
	if(is.eq.0) is=1
	if(l.eq.llast) cof1(idec,l)=cof1r 
	if(l.eq.llast) cof2(idec,l)=cof2r 
	if(l.eq.llast) go to 16 
	sdec(1,l)=sdecr
	cof1(1,l)=cof1r
	cof2(1,l)=cof2r
	idec=1
	mdec(l)=1
	llast=l 
	go to 16
20	close(4)
	write(6,*)'read decay times: lmin,lmax,lgra=',lmin,lmax,lgra  
	do 19 l=lmin,lmax
19	mdec(l)=mdec(l)+1
	do idep=1,40
cN
	de=28.+(0.001 - 28.)*real(idep-1)/39.
	depf(idep)=1.-de*depfac/bigr
	enddo
c-
c	dep1=1.-28.*depfac/bigr
c	dep2=1.-25.*depfac/bigr
c	dep3=1.-22.*depfac/bigr
c	dep4=1.-19.*depfac/bigr
c	dep5=1.-16.*depfac/bigr
c	dep6=1.-13.*depfac/bigr
c	dep7=1.-10.*depfac/bigr
c	dep8=1.-7.*depfac/bigr
c	dep9=1.-4.*depfac/bigr
c	dep10=1.-1.*depfac/bigr
c*
	write(6,*)'depth level of extra output Greens function (km)?'
	read(5,*) edep
cN
	depf(41)=1.-edep/bigr
c-
c	dep11=1.-edep/bigr
c*
c	open(8,file='vsph.u')
	open(2,file='vsph.out',form='unformatted')
	write(2) nd
	write(6,21) nd 
21	format(i2,3f10.3)
	do 50 l=lmin,lmax
	write(6,*)'l=',l
	depmax=1.-helast/bigr-3.0*(2.*pi/(real(l)+0.5))
	depmin=1.-helast/bigr+1.5*(2.*pi/(real(l)+0.5))
c	depmax=dimensionless radius 3.0 wavelengths into the earth.
	do 48 idec=1,mdec(l)-1
	s0=-sdec(idec,l)
	cof1l=cof1(idec,l)
	cof2l=cof2(idec,l)
	wo=0.d0
	fl21=dble(l*(l+1))
	fl12=dble((l-1)*(l+2))
	flp1=dble(l+1)
	om2=dble(w0*w0)
c	Start matrix propagation.
c	Give values at CMB those appropriate for homogeneous solid
c	in 0<r<rb(1).
	is=0
	do 30 i=1,n
	if(i.lt.2) go to 30
c	Start integration at most 2 wavelengths into the earth.
	if(rb(i).lt.depmax) go to 30
	if(rt(i).gt.depmin) go to 30
	  sfac=s0+dble(mu(i)/eta(i))
	tau1=dble(mu(i)/eta1(i))
	tau2=dble(mu2(i)/eta(i))
	mus=mu(i)*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
	lams=dble(kappa(i)-2.*mus/3.)
	bigde=(s0+tau1)*(s0+tau2)+mu(i)*s0/eta(i)
	pdds=(s0+tau1)+(s0+tau2)+mu(i)/eta(i)
c	dmds=dble((mu(i)-mupr(i))*(mu(i)/eta(i)))/(sfac*sfac)
	dmds=real(mus/s0 + mus/(s0+tau2) -mus*pdds/bigde)
c	dmds is the derivative of mus w.r.t. s.
	biga=lams+2.d0*mus
	r0=dble(rb(i))
c	Give boundary conditions at top of layer (i-1) appropriate for
c	homogeneous solid in 0 < r < rt(i-1).
	if(is.ne.0) go to 28
	is=1
	j=i-1
	call matra1(aj,j,1)
	do 27 k=1,4
	b1(k)=aj(k,1)
	b2(k)=aj(k,2)
27	continue
28	continue
c	Determine coefficient matrices d1 and d2 for layer i.
	do 33 m=1,2
	if(m.eq.1) call equal(b1,c,4)
	if(m.eq.2) call equal(b2,c,4)
	call matra1(aj,i,0)
	if(m.eq.1) call ainver(d1,aj,c,4)
	if(m.eq.2) call ainver(d2,aj,c,4)
33	continue
	nrad=50
29	dr=dble(rt(i)-rb(i))/dble(nrad)
c	Integrate equations of motion to propagate solutions up to
c	the surface.
	r=dble(rb(i))
	x1=0.d0
	do 31 j=1,nrad
	r2=r*r 
	r=r+dr
	r0=dble(r)  
	call matra1(aj,i,2)
c		write(6,*)'aj=',aj
	do 32 m=1,2
	if(m.eq.1) call prodr(aj,d1,4,c)
	if(m.eq.2) call prodr(aj,d2,4,c)
	if(m.eq.1) call equal(c,b1,4)
	if(m.eq.2) call equal(c,b2,4)
c		write(6,*)'m=',m
c		write(6,*)'b1=',b1
32	continue
	c1cof=b1(1)*cof1l+b2(1)*cof2l
	c2cof=b1(2)*cof1l+b2(2)*cof2l
	c3cof=b1(3)*cof1l+b2(3)*cof2l
	c4cof=b1(4)*cof1l+b2(4)*cof2l
	rf=real(r)
cN
	do idep=1,41
	if(rf.lt.depf(idep)) then
	us(idep)=c1cof/r 
	vs(idep)=c3cof/r  
	endif
	enddo
c-
	dc11=(b1(2)-(lams/r)*(2.d0*b1(1)-fl21*b1(3)))/biga 
	dc12=(b2(2)-(lams/r)*(2.d0*b2(1)-fl21*b2(3)))/biga 
	up=dc11*cof1l+dc12*cof2l
c	up is equal to dy1/dr.
	dc11=(1.d0/r)*(b1(3)-b1(1))+b1(4)/mus 
	dc12=(1.d0/r)*(b2(3)-b2(1))+b2(4)/mus 
	vp=dc11*cof1l+dc12*cof2l
c	vp is equal to dy3/dr.
cN
	do idep=1,41
	if(rf.lt.depf(idep)) then
	ups(idep)=up
	vps(idep)=vp
	endif
	enddo
c-
	y7=2.d0*c1cof-fl21*c3cof
	bigf=(2.*c1cof-fl21*c3cof)/r
	x1=x1+(third*(2.d0*up-bigf)**2+fl21*((r*vp-c3cof+c1cof)**2+
     &	fl12*c3cof**2)/r2)*r2*dmds
c	x1=x1+((r*c2cof+2.d0*mus*y7)**2/biga**2)*(dlds+2.d0*dmds)
c	x1=x1+(fl21*(r*c4cof)**2/mus**2-4.d0*up*y7+2.d0*(fl21*c3cof
c &	*(c1cof-c3cof)-c1cof*y7))*dmds 
31	continue
c	  if(idec.eq.1) write(6,*)'top layer #',i ,'r=',r 
c	  if(idec.eq.1) write(6,*)'c1cof,c2cof,c3cof,c4cof=',c1cof,c2cof,c3cof,c4cof
c	  if(i.gt.51) pause 
c	if(i.eq.54)  write(6,*)'b1=',b1
c	if(i.eq.54)  write(6,*)'b2=',b2
c	  if(idec.eq.1) write(6,*)'top layer #',i ,'x1=',x1,'dr=',dr,'IDEC=1'
c	  if(idec.eq.1) write(6,*)'top layer #',i ,'b1=',b1,'b2=',b2,'IDEC=1'
	wo=wo+x1*dr
c	  if(idec.eq.1) write(6,*)'top layer #',i ,'wo=',wo,'IDEC=1'
c	Make sure that matrix propagation  is done accurately.
	do 53 m=1,2
	if(m.eq.1) call equal(b1,c,4)
	if(m.eq.2) call equal(b2,c,4)
	call matra1(aj,i,1)
	if(m.eq.1) call prodr(aj,d1,4,c)
	if(m.eq.2) call prodr(aj,d2,4,c)
	if(m.eq.1) call equal(c,b1,4)
	if(m.eq.2) call equal(c,b2,4)
53	continue
46	bmag=abs(real(b1(1)))+real(l)*abs(real(b1(3)))
	if(bmag.lt.1.e+10) go to 30
	do 36 j1=1,4
	b1(j1)=b1(j1)*1.d-10
36	b2(j1)=b2(j1)*1.d-10
	wo=wo*1.d-20
cN
	do idep=1,41
	us(idep)=us(idep)*1.d-10
	vs(idep)=vs(idep)*1.d-10
	ups(idep)=ups(idep)*1.d-10
	vps(idep)=vps(idep)*1.d-10
	enddo
c-
30	continue
	  if(idec.eq.1) write(6,*)'TOP--IDEC=1. b1=',b1,'b2=',b2,'wo=',wo
	ur=b1(1)*cof1l+b2(1)*cof2l
	vr=b1(3)*cof1l+b2(3)*cof2l
	dc11=(b1(2)-(lams/1.d0)*(2.d0*b1(1)-fl21*b1(3)))/biga 
	dc12=(b2(2)-(lams/1.d0)*(2.d0*b2(1)-fl21*b2(3)))/biga 
	dur=dc11*cof1l+dc12*cof2l
c	dur is equal to dy1/dr at surface.
	dc11=(1.d0/1.d0)*(b1(3)-b1(1))+b1(4)/mus 
	dc12=(1.d0/1.d0)*(b2(3)-b2(1))+b2(4)/mus 
	dvr=dc11*cof1l+dc12*cof2l
c	dvr is equal to dy3/dr at surface.
	urf(idec,1)=real(ur)
	vrf(idec,1)=real(vr)
	durf(idec,1)=real(dur)
	dvrf(idec,1)=real(dvr)
cN
	ur=us(41)
	vr=vs(41)
	dur=ups(41)
	dvr=vps(41)
c-
c	NOTE: evaluate displacement coeff. at 0 and [edep] km depth.
	urf(idec,2)=real(ur)
	vrf(idec,2)=real(vr)
	durf(idec,2)=real(dur)
	dvrf(idec,2)=real(dvr)
	c1cof=b1(1)*cof1(idec,l)+b2(1)*cof2(idec,l)
	c2cof=b1(2)*cof1(idec,l)+b2(2)*cof2(idec,l)
	c3cof=b1(3)*cof1(idec,l)+b2(3)*cof2(idec,l)
	c4cof=b1(4)*cof1(idec,l)+b2(4)*cof2(idec,l)
c*
c	Test if an unstable mode is being considered.  If so, set
c	the eigenfunctions evaluated at the observation depth to zero.
	stre=dvrf(idec,1)-vrf(idec,1)+urf(idec,1)
cN
	umax=0.
	do idep=1,40
	if(dabs(us(idep)).gt.dble(umax)) umax=real(us(idep))
	enddo
	rat=abs(stre/umax)
c-
	if(rat.gt.5.e-1) then
	write(6,*)'UNSTABLE MODE: l,idec=',l,idec
	urf(idec,1)=0.
	vrf(idec,1)=0.
	durf(idec,1)=0.
	dvrf(idec,1)=0.
	urf(idec,2)=0.
	vrf(idec,2)=0.
	durf(idec,2)=0.
	dvrf(idec,2)=0.
	endif
c*	
	  write(6,*)'surface: U, R, V, S =', c1cof,c2cof,c3cof,c4cof
c	  write(6,*)'mus=',mus  
c	Determine residues of displacement functions.
c	const=sqrt((real(l)+0.5)/twopi)
c	facl=sqrt(real(l*(l+1)))
c	facl1=facl*sqrt(real((l-1)*(l+2)))
c	Determine Green's functions.
	kr=0
47	kr=kr+1
	if(kr.eq.41) go to 48 
	usf=real(us(kr))
	vsf=real(vs(kr))
	upsf=real(ups(kr))
	vpsf=real(vps(kr))
	 y1(idec,kr)=upsf/real(wo) 
	 y2(idec,kr)=(usf)/real(wo)
	 y3(idec,kr)=vpsf/real(wo)
	 y4(idec,kr)=(vsf)/real(wo) 
	md1=mdec(l)
41	format(i4,e13.6e2)
44	format(i4,3d22.15)
	go to 47
48	continue 
c	Solve for d0,a1,a3,b1,b2 for each decay time,
	nm1=md1-1
	fac=(2.464e-2)*(6371./bigr)**2
	write(6,*)'fac=',fac,'l,nm1=',l,nm1
c	fac=normalizing constant involving the units of mu,eta,and earth radius.
	do 65 i=1,nm1
c	Normalize geometrical factors w.r.t. respective decay time
c	and produce output units of cm displacement.  Output is in the form
c	l,d0(j)/s(j),a1(j)/s(j),a2(j)/s(j),b1(j)/s(j),b2(j)/s(j)
	sj=real(sdec(i,l))
cN	k from 1 to 40
	do 64 k=1,40
	y1(i,k)=y1(i,k)*fac/sj 
	y2(i,k)=y2(i,k)*fac/sj 
	y3(i,k)=y3(i,k)*fac/sj 
	y4(i,k)=y4(i,k)*fac/sj 
	write(2) l,y1(i,k),y2(i,k),y3(i,k),y4(i,k)
	dept=depfac*3.*real(11-k)
c	write(8,*) l,i,dept,y1(i,k),y2(i,k),y3(i,k),y4(i,k)
	  y1x=y1(i,k)*urf(i,1)*sj*1.e+6
	  if(k.eq.40) write(6,*) l,y1x
c	write(6,*) 'depth',bigr-bigr*depf(k),l,y1(i,k),y2(i,k),y3(i,k),y4(i,k)
43	format(i4,4e13.6e2)
64	continue 
	write(2) l,urf(i,1),vrf(i,1),durf(i,1),dvrf(i,1)
	write(2) l,urf(i,2),vrf(i,2),durf(i,2),dvrf(i,2)
c	write(8,*) "surface:"
c	write(8,*) l,i,dept,durf(i,1),urf(i,1),dvrf(i,1),vrf(i,1)
65	continue
50	continue
	end 
	subroutine matra1(aj,i,n)
c	Compute matrix elements for spheroidal modes (decay times).
c	n=0: evaluate at bottom of layer i
c	n=1: evaluate at top of layer i
c	Oherwise evaluate at r0.
	real*8 r0
	real*4 mu,mu2,kappa,mupr
	real*8 mus,lams,sfac,biga,r,flp,flm,r2,flp3,flm3,s0
	real*8 tau1,tau2
	real*8 a1p,a1m,a2p,a2m,a3p,a3m,a4p,a4m
	real*8 b1p,b1m,b2p,b2m,b3p,b3m,b4p,b4m,rlp,rlm  
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),
     &	eta(200),eta1(200),l,s0,r0,mupr(200)
	real*8 aj(4,4)
	r=r0
	tau1=dble(mu(i)/eta1(i))
	tau2=dble(mu2(i)/eta(i))
	mus=mu(i)*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
	lams=dble(kappa(i)-2.*mus/3.)
c		write(6,*)'mus,lams=',mus,lams
	biga=lams+2.d0*mus   
10	if(n.eq.1) r=dble(rt(i))
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
	aj(1,4)=mus*(flm+1.d0)*a1m+biga*b1m
	aj(2,2)=mus*(flp+1.d0)*a2p+biga*b2p
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
	return
	end 
	subroutine prodr(a,c,n,b)
c	Form matrix product b=(A)c
	real*8 b(n) 
	real*8 a(n,n),c(n)
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
	real*8 a(4),b(4,4),c(4),bsave,fac
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
	subroutine equal(a,b,n)
c	Set the n x 1 matrix (b) equal to (a).
	real*8 a(n),b(n)
	do 10 j=1,n
10	b(j)=a(j)
	return
	end 

