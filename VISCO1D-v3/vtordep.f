	program vtordep  
c*
c	Calculate Greens functions at surface and one other depth
c	specified in input. (differs from VTOR which does Greens
c	functions at surface and one other pre-specified depth within program).
c*
c	Read in singular values s=-s(j) of N-layer earth model
c	for toroidal modes.
c	Expect maximum
c	number of [maxmod] decay times.
c	Find coefficients which control geometric dependence
c	of displacement field. 
c	Write out these coefficients in the
c	same order as the decay times were read in.  Residues at
c	s=-s(j) are evaluated analytically.  
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
	parameter (maxmod=5)
	character*80 aread
	real*4 dens(200),kappa(200),mu,mu2,mupr 
cN	
	dimension depf(41)
	dimension wts(41),wtps(41)
c-
	real*8 mus,sfac,bigde,pdds  
	real*8 tau1,tau2
	real*8 sdecr,sdec1(maxmod,6500)
	dimension mdec(6500)
	real*8 aj(2,2),b(2),d(2),s0
	dimension y1(maxmod,40),y2(maxmod,40),wr(maxmod,2),
     &	dwr(maxmod,2)
	common/mat/rb(200),rt(200),mu(200),mu2(200),eta(200),eta1(200),
     &	l,s0,r,mupr(200)
	pi=3.14159265
	twopi=2.*3.14159265  
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
	nm1=nd-1
c ** *
	  etala=1.e+13
	do 10 j=1,n
	i=2*j-1
c
	read(4,14) aread 
	read(aread,18,err=27) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i),eta1(i)
	  etaval=min(eta(i),eta1(i))
	if(eta1(i).eq.0.0) go to 27
	go to 12
27	eta1(i)=1.e+13
	read(aread,17,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i)
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
	open(4,file='decay.out')
	rewind(4)
	read(4,21) ndd,bigrr,depfac
	if(ndd.ne.nd) write(6,*)'read incorrect #layers in decay.out'
	if(ndd.ne.nd) stop 
	write(6,21) nd,bigrr,depfac
c	Since the characteristic relaxation time is equal
c	to eta/(mu + mub) rather than eta/mu, rescale eta accordingly.
c	Note that mub=mupr*mu/(mu-mupr).
c----	do 18 i=1,n
c----18	eta(i)=eta(i)*(1.-mupr(i)/mu(i))
	lmin=6600  
	do 19 l=1,6500
	mdec(l)=0
19	continue
16	read(4,41,end=20) l,sdecr
	mdec(l)=mdec(l)+1
	kdec=mdec(l)
	sdec1(kdec,l)=sdecr
	write(6,41) l,sdecr
	lmax=l 
	if(l.lt.lmin) lmin=l 
	go to 16
20	close(4)
	  write(6,*)'lmax=',lmax
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
	open(2,file='vtor.out',form='unformatted')
	write(2) nd
21	format(i2,3f10.3)
	do 50 l=lmin,lmax 
	write(6,*)'l=',l
	depmax=1.-helast/bigr-3.0*(2.*pi/(real(l)+0.5))
	depmin=1.-helast/bigr+1.5*(2.*pi/(real(l)+0.5))
c	depmax=dimensionless radius 4.0 wavelengths into the earth.
	  nm1=mdec(l)
	do 48 idec=1,nm1
	s0=-sdec1(idec,l)
cc	Do perfectly elastic case last.
c	Start matrix propagation.
c	b(1)=1.d0
c	b(2)=0.d0
	is=0
	wo=0.
	do 30 i=1,n
	if(i.lt.2) go to 30
c	Start integration at most 4 wavelengths into the earth.
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
	r=rb(i)
	dr=(rt(i)-rb(i))/50.
	  sfac=s0+dble(mu(i)/eta(i))
	tau1=dble(mu(i)/eta1(i))
	tau2=dble(mu2(i)/eta(i))
	mus=mu(i)*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
	bigde=(s0+tau1)*(s0+tau2)+mu(i)*s0/eta(i)
	pdds=(s0+tau1)+(s0+tau2)+mu(i)/eta(i)
	x1=0.0
	call matra(aj,i,0)
	call ainver(d,aj,b,2)
	do 31 j=1,50
	r=r+dr 
	facl=real((l-1)*(l+2))
	call matra(aj,i,2)
	call prod(aj,d,2,b) 
cN
	do idep=1,41
	if(r.lt.depf(idep)) then
	wts(idep)=real(b(1))/r
	wtps(idep)=real(b(1))/r+real(b(2)/mus)
	endif
	enddo
c-
	x1=x1+(r*r*real(b(2)*b(2))/real(mus*mus)
     &	+facl*real(b(1)*b(1)))
31	continue
c	dmds=(mu(i)-mupr(i))*(mu(i)/eta(i))/real(sfac*sfac)
	dmds=real(mus/s0 + mus/(s0+tau2) -mus*pdds/bigde)
	wo=wo+x1*dr*dmds
c	  write(6,*)'top layer',i,'wo=',wo,'ws1',ws1
c	Note: wo is equal to (d/ds) I2 (Notes, p.11) needed to compute
c	residue at s=-s sub j.
	call matra(aj,i,1)
	call prod(aj,d,2,b)
	  bmag=abs(real(b(1))/r)+abs(real(b(2)/mus))
	if(bmag.lt.1.e+10) go to 30
	b(1)=b(1)*1.d-10
	b(2)=b(2)*1.d-10
	wo=wo*1.e-20
cN
	do idep=1,41
	wts(idep)=wts(idep)*1.e-10
	wtps(idep)=wtps(idep)*1.e-10
	enddo
c-
30	continue
	wr(idec,1)=real(b(1))
	dwr(idec,1)=real(b(1))/1.d0+real(b(2)/mus)
c	dwr=d(wr)/dr at surface.
cN
	wr(idec,2)=wts(41)
	dwr(idec,2)=wtps(41)
c-
c	NOTE: evaluate displ. coeff. at depths 0 and [edep] km. 
c	Test if an unstable mode is being considered.  If so, set
c	the eigenfunctions evaluated at the observation depth to zero.
	stre=dwr(idec,1)-wr(idec,1)
cN
        wmax=0.
        do idep=1,40
        if(abs(wts(idep)).gt.wmax) wmax=wts(idep)
        enddo
        rat=abs(stre/wmax)
	write(6,*)'stre,wmax=',stre,wmax
c-
	if(rat.gt.1.e-1) then
	write(6,*)'UNSTABLE MODE: l,idec=',l,idec
	wr(idec,1)=0.
	dwr(idec,1)=0.
	wr(idec,2)=0.
	dwr(idec,2)=0.
	endif
	  write(6,*)'at surface W, T=',b 
	wo=wo*real(l*(l+1))
43	format(i4,2e13.6e2)
c	Determine matrix elements for curve-fit in s-domain.
c	const=sqrt((real(l)+0.5)/twopi)
c	facl=sqrt(real(l*(l+1)))
c	facl1=facl*sqrt(real((l-1)*(l+2)))
	  write(6,*)'l,idec=',l,idec
c	  pause 
		write(6,*)'wo=',wo
	kr=0
47	kr=kr+1
	if(kr.eq.41) go to 48
	ws=wts(kr)
	wps=wtps(kr)
c	  write(6,*)'ws=',ws
c	  write(6,*)'wps=',wps
c	  write(6,*)'wr=',wr
c	  write(6,*)'wo=',wo
	y1(idec,kr)=(wps-ws)/wo
	y2(idec,kr)=ws/wo
	go to 47
41	format(i4,d17.10)
48	continue
	fac=(2.464e-2)*(6371./bigr)**2
c	fac=normalizing constant involving the units of mu,eta,and earth radius.
	  write(6,*)'nm1=',nm1
	do 65 i=1,nm1
c*	Normalize geometrical factors w.r.t. respective decay times
c	and produce units of cm displacement.  Output is in the
c	form: l,a1(j)/s(j),a2(j)/s(j),b1(j)/s(j),b2(j)/s(j) |j=1,..,N-1
	sj=real(sdec1(i,l))
c		write(6,*)'mode #',i,' sj=',sj
cN
	do 64 k=1,40
	y1(i,k)=y1(i,k)*fac/sj
	y2(i,k)=y2(i,k)*fac/sj
c	write(2,43) l,y1(i,k),y2(i,k)
	  y2x=y2(i,k)*wr(i,1)*sj*1.e+6
	  if(k.eq.7) write(6,*) l,y2x
c	write(6,*) 'depth',bigr-bigr*depf(k),l,y1(i,k),y2(i,k)
	write(2) l,y1(i,k),y2(i,k)
64	continue
	write(2) l,wr(i,1),dwr(i,1),wr(i,2),dwr(i,2)
c	write(6,*) l,wr(i,1),dwr(i,1),wr(i,2),dwr(i,2)
65	continue 
50	continue
	end 
	subroutine matra(aj,i,n)
c	Compute matrix elements for toroidal modes (decay times).
c	n=0: evaluate at bottom of layer i
c	n=1: evaluate at top of layer i
c	Oherwise evaluate at r.
	real*4 mu,mu2,mupr
	real*8 mus,s0,tau1,tau2
	common/mat/rb(200),rt(200),mu(200),mu2(200),eta(200),eta1(200),
     &	l,s0,r,mupr(200)
	real*8 aj(2,2)
	tau1=dble(mu(i)/eta1(i))
	tau2=dble(mu2(i)/eta(i))
	mus=mu(i)*s0*(s0+tau2)/((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
c		write(6,*)'mus( layer',i,')=',mus
c		write(6,*)'mus(',i,')=',mus
c		write(6,*)'eta1,eta=',eta1(i),eta(i)
c		write(6,*)'tau1,tau2=',tau1,tau2
c		write(6,*)'denominator=',((s0+tau2)*(s0+tau1)+mu(i)*s0/eta(i))
c		write(6,*)'------------'
	if(n.eq.1) r=rt(i)
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


