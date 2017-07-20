	subroutine getmin(n,radmin)
c	Using the Earth model and l-value, frequency determined by RWAVESM,
c	evaluate the minors at every radius level from rb(2) to
c	Earth's surface with roughly 10 depth samples per wavelength.
c	---------------------------
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
	real*4 kappa,mu,mu2,mupr
	real*8 htot,eps,rsav,hdid,hnext,yscal(5),dydr(5)
	real*8 r,dr
	real*8 aj4(4,4)
	real*8 b14(4),b24(4),bm(5)
	real*8 s0,r0
		real*8 fl
c	real*4 dens1(200),kapp1(200),mu1(200)
c*
c		dimension rw(500)
	real*8 a1b,a2b,a3b,a4b,c5
	real*8 bmf1,bmf2,bmf3,bmf4,bmf5,bmftop
c*
	common/botcof/a1b(2),a2b(2),a3b(2),a4b(2),c5
c---
cOLD	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),l,w0,r0
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0,r0
c---
	common/jfac0/j01,j02,jm1
	common/jdif/jdif1,jdif3,jdif4,jdif6
	common/jfac/jfac1,jfac2,jfac3,jfac4
c*
	common/minfn/bmf1(201,51),bmf2(201,51),bmf3(201,51),
     &	bmf4(201,51),bmf5(201,51),bmftop
c*

	pi=3.14159265
	bigr=6371.
	fl=dble(l)
c ** *
c		write(6,*)'GETMIN -- eta1(54)=',eta1(54)
c		write(6,*)'GETMIN: entering do 121'
	is=0
	do 121 i=1,n
c		write(6,*)'GETMIN: i=',i
	if(rb(i).lt.radmin) go to 121
c	Use numerical solution above radmin.
	rbot=rb(i)
	rtop=rt(i)
	if(is.ne.0) go to 28
	is=1
	call matra(aj4,i,0)
c		write(6,*)'GETMIN: OUT OF MATRA: aj4=',aj4
	do k=1,4
	b14(k)=aj4(k,1)
	b24(k)=aj4(k,2)
	enddo
c	form second order minors corresponding to chosen b's above.
	bm(1)=b14(1)*b24(2)-b14(2)*b24(1)
	bm(2)=b14(1)*b24(3)-b14(3)*b24(1)
	bm(3)=b14(1)*b24(4)-b14(4)*b24(1)
	bm(4)=b14(2)*b24(3)-b14(3)*b24(2)
	bm(5)=b14(2)*b24(4)-b14(4)*b24(2)
c		write(6,*)'GETMIN: initial b14=',b14,'b24=',b24
28	continue
c**
c*
	nrad=50
c		write(6,*)'layer #',i,'thickness=',(rt(i)-rbot)*bigr,' km'
c		write(6,*)'rt(i)=',rt(i),'rbot=',rbot,'bigr=',bigr
c		write(6,*)'wavnum=',wavnum
c		write(6,*)'        wavelength=',wavlen,' km.  nrad=',nrad
c		pause
	dr=dble(rt(i)-rbot)/dble(nrad)
c		write(8,*)'NUMERICAL: rbot,rtop=',rbot,rt(i)
c	Integrate equations of motion to propagate solutions up to
c	the surface.
c		if(i.eq.2) write(6,*)'bottom layer',i,'bm(5)=',bm(5)
c36	continue
	r=dble(rbot)
c	Do bm using bsstep.
c		write(6,*)'nrad=',nrad
c	iwk=0
c	Have index ,1 correspond to bottom of layer.
	bmf1(i,1)=bm(1)
	bmf2(i,1)=bm(2)
	bmf3(i,1)=bm(3)
	bmf4(i,1)=bm(4)
	bmf5(i,1)=bm(5)
	do 51 j=1,nrad
	rsav=r
c		write(6,*)'entering derivs5'
	call derivs5(i,r,bm,dydr)
	htot=dr
	eps=3.d-7
	yscal(1)=dabs(bm(1))+fl*dabs(bm(3))+dabs(bm(2))*fl**2+dabs(bm(4))*fl
     &		 +dabs(bm(5))
	yscal(2)=yscal(1)/fl**2
	yscal(3)=yscal(1)/fl
	yscal(4)=yscal(3)
	yscal(5)=yscal(1)
c		write(6,*)'entering bsstep fine'
c		write(6,*)'bm=',bm
c		write(6,*)'dydr=',dydr
c		write(6,*)'r,htot=',r,htot
c		write(6,*)'eps=',eps
c		write(6,*)'yscal=',yscal
	call bsstep(bm,dydr,5,r,htot,eps,i,yscal,hdid,hnext)
	bmf1(i,j+1)=bm(1)
	bmf2(i,j+1)=bm(2)
	bmf3(i,j+1)=bm(3)
	bmf4(i,j+1)=bm(4)
	bmf5(i,j+1)=bm(5)
c		write(6,*)'GETMIN: i,j+1,bmf3=',i,j+1,bmf3(i,j+1)
c		write(6,*)'GETMIN: r=',r,'bm(3)=',bm(3)
c			val=dlog(dabs(bm(4))+(1.d0/l**2)*dabs(bm(1)))-dlog(dabs(bm(2))+
c&			(1.d0/l**2)*dabs(bm(4)))
c			write(8,*) val,dlog(l)
c	if((hdid-htot)/htot.ne.0.d0) iwk=1
	if((hdid-htot)/htot.ne.0.d0) pause "(hdid-htot)/htot.ne.0.d0 in getmin"
c		write(6,*)'r=',r,'bm(5)=',bm(5)
52	r=rsav+dr    
51	continue 
	if(i.eq.n) bmftop=dabs(bm(5))/(dabs(bm(1))+fl*dabs(bm(3))
     &	+dabs(bm(2))*fl**2+dabs(bm(4))*fl)
c		write(6,*)'out of do 51, top layer',i
		if(i.eq.n)write(6,*)'SURFACE: finer (but less accurate!)bm=',bm
c	Use rougher integration step to re-evaluate minor vector at
c	top of layer i.  This rougher step is needed for accuracy because
c	it was used originally in DECAY4M to derive the lvalue.
c	first reset the minor vector to its value at the bottom of layer i.
	bm(1)=bmf1(i,1)
	bm(2)=bmf2(i,1)
	bm(3)=bmf3(i,1)
	bm(4)=bmf4(i,1)
	bm(5)=bmf5(i,1)
c		 write(6,*)'GETMIN: bottom layer',i,'bm=',bm
	nrad=2
	dr=dble(rt(i)-rbot)/dble(nrad)
36	continue
	r=dble(rbot)
c	mus=mu(i)
c	lams=kappa(i)-2.*mu(i)/3.
c	biga=lams+2.*mus 
c	  write(6,*)'entering do 31: l,nrad=',l,nrad
c	  pause
c	Do bm using bsstep.
c		write(6,*)'nrad=',nrad
	iwk=0
c		write(6,*)'entering do 31 (rough integration)'
	do 31 j=1,nrad
	rsav=r
c		write(6,*)'entering derivs5'
	call derivs5(i,r,bm,dydr)
	htot=dr
	eps=3.d-7
	yscal(1)=dabs(bm(1))+fl*dabs(bm(3))+dabs(bm(2))*fl**2+dabs(bm(4))*fl
     &		 +dabs(bm(5))
	yscal(2)=yscal(1)/fl**2
	yscal(3)=yscal(1)/fl
	yscal(4)=yscal(3)
	yscal(5)=yscal(1)
	call bsstep(bm,dydr,5,r,htot,eps,i,yscal,hdid,hnext)
	if((hdid-htot)/htot.ne.0.d0) iwk=1
c		write(6,*)'r=',r,'bm(5)=',bm(5)
32	r=rsav+dr    
31	continue 
c		write(6,*)'GETMIN: rbot,rtop=',rbot,rtop,'nrad=',nrad
c		write(6,*)'GETMIN: top layer',i,' bm=',bm
	if(iwk.ne.0) then
		write(6,*)'*** iwk=1 ***: redoing layer ',i
		bm(1)=bmf1(i,1)
		bm(2)=bmf2(i,1)
		bm(3)=bmf3(i,1)
		bm(4)=bmf4(i,1)
		bm(5)=bmf5(i,1)
		write(6,*)'old nrad=',nrad
		nrad=2*nrad
		write(6,*)'new nrad=',nrad
		dr=dr/2.d0
		go to 36
	endif
38	bmag=real(dabs(bm(1)))+real(l)*real(dabs(bm(3)))
	if(bmag.lt.1.e+5) go to 121
c		write(6,*)'bmag=',bmag,'in layer',i,'RESCALING','bm=',bm
	do 37 j1=1,5
37	bm(j1)=bm(j1)*1.d-10
c-	do 39 ib1=1,201
c-	do 40 ib2=1,51
c-	bmf1(ib1,ib2)=bmf1(ib1,ib2)*1.d-10
c-	bmf2(ib1,ib2)=bmf2(ib1,ib2)*1.d-10
c-	bmf3(ib1,ib2)=bmf3(ib1,ib2)*1.d-10
c-	bmf4(ib1,ib2)=bmf4(ib1,ib2)*1.d-10
c-40	bmf5(ib1,ib2)=bmf5(ib1,ib2)*1.d-10
c-39	continue
	go to 38
121	continue
c		pause
c		i=1
c		if(i.eq.1) stop
c		write(6,*)'GETMIN: bmf3(19,2)=',bmf3(19,2)
	return
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

	subroutine update(radmin,rtop,rbot,nrad,x1)
c	SPHEROIDAL MODES
c	Using integration done in "do 31" loop of main program,
c	interpolate finely the values of eigenfunctions stored in
c	that loop to evaluate the layer integral x1 and to update
c	the recorded eigenfunction values at radii in rw-array.
	real*4 lams,mus,kappa,mu,mu2,mupr
	real*8 dr,r,fnrad,ypoint,evaleq
	real*8 s0,r0
	real*8 sl1(51),sl2(51),sl3(51),sl4(51)
c
	real*8 ya1,ya2,ya3,ya4
cOLD	common/mat/rb(200),rt(200),kappa(200),mu(200),dens(200),l,w0,r0
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0,r0
	common/upd1/rw(500),ur(500),dur(500),vr(500),dvr(500),densi(500),
     &	sur(500)
	common/upd/ya1(51),ya2(51),ya3(51),ya4(51)
c
	third=1./3.
	s0r=real(s0)
	fl21=real(l*(l+1.d0))
	fl12=real((l-1)*(l+2))
	fnrad=dble(nrad)
c	first determine spline interpolation coefficients.
	nrad1=nrad+1
	call splneq(nrad1,ya1,sl1)
	call splneq(nrad1,ya2,sl2)
	call splneq(nrad1,ya3,sl3)
	call splneq(nrad1,ya4,sl4)
c
c	Update eigenfunctions.
	do 38 k=500,1,-1
	if(rw(k).lt.rbot.or.rw(k).gt.rtop) go to 38
	rf=rw(k)
	ypoint=1.d0+dble(rw(k)-rbot)*fnrad/dble(rtop-rbot)
	b1w=real(evaleq(ypoint,nrad1,ya1,sl1))
	b2w=real(evaleq(ypoint,nrad1,ya2,sl2))
	b3w=real(evaleq(ypoint,nrad1,ya3,sl3))
	b4w=real(evaleq(ypoint,nrad1,ya4,sl4))
	r0=dble(rw(k))
	if(r0.lt.dble(radmin)) r0=dble(radmin)	
	call interp(i,mus,lams,biga)
c
	ur(k)=b1w
	dur(k)=(b2w-(lams/rf)*(2.*b1w-fl21*b3w))/biga
	vr(k)=b3w
	dvr(k)=(1./rf)*(b3w-b1w)+b4w/mus 
	sur(k)=b2w
38	continue  
c
	x1=0.
	nrad10=10*nrad
        dr=-dble(rtop-rbot)/dble(nrad10)
	r=dble(rtop)-dr/2.d0
	do 31 j=1,nrad10
	r=r+dr
	ypoint=1.d0+(r-dble(rbot))*fnrad/dble(rtop-rbot)
	b1w=real(evaleq(ypoint,nrad1,ya1,sl1))
	b2w=real(evaleq(ypoint,nrad1,ya2,sl2))
	b3w=real(evaleq(ypoint,nrad1,ya3,sl3))
	b4w=real(evaleq(ypoint,nrad1,ya4,sl4))
	rf=real(r)
	r2=rf*rf
	r0=r
	if(r0.lt.dble(radmin)) r0=dble(radmin)	
	call interp(i,mus,lams,biga)
	tau1=mu(i)/eta1(i)
	tau2=mu2(i)/eta(i)
c	mus=mu(i)*s0r*(s0r+tau2)/((s0r+tau2)*(s0r+tau1)+mu(i)*s0r/eta(i))
c	lams=dble(kappa(i)-2.*mus/3.)
	bigde=(s0r+tau1)*(s0r+tau2)+mu(i)*s0r/eta(i)
	pdds=(s0r+tau1)+(s0r+tau2)+mu(i)/eta(i)
c	dmds=dble((mu(i)-mupr(i))*(mu(i)/eta(i)))/(sfac*sfac)
	dmds=mus/s0r + mus/(s0r+tau2) -mus*pdds/bigde
c	biga=lams+2.*mus
c
	up=(b2w-(lams/rf)*(2.*b1w-fl21*b3w))/biga
	vp=(1./rf)*(b3w-b1w)+b4w/mus
	bigf=(1./rf)*(2.*b1w-fl21*b3w)
	x1=x1+(third*(2.*up-bigf)**2+fl21*((rf*vp-b3w+b1w)**2+
     &	fl12*b3w**2)/r2)*r2*dmds
c		write(6,*)'r=',r,'x1/u**2=',x1/(b1w**2+fl21*b3w**2),'dmds=',dmds
31	continue	
c	x1 must be multiplied by dr, which is a factor of 10 smaller
c	here than in main program.  Main program multiplies by
c	the larger dr, so divide x1 by 10 in this subroutine.
	x1=x1/10.
	return
	end

	subroutine interp(i,mus,lams,biga)
c	Give back mus,lams,biga at radius r0 interpolated between layers.
c	Also identify the layer #i that r0 occupies.
	real*4 mus,lams,kappa,mu,mu2,mupr
	real*8 s0,r0
	common/mat/rb(200),rt(200),kappa(200),mu(200),mu2(200),eta(200),
     &	dens(200),eta1(200),l,mupr(200)
	common/s0val/s0,r0
	s0r=real(s0)
	j=0
5	j=j+1
	if(rt(j).lt.r0) go to 5
	tau1=mu(j)/eta1(j)
	tau2=mu2(j)/eta(j)
	mus=mu(j)*s0r*(s0r+tau2)/((s0r+tau2)*(s0r+tau1)+mu(j)*s0r/eta(j))
	lams=kappa(j)-2.*mus/3.
	biga=lams+2.*mus   
	i=j
	return
	end

