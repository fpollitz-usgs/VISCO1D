c	strainA.f 
c***	THIS VERSION: 3 April 2005.  Response Greens functions
c	at appropriate depths are evaluated with spherical harmonic sum
c	ahead of time, then stored, leading to extremely rapid computation!
c
c	Determine displacement 
c	and horizontal strain fields for summed spheroidal and toroidal
c	modes (post-seismic).
c	Both the vertical and longitudinal components are written out.
c	Read in decay times and displacement coefficients from
c	input files 'decay4.out' and 'vsph.out' and 'decay.out' and 'vtor.out'.
c	Then read in the event and
c	time after earthquake.  
c	This program handles finite faults.  The dip, strike, and rake follow
c	the convention of Ben Menahem and Singh...
c	If the rake is >181 (<-181) deg., a tensile (compressional) dislocation
c	is assumed.
c	The input location of each fault segment 
c	is at the lowermost corner of the fault segment
c	rectangle closest to its strike direction.
c	*** NOTE ***
c	Handles dipping faults with any shear dislocation.
c	This program computes displacements and strains at any number
c	of observation points (up to 15600), and handles faults composed
c	of several segments of varying length and strike (up to
c	2700 segments).  The dip and depth range of each fault segment
c	are fixed.  
c**
c	Output strains in units of 10**(-6) .
c	Second option is to compute strain rate, which is
c	output in units of 10**(-6)/yr (microstrain per year).
c**
c	2.28.97	MODIFIED
c	when strain rate or acceleration is requested (israte=1 or 2),
c	then displ rate or acceleration is also computed.
c**
	parameter (maxmod=21)
	character*80 ahead
	real*4 mrr,mtt,mpp,mrt,mrp,mtp  
cnew
	real*4 rw(40),mudep(40),lamdep(40)
	dimension ur(maxmod,6500),vr(maxmod,6500),dur(maxmod,6500),
     &	dvr(maxmod,6500)
	dimension wr(maxmod,6500),dwr(maxmod,6500)
	real*8 cl,sl,xl0(6501),dxl0(6501),xl1(6501),dxl1(6501)
	real*8 xl2(6501),dxl2(6501)
	real*8 sdecr,cof1r,cof2r
c--
	dimension y1(40*maxmod,6500),y2(40*maxmod,6500),y3(40*maxmod,6500)
	dimension y4(40*maxmod,6500),mdec(6500)
	dimension sdec(maxmod,6500)
c--
	dimension y5(40*maxmod,6500),y6(40*maxmod,6500)
	dimension mdec1(6500),sdec1(maxmod,6500)
c--
	dimension fm(6)
	dimension x1s(6500),x2s(6500),x3s(6500)
	dimension dx1s(6500),dx2s(6500),dx3s(6500)
	dimension ddx1s(6500),ddx2s(6500),ddx3s(6500)
	dimension olat(15600),olon(15600),flat(15600),flon(15600),
     &	fbigl(15600),fstr(15600)
cNEW
	dimension deltr(80)
	dimension rd1S(6,40,80),r1hS(6,40,80),r1tS(6,40,80)
	dimension rttS(6,40,80),rppS(6,40,80),rtrS(6,40,80)
	dimension rtpS(6,40,80),rprS(6,40,80),rrrS(6,40,80)
	dimension romT(6,40,80)
	dimension r1hT(6,40,80),r1tT(6,40,80)
	dimension rttT(6,40,80),rppT(6,40,80),rtrT(6,40,80)
	dimension rtpT(6,40,80),rprT(6,40,80)
	dimension vd1S(6),v1hS(6),v1tS(6),vttS(6),vppS(6)
	dimension vtrS(6),vtpS(6),vprS(6),vrrS(6)
	dimension vomT(6)
	dimension v1hT(6),v1tT(6),vttT(6),vppT(6)
	dimension vtrT(6),vtpT(6),vprT(6)
	real*8 yd1S(80),y1hS(80),y1tS(80),yttS(80),yppS(80)
	real*8 ytrS(80),ytpS(80),yprS(80),yrrS(80)
	real*8 yomT(80)
	real*8 y1hT(80),y1tT(80),yttT(80),yppT(80)
	real*8 ytrT(80),ytpT(80),yprT(80)
	real*8 sarr(80)
	real*8 sd1S(6,40,80),s1hS(6,40,80),s1tS(6,40,80),sttS(6,40,80),
     &	sppS(6,40,80)
	real*8 strS(6,40,80),stpS(6,40,80),sprS(6,40,80),srrS(6,40,80)
	real*8 somT(6,40,80)
	real*8 s1hT(6,40,80),s1tT(6,40,80),sttT(6,40,80),sppT(6,40,80)
	real*8 strT(6,40,80),stpT(6,40,80),sprT(6,40,80)
c---
	dimension fwt(15600),frake(15600)
	real*8 yi1(40),yi2(40),yi3(40),yi4(40),yi5(40),yi6(40)
	real*8 sl1(40),sl2(40),sl3(40),sl4(40),sl5(40),sl6(40),dep,rad1,rad10
	real*8 ypoint,evaleq 
	real*8 wkv1,wkv2,wkv3,wkv4
	real*8 rmin,rmax,deltp
	real*4 len
	common/maxl/lmax
cnew
	common/mpar/u1,aleng,dmax,dmin,sdip,cdip,s2dip,c2dip,frak,srak,crak,
     &	sstr,cstr,s2str,c2str
	common/mtens/mrr,mtt,mpp,mrt,mrp,mtp
	parameter (dstmax=2700.)
c	dstmax has the maximum fault-observation pt. distance. 
	pi=3.1415926  
	twopi=2.*3.1415926
c	ethr=6371.
c	ethr=earth's radius in km.
c	The displacement coefficients are evaluated at depths
c	(1,4,7,...,28)* (depfac) km.
c	bigr=673.1
c	bigr=earth radius in 10**6 cm.
	rad=180./3.1415926
c--
	  write(6,*)'reading spheroidal decay times'
	open(4,file='decay4.out')
	rewind(4)
	read(4,21) nd,ethr,depfac
	write(6,21) nd,ethr,depfac
c	Faulting below depth DEPFAC will not be evaluated.
	bigr=ethr/10.
	idec=0
	is=0 
16	idec=idec+1 
	read(4,44,end=20) l,sdecr,cof1r,cof2r
	if(is.eq.0) llast=l
	lmax=l
	if(l.eq.llast) mdec(l)=idec
	if(l.ne.llast) mdec(l)=1
	if(l.eq.llast) sdec(idec,l)=real(sdecr)
	if(is.eq.0) lmin=l  
	if(is.eq.0) is=1
	if(l.eq.llast) go to 16 
	sdec(1,l)=real(sdecr)
	idec=1
	llast=l 
	go to 16
20	close(4)
	write(6,*)'read spheroidal decay times: lmin,lmax=',lmin,lmax
c--
	  write(6,*)'reading toroidal decay times'
	open(4,file='decay.out')
	rewind(4)
	read(4,21) nd,ethr,depfac
	write(6,21) nd,ethr,depfac
c	Faulting below depth DEPFAC will not be evaluated.
	bigr=ethr/10.
21	format(i2,2f10.3)
	idec=0
	is=0 
26	idec=idec+1 
	read(4,41,end=25) l,sdecr
	if(is.eq.0) llast=l
c	write(6,41) l,sdecr 
	lmax=l
	if(l.eq.llast) mdec1(l)=idec
	if(l.ne.llast) mdec1(l)=1
	if(l.eq.llast) sdec1(idec,l)=real(sdecr)
	if(is.eq.0) lmin=l  
	if(is.eq.0) is=1
	if(l.eq.llast) go to 26 
	sdec1(1,l)=real(sdecr)
	idec=1
	llast=l 
	go to 26
25	close(4)
	write(6,*)'read spheroidal decay times: lmin,lmax=',lmin,lmax
c--
	rad1=dble(1.-28.*depfac/ethr)
	rad10=dble(1.-1.*depfac/ethr)
		write(6,*)'rad1,rad10=',rad1,rad10
	read(5,8) ahead
8	format(a80)
	write(6,*)'max depth, min depth (km), dip(deg.)?'
	read(5,*) dmax,dmin,dip
	cdip=cos(dip/rad)
	sdip=sin(dip/rad)
	c2dip=cdip*cdip-sdip*sdip
	dip=cos(dip/rad)/sin(dip/rad)
	rmin=1.d0-dble(dmax/ethr)
	rmax=1.d0-dble(dmin/ethr)
	write(6,*)'year of earthquake, year obs.#1, year obs.#2 (yrs) ?'
	write(6,*)'same line--viscosity multiplier?'
	read(5,*) tm0,tm1,tm2,vmult
c	vmult equivalent to viscosity=vmult*(8*10**18 Pa-s).
	tm=tm1-tm0+0.5*(tm2-tm1)
c	If strain rate or strain acceleration is specified, then
c	it will be evaluated at time tm1-tm0+0.5*(tm2-tm1).  
	tm1=(tm1-tm0)/3.17
	tm2=(tm2-tm0)/3.17
	tm=tm/3.17
	tm5=5./3.17
c 
	write(6,*)'finite fault with iseg # segments.  iseg=?'
	read(5,*) iseg
c	iobs=0
	fleng=0.
	do 39 i=1,iseg
c	write(6,*)'segment #',i,'lat,lon(deg.),length(km)'
c	write(6,*) 'strike(deg.),slip(cm)?'
	read(5,*) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
c	write(6,*) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
	fleng=fleng+fbigl(i)
	flat(i)=(pi/2.-flat(i)/rad)
	flon(i)=flon(i)/rad
	fstr(i)=fstr(i)/rad 
	frake(i)=frake(i)/rad 
39	continue
	write(6,*)'total fault length =',fleng,'km'
c		pause
c 
	write(6,*)'number of observation points?'
	read(5,*) ipts
	do 40 i=1,ipts 
c	write(6,*)'point #',i,'lat,lon(deg.)?'
	read(5,*) olat(i),olon(i)
	olat(i)=(pi/2.-olat(i)/rad)
	olon(i)=olon(i)/rad
40	continue
c 
	write(6,*)'compute strain rate?(yes=1); strain acceleration=2'
	read(5,*) israte 
	write(6,*)'surface (0) or EDEP km observation(.ne.0)?'
	write(6,*)'(EDEP originally specified in input to VSPHDEP)'
	read(5,*) iobs
		write(6,*)'reading in SPHEROIDAL MOTION response functions'
	open(2,file='vsph.out',form='unformatted')
	rewind(2)
41	format(i4,d17.10)
44	format(i4,3d22.15)
c	Below line is old format.
cOLD44	format(i4,3d20.13)
	read(2) nd
	do 50 l=lmin,lmax
	m10=mdec(l)*40
	idec=0
	do 55 i=1,m10
	read(2) ldum,y1(i,l),y2(i,l),y3(i,l),y4(i,l)
c	  write(6,*) ldum,y1(i,l),y2(i,l),y3(i,l),y4(i,l)
	k=i-40*(i/40)
	if(k.ne.0) go to 55
	idec=idec+1
c	read(2) l,ur(idec,l),vr(idec,l),dur(idec,l),dvr(idec,l)
	read(2) ldum,urdum,vrdum,durdum,dvrdum
c		write(6,*)'l,idec=',l,i,'0 km depth: ',ldum,urdum,vrdum,durdum,dvrdum
	read(2) ldum,ur(idec,l),vr(idec,l),dur(idec,l),dvr(idec,l)
c		write(6,*) 'EFAC km depth: ',ldum,ur(idec,l),vr(idec,l),dur(idec,l),dvr(idec,l)
c	EFAC km observation depth in above read statements.
	if(iobs.eq.0) ur(idec,l)=urdum
	if(iobs.eq.0) vr(idec,l)=vrdum
	if(iobs.eq.0) dur(idec,l)=durdum
	if(iobs.eq.0) dvr(idec,l)=dvrdum
55	continue
50	continue
	close(2)
	  write(6,*)'finished reading in SPHEROIDAL MOTION response functions'
c		lmin=20
c		lmax=20
c		i=1
c		if(i.eq.1) stop
c--
		write(6,*)'reading in TOROIDAL MOTION response functions'
	open(2,file='vtor.out',form='unformatted')
	rewind(2)
	read(2) nd 
		write(6,*)'lmin,lmax=',lmin,lmax
	do 53 l=lmin,lmax
	  nm1=mdec1(l)
	nm10=nm1*40
	idec=0
	  kl=l-100*(l/100)
	do 54 i=1,nm10
	read(2) ldum,y5(i,l),y6(i,l)
c	  if(kl.eq.0) write(6,*) l,y5(i,l),y6(i,l)
cN
	k=i-40*(i/40)
	if(k.ne.0) go to 54
	idec=idec+1
	read(2) ldum,wrdum,dwrdum,wr(idec,l),dwr(idec,l)
c	16*depfac km observation depth in above read statement.
	if(iobs.eq.0) wr(idec,l)=wrdum
	if(iobs.eq.0) dwr(idec,l)=dwrdum
	if(l.eq.100) write(6,*)'idec,l=',idec,l
	if(l.eq.100) write(6,*) l,wr(idec,l),dwr(idec,l)
	if(l.eq.100) write(6,*)'----'
54	continue
53	continue
	close(2)
	  write(6,*)'finished reading in TOROIDAL MOTION response functions'
	write(6,*)'finished reading displacement coefficients'
cNOTE
c		lmin=100
c		lmax=100
c	Divide moment tensor by nmesh (# elements in fault subdivision).
c	nmesh=21
c	nmesh=8
c	nmesh=2
c	nmesh1=8
c	Always use nmesh1=40 in this program!
cUSED FOR MOST SUMATRA RUNS	nmesh1=21
c	nmesh1=2
	nmesh1=21
c	nmesh1=11 (use this value for most examples)
c	nmesh1=2
c	nmesh1=13
c * * *
c       NOTE: because of the overwriting of the y1,y2-arrays done in
c       calculating the spline interpolation functions, do not use nmesh1
c       greater than 41.
c * * *
c	nmesh2=4 (for Landers modeling)
c	nmesh2=25 (used this for 1868 and 1857 events)
cUSED FOR MOST SUMATRA RUNS	nmesh2=25
c	nmesh2=18
	nmesh2=25
c	nmesh2=4
c	nmesh2=3
c	nmesh2=5 (for 1892 event)
c	nmesh2=2
c	nmesh2=2
c	nmesh2=16 (use this value for most examples)
c	nmesh2=2
c	nmesh2=11
c	nmesh2=33
c	nmesh2=35
c	nmesh2=8
c	nmesh2=21
c	nmesh1-1=# depth points
c	nmesh2-1=# horizontal length points 
	nmesh0=nmesh1
c	if(nmesh1.eq.2) nmesh0=41
	dh=(rmax-rmin)*dip/real(nmesh0-1)
		write(6,*)'rmax,rmin,dh=',rmax,rmin,dh
c*****	
	write(6,*)'calculating spline interpolation coefficient'
	write(6,*)'and evaluation functions ahead of time'
c	Calculate the required spheroidal motion
c	coefficients at interpolated depths.  It saves a lot of time to
c	do this ahead of time and retrieve these values later.
		write(6,*)'nmesh0=',nmesh0
	do 150 l=lmin,lmax
	do 275 idec=1,mdec(l)
c	Form arrays yi1--yi4  containing eigenfunction
c	values at discrete radius points for use in spline interpolation.
	do j=1,40
	k=40*idec-40+j
	yi1(j)=dble(y1(k,l))
	yi2(j)=dble(y2(k,l))
	yi3(j)=dble(y3(k,l))
	yi4(j)=dble(y4(k,l))
	enddo
	call splneq(40,yi1,sl1)
	call splneq(40,yi2,sl2)
	call splneq(40,yi3,sl3)
	call splneq(40,yi4,sl4)
	bigh=-dh/2.
c		write(6,*)'doing spline coeff., nmesh0=',nmesh0
	idip=0
181	idip=idip+1
	if(idip.eq.nmesh0) go to 275
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	bigh=bigh+dh
	dep=dble(rmin+bigh/dip)
cnew
	rw(idip)=rmin+bigh/dip
	if(dep.lt.rad1.or.dep.gt.rad10) go to 272
c	Note, in above line, that rad1 corresponds to the deepest depth
c	and rad10 correspnds to the shallowest depth of the
c	eigenfunction write-out.  This is opposite
c	to the convention in STAT1A.  This difference only needs to
c	be noted here.
	ypoint=1.d0+(dep-rad1)*39.d0/(rad10-rad1)
c		ypoint=1.d0
c		if(l.eq.10.and.idec.eq.1) write(6,*)'bigh,rmin,dep,rad1,rad10=',bigh,rmin,dep,rad1,rad10	
c		if(l.eq.10.and.idec.eq.1) write(6,*)'kdip,ypoint=',kdip,ypoint
c	write over existing y1,y2,y3,y4,bfacl,bfacm arrays.
	y1(kdip,l)=real(evaleq(ypoint,40,yi1,sl1))
	y2(kdip,l)=real(evaleq(ypoint,40,yi2,sl2))
	y3(kdip,l)=real(evaleq(ypoint,40,yi3,sl3))
	y4(kdip,l)=real(evaleq(ypoint,40,yi4,sl4))	
	go to 181
c	NOTE--with 272--
c	If part of the fault depth range lies outside of the depth range
c	of the displacement eigenfunctions calculated by VSPHDEP or VSPHG,
c	then do not include that fault depth.
272	continue
	y1(kdip,l)=0.
	y2(kdip,l)=0.
	y3(kdip,l)=0.
	y4(kdip,l)=0.
	y5(kdip,l)=0.
	y6(kdip,l)=0.
	go to 181
c	
275	continue
c
	do 325 idec=1,mdec1(l)
c	Form arrays yi1--yi2 containing eigenfunction
c	values at discrete radius points for use in spline interpolation.
	do j=1,40
	k=40*idec-40+j
	yi5(j)=dble(y5(k,l))
	yi6(j)=dble(y6(k,l))
	enddo
	call splneq(40,yi5,sl5)
	call splneq(40,yi6,sl6)
	bigh=-dh/2.
c		write(6,*)'doing spline coeff., nmesh0=',nmesh0
	idip=0
281	idip=idip+1
	if(idip.eq.nmesh0) go to 325
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	bigh=bigh+dh
	dep=dble(rmin+bigh/dip)
cnew
	rw(idip)=rmin+bigh/dip
	if(dep.lt.rad1.or.dep.gt.rad10) go to 372
c	Note, in above line, that rad1 corresponds to the deepest depth
c	and rad10 correspnds to the shallowest depth of the
c	eigenfunction write-out.  This is opposite
c	to the convention in STAT1A.  This difference only needs to
c	be noted here.
	ypoint=1.d0+(dep-rad1)*39.d0/(rad10-rad1)
c	write over existing y1,y2 arrays.
	y5(kdip,l)=real(evaleq(ypoint,40,yi5,sl5))
	y6(kdip,l)=real(evaleq(ypoint,40,yi6,sl6))
	go to 281
c	NOTE--with 372--
c	If part of the fault depth range lies outside of the depth range
c	of the displacement eigenfunctions calculated by stat0A,
c	then do not include that fault depth.
372	continue
	y5(kdip,l)=0.
	y6(kdip,l)=0.
	go to 281
c	
325	continue
150	continue
cnew
c*
	call eladep(rw,nmesh0,mudep,lamdep)
c *****
c	Next, calculate response Greens functions for m=0, 1, and 2 at 
c	specific values of DELTA.
	write(6,*)'calculate response Greens functions for m=0, 1, and 2 at' 
	write(6,*)'specific values of DELTA'
	do i=1,80
	deltr(i)=dstmax**(0.0125*real(i))/ethr
	enddo
	do 105 i=1,80
	deltpf=deltr(i)
	write(6,*)'i=',i,' DELTA=',deltpf
	cotd=cos(deltpf)/sin(deltpf)
	slf=sin(deltpf)
	call lgndr0(deltpf,xl0,dxl0)
	call lgndr1(deltpf,xl1,dxl1)
	call lgndr2(deltpf,xl2,dxl2) 
c	Store Legendre functions for later use.
	do 364 lk=lmin,lmax+1
	faclk=sqrt(real(lk*(lk+1)))
cc	Note multiplication by [lste] because of the stepped l-summation.
	x1=real(xl0(lk))
	x2=real(xl1(lk))
	x3=real(xl2(lk))
	dx1=real(dxl0(lk))
	dx2=real(dxl1(lk))
	dx3=real(dxl2(lk))
	ddx1=(-faclk*faclk)*x1-cotd*dx1
	ddx2=(1./(slf*slf)-faclk*faclk)*x2-cotd*dx2 
	ddx3=(4./(slf*slf)-faclk*faclk)*x3-cotd*dx3
	x1s(lk)=x1
	x2s(lk)=x2
	x3s(lk)=x3
	dx1s(lk)=dx1
	dx2s(lk)=dx2
	dx3s(lk)=dx3
	ddx1s(lk)=ddx1
	ddx2s(lk)=ddx2
	ddx3s(lk)=ddx3
364	continue 
	do 107 im=1,6
	mrr=0.
	mtt=0.
	mpp=0.
	mrt=0.
	mrp=0.
	mtp=0.
	if(im.eq.1) mrr=1.
	if(im.eq.2) then
	  mtt=0.5
	  mpp=0.5
	endif
	if(im.eq.3) mrt=1.
	if(im.eq.4) mrp=1.
	if(im.eq.5) then
	  mtt=0.5
	  mpp=-0.5
	endif
	if(im.eq.6) mtp=1.
	do 106 idip=1,nmesh0-1
	rd1S(im,idip,i)=0.
	r1hS(im,idip,i)=0.
	r1tS(im,idip,i)=0.
	rttS(im,idip,i)=0.
	rppS(im,idip,i)=0.
	rtrS(im,idip,i)=0.
	rtpS(im,idip,i)=0.
	rprS(im,idip,i)=0.
	rrrS(im,idip,i)=0.
	romT(im,idip,i)=0.
	r1hT(im,idip,i)=0.
	r1tT(im,idip,i)=0.
	rttT(im,idip,i)=0.
	rppT(im,idip,i)=0.
	rtrT(im,idip,i)=0.
	rtpT(im,idip,i)=0.
	rprT(im,idip,i)=0.
c		if(im.ne.6) go to 106
	cp=1.
	sp=1.
	cp2=1.
	sp2=1.
c		write(6,*)'entering do 71: lmin,lmax=',lmin,lmax
	do 371 l=lmin,lmax
c		if(l.lt.325.or.l.gt.350) go to 71
	const=sqrt((real(l)+0.5)/twopi)
	facl=sqrt(real(l*(l+1)))
	facl1=facl*sqrt(real((l-1)*(l+2))) 
c		write(6,*)'entering do 375: mdec(l)=',mdec(l)
	do 375 idec=1,mdec(l)
	s0=sdec(idec,l)
c	Do not include very long period modes since these are not all
c	included in the data set.
c	dfac=1.-exp(-s0*tm)
c	Form difference between observation #2 and observation #1.
	ax1=-s0*tm1/vmult
	ax2=-s0*tm2/vmult
	if(ax1.lt.-20.0) ax1=-20.
	if(ax2.lt.-20.0) ax2=-20.
	dfac=exp(ax1)-exp(ax2)
	dfac1=dfac 
	ax3=-s0*tm/vmult
	if(ax3.lt.-20.0) ax3=-20.
	if(israte.eq.1) dfac1=(s0/vmult)*exp(ax3)/3.17
	if(israte.eq.2) dfac1=-((s0/vmult)/3.17)*
     &	((s0/vmult)/3.17)*exp(ax3)
	dfac=dfac1
c
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	wkv1=dble(y1(kdip,l))
	wkv2=dble(y2(kdip,l))
	wkv3=dble(y3(kdip,l))
	wkv4=dble(y4(kdip,l))
c
c * *	Moment tensor excitation
	z1=const*(mrr*real(wkv1)+0.5*(mpp+mtt)*
     &	real(2.*wkv2-facl*facl*wkv4))
	z2=-const*facl*mrt*real(wkv3+(wkv2-wkv4))
	z3=-const*facl*mrp*real(wkv3+(wkv2-wkv4))
	z4=const*facl1*0.5*(mtt-mpp)*real(wkv4)
	z5=const*facl1*mtp*real(wkv4)
c * *
c		write(6,*)'z1,z2,z3,z4,z5=',z1,z2,z3,z4,z5
	dz1h=z1*dvr(idec,l)
	dz2h=z2*dvr(idec,l)
	dz3h=z3*dvr(idec,l)
	dz4h=z4*dvr(idec,l)
	dz5h=z5*dvr(idec,l)
	dz1=z1*dur(idec,l)
	dz2=z2*dur(idec,l)
	dz3=z3*dur(idec,l)
	dz4=z4*dur(idec,l)
	dz5=z5*dur(idec,l)
	z1h=z1*vr(idec,l)
	z2h=z2*vr(idec,l)
	z3h=z3*vr(idec,l)
	z4h=z4*vr(idec,l)
	z5h=z5*vr(idec,l)
	z1=z1*ur(idec,l)
	z2=z2*ur(idec,l)
	z3=z3*ur(idec,l)
	z4=z4*ur(idec,l)
	z5=z5*ur(idec,l)
	x1=x1s(l)
	x2=x2s(l)
	x3=x3s(l)
	dx1=dx1s(l)
	dx2=dx2s(l)
	dx3=dx3s(l)
	ddx1=ddx1s(l)
	ddx2=ddx2s(l)
	ddx3=ddx3s(l)
c 
	rd1S(im,idip,i)=rd1S(im,idip,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac
	r1hS(im,idip,i)=r1hS(im,idip,i)+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac
	r1tS(im,idip,i)=r1tS(im,idip,i)+(-z2h*x2*sp
     &	+z3h*x2*cp -2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac/slf 
	rttS(im,idip,i)=rttS(im,idip,i)+(z1h*ddx1 +z2h*ddx2*cp
     &	+z3h*ddx2*sp +z4h*ddx3*cp2 +z5h*ddx3*sp2)*dfac1/bigr 
	rttS(im,idip,i)=rttS(im,idip,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	rppS(im,idip,i)=rppS(im,idip,i)+(-z2h*x2*cp -z3h*x2*sp 
     &	-4.*z4h*x3*cp2 -4.*z5h*x3*sp2)*(dfac1/slf)/(slf*bigr) 
	rppS(im,idip,i)=rppS(im,idip,i)+(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1*cotd/bigr 
	rppS(im,idip,i)=rppS(im,idip,i)+(z1*x1 +z2*x2*cp
     &	+z3*x2*sp +z4*x3*cp2 +z5*x3*sp2)*dfac1/bigr 
	rtrS(im,idip,i)=rtrS(im,idip,i)+(dz1h*dx1 +dz2h*dx2*cp
     &	+dz3h*dx2*sp +dz4h*dx3*cp2 +dz5h*dx3*sp2)*dfac1/(2.*bigr) 
	rtrS(im,idip,i)=rtrS(im,idip,i)-(z1h*dx1 +z2h*dx2*cp
     &	+z3h*dx2*sp +z4h*dx3*cp2 +z5h*dx3*sp2)*dfac1/(2.*bigr) 
	rtrS(im,idip,i)=rtrS(im,idip,i)+(z1*dx1 +z2*dx2*cp
     &	+z3*dx2*sp +z4*dx3*cp2 +z5*dx3*sp2)*dfac1/(2.*bigr) 
	rtpS(im,idip,i)=rtpS(im,idip,i)+(-z2h*dx2*sp
     &	+z3h*dx2*cp -2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*bigr*slf) 
	rtpS(im,idip,i)=rtpS(im,idip,i)-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rtpS(im,idip,i)=rtpS(im,idip,i)+(-z2h*dx2*sp +z3h*dx2*cp 
     &	-2.*z4h*dx3*sp2 +2.*z5h*dx3*cp2)*dfac1/(2.*slf*bigr)
	rprS(im,idip,i)=rprS(im,idip,i)+(-z2*x2*sp
     &	+z3*x2*cp -2.*z4*x3*sp2 +2.*z5*x3*cp2)*dfac1/(2.*slf*bigr) 
	rprS(im,idip,i)=rprS(im,idip,i)+(-dz2h*x2*sp +dz3h*x2*cp 
     &	-2.*dz4h*x3*sp2 +2.*dz5h*x3*cp2)*dfac1/(2.*bigr*slf) 
	rprS(im,idip,i)=rprS(im,idip,i)-(-z2h*x2*sp +z3h*x2*cp 
     &	-2.*z4h*x3*sp2 +2.*z5h*x3*cp2)*dfac1/(2.*bigr*slf) 
	rrrS(im,idip,i)=rrrS(im,idip,i)+(dz1*x1 +dz2*x2*cp
     &	+dz3*x2*sp +dz4*x3*cp2 +dz5*x3*sp2)*dfac1/bigr 
c	Done with Spheroidal motion component of degree l.
375	continue
c
	do 75 idec=1,mdec1(l)
cNOTE
c--
	s0=sdec1(idec,l)
c	dfac=1.-exp(-s0*tm)
c	Form difference between observation #2 and observation #1.
	ax1=-s0*tm1/vmult
	ax2=-s0*tm2/vmult
	if(ax1.lt.-20.0) ax1=-20.
	if(ax2.lt.-20.0) ax2=-20.
	dfac=exp(ax1)-exp(ax2)
	dfac1=dfac 
	ax3=-s0*tm/vmult
	if(ax3.lt.-20.0) ax3=-20.
	if(israte.eq.1) dfac1=(s0/vmult)*exp(ax3)/3.17
	if(israte.eq.2) dfac1=-((s0/vmult)/3.17)*
     &	((s0/vmult)/3.17)*exp(ax3)
	dfac=dfac1
c
	kdip=(nmesh0-1)*idec-(nmesh0-1)+idip
	wkv1=dble(y5(kdip,l))
	wkv2=dble(y6(kdip,l))
c
	z1=-const*facl*mrp*real(wkv1)
	z2=const*facl*mrt*real(wkv1)
	z3=const*facl1*mtp*real(wkv2)
	z4=-const*facl1*0.5*(mtt-mpp)*real(wkv2)
c----
	dz1=z1*dwr(idec,l)
	dz2=z2*dwr(idec,l)
	dz3=z3*dwr(idec,l)
	dz4=z4*dwr(idec,l)
	z1=z1*wr(idec,l)
	z2=z2*wr(idec,l)
	z3=z3*wr(idec,l)
	z4=z4*wr(idec,l)
	x2=x2s(l)
	x3=x3s(l)
	dx2=dx2s(l)
	dx3=dx3s(l)
	ddx2=ddx2s(l)
	ddx3=ddx3s(l)
c
	romT(im,idip,i)=romT(im,idip,i)-facl*facl*((-z1*cp-z2*sp)*x2
     &	-(z3*cp2+z4*sp2)*x3)*dfac/(2.*bigr)
	r1tT(im,idip,i)=r1tT(im,idip,i)+(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac
	r1hT(im,idip,i)=r1hT(im,idip,i)-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac/slf 
	rttT(im,idip,i)=rttT(im,idip,i)-(-z1*dx2*sp +z2*dx2*cp
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf) 
	rttT(im,idip,i)=rttT(im,idip,i)+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rppT(im,idip,i)=rppT(im,idip,i)+(-z1*dx2*sp +z2*dx2*cp 
     &	-2.*z3*dx3*sp2 +2.*z4*dx3*cp2)*dfac1/(bigr*slf)
	rppT(im,idip,i)=rppT(im,idip,i)-(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1*cotd/(bigr*slf) 
	rtpT(im,idip,i)=rtpT(im,idip,i)+(z1*ddx2*cp +z2*ddx2*sp
     &	+z3*ddx3*cp2 +z4*ddx3*sp2)*dfac1/(2.*bigr)
	rtpT(im,idip,i)=rtpT(im,idip,i)-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1*cotd/(2.*bigr)
	rtpT(im,idip,i)=rtpT(im,idip,i)-(-z1*x2*cp -z2*x2*sp
     &	-4.*z3*x3*cp2 -4.*z4*x3*sp2)*(dfac1/slf)/(2.*bigr*slf)
	rtrT(im,idip,i)=rtrT(im,idip,i)-(-dz1*x2*sp +dz2*x2*cp
     &	-2.*dz3*x3*sp2 +2.*dz4*x3*cp2)*dfac1/(2.*bigr*slf) 
	rtrT(im,idip,i)=rtrT(im,idip,i)+(-z1*x2*sp +z2*x2*cp
     &	-2.*z3*x3*sp2 +2.*z4*x3*cp2)*dfac1/(2.*bigr*slf) 
	rprT(im,idip,i)=rprT(im,idip,i)+(dz1*dx2*cp +dz2*dx2*sp
     &	+dz3*dx3*cp2 +dz4*dx3*sp2)*dfac1/(2.*bigr)
	rprT(im,idip,i)=rprT(im,idip,i)-(z1*dx2*cp +z2*dx2*sp
     &	+z3*dx3*cp2 +z4*dx3*sp2)*dfac1/(2.*bigr)
c	Done with Toroidal motion component of degree l.
75	continue
c
371	continue
106	continue
107	continue
105	continue
	write(6,*)'response Greens functions determined'
c *****
c	Set up spline interpolation arrays for the response functions.
	do 111 im=1,6
	write(6,*)'setting up spline interpolation arrays for response functions'
	do 112 idip=1,nmesh0-1
	do 113 i=1,80
	yd1S(i)=dble(rd1S(im,idip,i))
	y1hS(i)=dble(r1hS(im,idip,i))
	y1tS(i)=dble(r1tS(im,idip,i))
	yttS(i)=dble(rttS(im,idip,i))
	yppS(i)=dble(rppS(im,idip,i))
	ytrS(i)=dble(rtrS(im,idip,i))
	ytpS(i)=dble(rtpS(im,idip,i))
	yprS(i)=dble(rprS(im,idip,i))
	yrrS(i)=dble(rrrS(im,idip,i))
	yomT(i)=dble(romT(im,idip,i))
	y1hT(i)=dble(r1hT(im,idip,i))
	y1tT(i)=dble(r1tT(im,idip,i))
	yttT(i)=dble(rttT(im,idip,i))
	yppT(i)=dble(rppT(im,idip,i))
	ytrT(i)=dble(rtrT(im,idip,i))
	ytpT(i)=dble(rtpT(im,idip,i))
	yprT(i)=dble(rprT(im,idip,i))
113	continue
c		if(im.eq.1) write(6,*)'y1hS=',y1hS
c
	call splneq(80,yd1S,sarr)
	do i=1,80
	sd1S(im,idip,i)=sarr(i)
	enddo
	mu=real(evaleq(ypoint,16,ymu,smu))
	call splneq(80,y1hS,sarr)
	do i=1,80
	s1hS(im,idip,i)=sarr(i)
c		if(im.eq.1) write(6,*)'i=',i,'s1hS=',s1hS(im,idip,i)
	enddo
	call splneq(80,y1tS,sarr)
	do i=1,80
	s1tS(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yttS,sarr)
	do i=1,80
	sttS(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yppS,sarr)
	do i=1,80
	sppS(im,idip,i)=sarr(i)
	enddo
	call splneq(80,ytrS,sarr)
	do i=1,80
	strS(im,idip,i)=sarr(i)
	enddo
	call splneq(80,ytpS,sarr)
	do i=1,80
	stpS(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yprS,sarr)
	do i=1,80
	sprS(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yrrS,sarr)
	do i=1,80
	srrS(im,idip,i)=sarr(i)
	enddo
c
	call splneq(80,yomT,sarr)
	do i=1,80
	somT(im,idip,i)=sarr(i)
	enddo
	call splneq(80,y1hT,sarr)
	do i=1,80
	s1hT(im,idip,i)=sarr(i)
c		if(idip.eq.1) write(6,*)'s1hT(',im,idip,i,')=',s1hT(im,idip,i)
	enddo
	call splneq(80,y1tT,sarr)
	do i=1,80
	s1tT(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yttT,sarr)
	do i=1,80
	sttT(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yppT,sarr)
	do i=1,80
	sppT(im,idip,i)=sarr(i)
	enddo
	call splneq(80,ytrT,sarr)
	do i=1,80
	strT(im,idip,i)=sarr(i)
	enddo
	call splneq(80,ytpT,sarr)
	do i=1,80
	stpT(im,idip,i)=sarr(i)
	enddo
	call splneq(80,yprT,sarr)
	do i=1,80
	sprT(im,idip,i)=sarr(i)
	enddo
c
112	continue
111	continue
	write(6,*)'done determining spline interpolations of response functions'
c		pause
c
	open(2,file='strainA.out')
c	Begin sweep over observation points.
	do 60 ip=1,ipts 
	dispx=0.
	dispy=0.
	dispz=0.
	exx=0.
	eyy=0.
	exy=0.
	exz=0.
	eyz=0.
	ezz=0.
	omxy=0.
c	Compute displacements, strains and gravity anomaly 
c	at points specified by
c	arrays olat and olon.  Begin summation over fault segments.
	jf=0
65	jf=jf+1
	if(jf.gt.iseg) go to 84
c		write(6,*)'jf=',jf
	write(6,*) flat(jf),flon(jf),fbigl(jf),fstr(jf)
	sstr=sin(fstr(jf))
	cstr=cos(fstr(jf))
	s2str=2.*sstr*cstr
	c2str=cstr*cstr-sstr*sstr 
	frak=frake(jf)
	srak=sin(frake(jf))
	crak=cos(frake(jf))
	dlen=(fbigl(jf)/ethr)/real(nmesh2-1)
	dlon=olon(ip)-flon(jf)
c	Find angular distance and azimuth to observation point from
c	the initial point on the fault segment.
c	Note: delta is the angular distance from the earthquake.
c	Note: phi is the azimuth meazsured positive counterclockwise
c	from south.
	cdelt=cos(flat(jf))*cos(olat(ip))+sin(flat(jf))*
     &	sin(olat(ip))*cos(dlon)
	if (cdelt.le.0.9999) delta=acos(cdelt)
	if(cdelt.gt.0.9999) delta=sqrt((flat(jf)-olat(ip))**2+
     &	(dlon*sin(flat(jf)))**2)
	spsi=sin(dlon)*sin(olat(ip))/sin(delta)
	cpsi=(cos(olat(ip))-cos(flat(jf))*cdelt)/(sin(flat(jf))*sin(delta))
	phi=pi-atan2(spsi,cpsi) 
	dwrite=rad*delta
	pwrite=rad*phi
	  write(6,*)'delta,phi=',dwrite,'deg.',pwrite 
	cphi=cos(phi)
	sphi=sin(phi)
cNEW
	spsi=sin(dlon)*sin(flat(jf))/sin(delta)
	cpsi=(cos(flat(jf))-cos(olat(ip))*cdelt)/(sin(olat(ip))*sin(delta))
	peps=atan2(spsi,cpsi)
	cphio=cos(peps)
	sphio=sin(peps)
c--
	if(jf.gt.1) go to 59
	cphi1=cphi
	sphi1=sphi
	delt1=delta 
59	delta=delta*(6371./ethr)
c	Angular distances are a factor of (6371./ethr) larger on the
c	earth with radius ethr km.
	u1=fwt(jf)
	  write(6,*)'after 59, delta=',delta
c	u1=magnitude of slip.
c	Integrate displacements over fault elements.
	bigh=-dh/2.
	idip=0
81	idip=idip+1
	if(idip.eq.nmesh0) go to 65
c	DO DEPTH INTEGRAL.
c	  write(6,*)'idip=',idip
	bigh=bigh+dh
	dep=dble(rmin+bigh/dip)
	ypoint=1.d0+(dep-rad1)*39.d0/(rad10-rad1)
c	  write(6,*)'dep,rad1,rad10=',dep,rad1,rad10
c	  write(6,*)'ypoint=',ypoint
c ***	Determine moment tensor elements for fault element (ilen,idip).
c	Units of moment tensor are 10**20 N-m.
c	First get [amu] and [alam] for fault element at radius [dep].
c	Determine moment tensor for this fault patch.
	amu=mudep(idip)
	alam=lamdep(idip)
	aleng=fbigl(jf)
	call momten(amu,alam,nmesh1,nmesh2)
c	  write(6,*)'amu,alam=',amu,alam
c	  write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp='
c	  write(6,*) mrr,mtt,mpp,mrt,mrp,mtp
c ***
c	Determine multiplying factors for the cases m=0,1,2 
	fm(1)=mrr
	fm(2)=mtt+mpp
	fm(3)=mrt
	fm(4)=mrp
	fm(5)=mtt-mpp
	fm(6)=mtp
c
	len=-dlen/2.
	ilen=0
83	ilen=ilen+1
	if(ilen.eq.nmesh2) go to 81
c	  write(6,*)'ilen=',ilen
	len=len+dlen
	bigx=delta*sphi+len*sstr+bigh*cstr
	bigy=-delta*cphi+len*cstr-bigh*sstr
	deltp=dble(sqrt(bigx*bigx+bigy*bigy))
	phip=atan2(bigy,bigx)+pi/2.	
	  deltpw=real(deltp)*rad
	  phipw=real(phip)*rad 
c	  write(6,*)'bigx,bigy,deltp,phip=',bigx,bigy,deltp,phipw 
	cp=cos(phip)
	sp=sin(phip)
	cp2=cos(2.*phip)
	sp2=sin(2.*phip)
	cl=dcos(deltp)
	sl=dsin(deltp)
	cotd=real(cl/sl)
	slf=real(sl)
	deltpf=real(deltp)
c*****
c	Now interpolate response functions at this deltpf.
	ypoint=dble(80.*log(deltpf*ethr)/log(dstmax))
	do 114 im=1,6
	do i=1,80
	sarr(i)=sd1S(im,idip,i)
	yd1S(i)=dble(rd1S(im,idip,i))
	enddo
	vd1S(im)=real(evaleq(ypoint,80,yd1S,sarr))
	do i=1,80
	sarr(i)=s1hS(im,idip,i)
	y1hS(i)=dble(r1hS(im,idip,i))
	enddo
	v1hS(im)=real(evaleq(ypoint,80,y1hS,sarr))
	do i=1,80
	sarr(i)=s1tS(im,idip,i)
	y1tS(i)=dble(r1tS(im,idip,i))
	enddo
	v1tS(im)=real(evaleq(ypoint,80,y1tS,sarr))
	do i=1,80
	sarr(i)=sttS(im,idip,i)
	yttS(i)=dble(rttS(im,idip,i))
	enddo
	vttS(im)=real(evaleq(ypoint,80,yttS,sarr))
	do i=1,80
	sarr(i)=sppS(im,idip,i)
	yppS(i)=dble(rppS(im,idip,i))
	enddo
	vppS(im)=real(evaleq(ypoint,80,yppS,sarr))
	do i=1,80
	sarr(i)=strS(im,idip,i)
	ytrS(i)=dble(rtrS(im,idip,i))
	enddo
	vtrS(im)=real(evaleq(ypoint,80,ytrS,sarr))
	do i=1,80
	sarr(i)=stpS(im,idip,i)
	ytpS(i)=dble(rtpS(im,idip,i))
	enddo
	vtpS(im)=real(evaleq(ypoint,80,ytpS,sarr))
	do i=1,80
	sarr(i)=sprS(im,idip,i)
	yprS(i)=dble(rprS(im,idip,i))
	enddo
	vprS(im)=real(evaleq(ypoint,80,yprS,sarr))
	do i=1,80
	sarr(i)=srrS(im,idip,i)
	yrrS(i)=dble(rrrS(im,idip,i))
	enddo
	vrrS(im)=real(evaleq(ypoint,80,yrrS,sarr))
c
	do i=1,80
	sarr(i)=somT(im,idip,i)
	yomT(i)=dble(romT(im,idip,i))
	enddo
cNOTE Below and elsewhere -- 0.000 factor in order to test
c	spheroidal-only response
	vomT(im)=real(evaleq(ypoint,80,yomT,sarr))
	do i=1,80
	sarr(i)=s1hT(im,idip,i)
	y1hT(i)=dble(r1hT(im,idip,i))
	enddo
	v1hT(im)=real(evaleq(ypoint,80,y1hT,sarr))
	do i=1,80
	sarr(i)=s1tT(im,idip,i)
	y1tT(i)=dble(r1tT(im,idip,i))
	enddo
	v1tT(im)=real(evaleq(ypoint,80,y1tT,sarr))
	do i=1,80
	sarr(i)=sttT(im,idip,i)
	yttT(i)=dble(rttT(im,idip,i))
	enddo
	vttT(im)=real(evaleq(ypoint,80,yttT,sarr))
	do i=1,80
	sarr(i)=sppT(im,idip,i)
	yppT(i)=dble(rppT(im,idip,i))
	enddo
	vppT(im)=real(evaleq(ypoint,80,yppT,sarr))
	do i=1,80
	sarr(i)=strT(im,idip,i)
	ytrT(i)=dble(rtrT(im,idip,i))
	enddo
	vtrT(im)=real(evaleq(ypoint,80,ytrT,sarr))
	do i=1,80
	sarr(i)=stpT(im,idip,i)
	ytpT(i)=dble(rtpT(im,idip,i))
	enddo
	vtpT(im)=real(evaleq(ypoint,80,ytpT,sarr))
	do i=1,80
	sarr(i)=sprT(im,idip,i)
	yprT(i)=dble(rprT(im,idip,i))
	enddo
	vprT(im)=real(evaleq(ypoint,80,yprT,sarr))
114	continue
c*****
cOLD	  gamma1=sp
cOLD	  gamma2=-cp 
cNEW
	bigx=delta*sphio+len*sstr+bigh*cstr
	bigy=-delta*cphio+len*cstr-bigh*sstr
	deltp=dble(sqrt(bigx*bigx+bigy*bigy))
	phip=atan2(bigy,bigx)+pi/2.	
	  gamma1=sin(phip)
	  gamma2=-cos(phip)
c--
	omxy1=vomT(4)*fm(4)*cp+vomT(3)*fm(3)*sp+vomT(6)*fm(6)*cp2
     &	+vomT(5)*fm(5)*sp2
	disp1=vd1S(1)*fm(1)+vd1S(2)*fm(2)+vd1S(3)*fm(3)*cp
     &	+vd1S(4)*fm(4)*sp+vd1S(5)*fm(5)*cp2+vd1S(6)*fm(6)*sp2
	disp1h=v1hS(1)*fm(1)+v1hS(2)*fm(2)+(v1hS(3)+v1hT(3))*fm(3)*cp
     &	+(v1hS(4)+v1hT(4))*fm(4)*sp+(v1hS(5)+v1hT(5))*fm(5)*cp2
     &	+(v1hS(6)+v1hT(6))*fm(6)*sp2
	disp1t=(v1tS(3)+v1tT(3))*fm(3)*sp
     &	+(v1tS(4)+v1tT(4))*fm(4)*cp+(v1tS(5)+v1tT(5))*fm(5)*sp2
     &	+(v1tS(6)+v1tT(6))*fm(6)*cp2
	ett=vttS(1)*fm(1)+vttS(2)*fm(2)+(vttS(3)+vttT(3))*fm(3)*cp
     &	+(vttS(4)+vttT(4))*fm(4)*sp+(vttS(5)+vttT(5))*fm(5)*cp2
     &	+(vttS(6)+vttT(6))*fm(6)*sp2
	epp=vppS(1)*fm(1)+vppS(2)*fm(2)+(vppS(3)+vppT(3))*fm(3)*cp
     &	+(vppS(4)+vppT(4))*fm(4)*sp+(vppS(5)+vppT(5))*fm(5)*cp2
     &	+(vppS(6)+vppT(6))*fm(6)*sp2
	etr=vtrS(1)*fm(1)+vtrS(2)*fm(2)+(vtrS(3)+vtrT(3))*fm(3)*cp
     &	+(vtrS(4)+vtrT(4))*fm(4)*sp+(vtrS(5)+vtrT(5))*fm(5)*cp2
     &	+(vtrS(6)+vtrT(6))*fm(6)*sp2
	etp=(vtpS(3)+vtpT(3))*fm(3)*sp
     &	+(vtpS(4)+vtpT(4))*fm(4)*cp+(vtpS(5)+vtpT(5))*fm(5)*sp2
     &	+(vtpS(6)+vtpT(6))*fm(6)*cp2
	epr=(vprS(3)+vprT(3))*fm(3)*sp
     &	+(vprS(4)+vprT(4))*fm(4)*cp+(vprS(5)+vprT(5))*fm(5)*sp2
     &	+(vprS(6)+vprT(6))*fm(6)*cp2
	err=vrrS(1)*fm(1)+vrrS(2)*fm(2)+vrrS(3)*fm(3)*cp
     &	+vrrS(4)*fm(4)*sp+vrrS(5)*fm(5)*cp2+vrrS(6)*fm(6)*sp2
c
c	Rotate strain tensor into x-y-z Cartesian coordinates.
c	Also rotate displacements into x-y coordinates.
	dispx1=gamma1*disp1h-gamma2*disp1t
	dispy1=gamma2*disp1h+gamma1*disp1t
	dispz1=disp1
	exx1=gamma1**2*ett+gamma2**2*epp-2.*gamma1*gamma2*etp
	exy1=gamma1*gamma2*(ett-epp)+(gamma1**2-gamma2**2)*etp
	eyy1=gamma2**2*ett+gamma1**2*epp+2.*gamma1*gamma2*etp
	exz1=etr*gamma1-epr*gamma2
	eyz1=etr*gamma2+epr*gamma1
	ezz1=err
70	dispx=dispx+dispx1
	dispy=dispy+dispy1
	dispz=dispz+dispz1
	exx=exx+exx1
	eyy=eyy+eyy1
	exy=exy+exy1
	exz=exz+exz1
	eyz=eyz+eyz1
	ezz=ezz+ezz1
	omxy=omxy+omxy1
	go to 83
84	continue
c		write(6,*)'B omxy,dispz,flste=',omxy,dispz,flste
c	Find cartesian corrdinates of observation point relative
c	to fault segment #1.
	xc=(twopi*6371./360.)*rad*(delt1)*sphi1  
	yc=-(twopi*6371./360.)*rad*(delt1)*cphi1  
c	xc measured in east direction; yc measured in north direction (km).
c	displacements in cm.  gravity anomaly in mGal. 
	write(2,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy  
	write(6,76) xc,yc,dispx,dispy,dispz,exx,eyy,exy,exz,eyz,ezz,omxy
76	format(2f10.3,10e13.6e2)
77	format(2f10.3,2e13.6e2)
60	continue
	close(2)
	write(6,*)'end of strainA'
	end  

	subroutine eladep(rw,nmesh0,mudep,lamdep)
	character*80 aread
	real*4 mudep(40),lamdep(40)
	real*4 kappa(200),mu,mupr
	dimension dens(200),rb(200),rt(200),mu(200),eta(200),mupr(200)
	dimension rw(40)
c	Read in earth model
	open(4,file='earth.model')
	rewind(4)
	read(4,5) n,nd,bigr,depfac 
c	write(6,5) n,nd,bigr,depfac
5	format(2i2,2f10.3)
c ** *
	do 10 j=1,n
	i=2*j-1
c
	read(4,14) aread 
	read(aread,16,err=11) rb(i),rt(i),dens(i),kappa(i),mu(i),
     &	mupr(i),eta(i)
	go to 12
11	continue
		write(6,14) aread
	read(aread,15) rb(i),rt(i),dens(i),kappa(i),mu(i),eta(i)
	mupr(i)=0.
12	continue
c	write(6,16) rb(i),rt(i),dens(i),kappa(i),mu(i),mupr(i),eta(i)
c***	Above lines: specify [mu prime] for standard linear solid
c***	in same units as mu.
c***	mu prime = 0 corresponds to Maxwell material.
c***	Reference: Cohen, Time dependent deformation following an earthquake,
c***	JGR, vol. 87, pp. 5414-5421 (1982).
c***
	dens(i+1)=dens(i)
	mu(i+1)=mu(i)
	mupr(i+1)=mupr(i)
	eta(i+1)=eta(i)
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
16	format(6f9.3,e13.6e2)
	n=2*n 
	close(4)
c ** *
	do 30 j=1,nmesh0-1
	do 20 i=1,n
c	if(j.eq.1) write(6,*)'rb(',i,')=',rb(i),'rt(',i,')=',rt(i),'rw(',j,')=',rw(j)
	if(rb(i).lt.rw(j).and.rt(i).ge.rw(j)) then
	  mudep(j)=mu(i)
	  lamdep(j)=kappa(i)-(2./3.)*mu(i)
	endif
20	continue
30	continue
c		write(6,*)'mudep=',mudep
c		write(6,*)'lamdep=',lamdep
c		pause
	return
	end
 
	subroutine momten(amu,alam,nmesh1,nmesh2)
	real*4 mrr,mtt,mpp,mrt,mrp,mtp  
	common/mpar/u1,aleng,dmax,dmin,sdip,cdip,s2dip,c2dip,frak,srak,crak,
     &	sstr,cstr,s2str,c2str
	common/mtens/mrr,mtt,mpp,mrt,mrp,mtp
c**	Determine moment tensor using depth-dependent elastic moduli
c	Use mu=amu and lambda=alam*10**10 Pa.
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
		write(6,*)'strainA-momten: abs(frak)=',abs(frak)
	if(abs(frak).gt.3.159046) go to 49
c	Next line is shear moment.
	shrm1=u1*aleng*((dmax-dmin)/sdip)*amu*1.e-6
	p1=srak*sdip*cdip*s2str+crak*sdip*c2str
	p2=crak*sdip*s2str-srak*sdip*cdip*c2str
	p3=-crak*cdip*sstr+srak*c2dip*cstr
	p4=srak*c2dip*sstr+crak*cdip*cstr
	p5=srak*sdip*cdip
	mrr=shrm1*2.*p5
	mtt=-shrm1*(p2+p5)
	mpp=shrm1*(p2-p5)
	mrt=-shrm1*p4
	mrp=-shrm1*p3
	mtp=-shrm1*p1
	go to 51
c *
49	erak=1.
cOLD	if(srak.lt.-1.74526E-02) erak=-1.
	if(frak.lt.-3.159046) erak=-1.
c	  if(erak.eq.1.0) write(6,*)'tensile dislocation'
c	  if(erak.eq.-1.0) write(6,*)'compressional dislocation'
        shrm1=2.*u1*aleng*((dmax-dmin)/sdip)*erak*amu*1.e-6
	mrr=shrm1*cdip**2
	mtt=shrm1*(sdip*sstr)**2
	mpp=shrm1*(sdip*cstr)**2
	mrt=shrm1*sdip*cdip*sstr
	mrp=shrm1*sdip*cdip*cstr
	mtp=shrm1*sdip**2*sstr*cstr
	shrm2=0.5*shrm1*(alam/amu)
	mrr=mrr+shrm2
	mtt=mtt+shrm2
	mpp=mpp+shrm2
c *
51	continue
	  write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp='
	  write(6,*) mrr,mtt,mpp,mrt,mrp,mtp
	mrr=mrr/real((nmesh1-1)*(nmesh2-1))
	mtt=mtt/real((nmesh1-1)*(nmesh2-1))
	mpp=mpp/real((nmesh1-1)*(nmesh2-1))
	mrt=mrt/real((nmesh1-1)*(nmesh2-1))
	mrp=mrp/real((nmesh1-1)*(nmesh2-1))
	mtp=mtp/real((nmesh1-1)*(nmesh2-1))
c**
c		mrr=0.
c		mtt=0.
c		mpp=1.
c		mrt=0.
c		mrp=0.
c		mtp=0.
	return
	end

