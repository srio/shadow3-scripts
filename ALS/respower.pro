;+
;
;  NAME:
; 	RESPOWER
;  PURPOSE:
; 	to compute the resolving power E/DE from the dispersion in a plane
;       (usually the exit slit plane).
;  CATEGORY:
;        SHADOW's utilities
;  CALLING SEQUENCE:
; 	respower,'myfile',col_E, col_D
;  INPUTS:
; 	myfile  name of the file with data (between quotes) 
;                (it can also be a structure like idl_var)
;       col_E:  the SHADOW column with the energy information: 
;                   11: Energy in eV
;                   19: Wavelength in A
;       col_D:  the SHADOW column with the legth in the dispersion direction (1 or 3)
;  KEYWORD PARAMETERS:
; 		ERANGE=[emin,emax] range of the energy variable for the TOP plot
; 		ZRANGE=[Zmin,Zmax]              Z
; 		E2RANGE=[emin,emax] range of the energy variable for the BOTTOM plot
; 		Z2RANGE=[Zmin,Zmax]              Z
;               NBINS              number of bins for the histogram
; 		HLIMIT             the normalized heigh at which the histogram 
;                                  limits taken for calculating the resolving power. 
;                                  Default: 0.1
; 		TITLE='top_title' title to be written at the top
; 		NOPLOT = When set, inhibits the graph.
;
;  OUTPUTS:
; 	a plot
;  OPTIONAL OUTPUT PARAMETERS:
; 	resolvingPower: Set this keyword to a named variable to receive the Resolving Power
;       deltaEmin: Set this keyword to a named variable to receive the minimum DE 
;                  (at zero exist slit opening)
;       Pendent: Set this keyword to a named variable to receive P
;
;                Note that for any exit slit aperture  DE = deltaE + DZ/P
;  COMMON BLOCKS:
; 	None.
;  SIDE EFFECTS:
; 	None.
;  RESTRICTIONS:
; 	None.
;  PROCEDURE:
; 	Compute the dispersion plot (Z vs DE, top graph) then calculate the dispersion
;       band by 
;         i) Removing the regression line (bottom: plot with regression removed)
;        ii) Calculate Z histogram
;       iii) Calculate the histogram limits at HLIMIT height (bottom plot), 
;        iv) Plot these limits in the upper graphic to define the dispersion band,
;         v) compute related parameters:
;          The "pendent" P from the regression fit gives the disperion DZ/DE
;          The "origin" is the Eo value at the upper boundary of the dispersion band
;          The Resolving Power is  Eo/P 
;          The "Delta" value is the DE corresponding to the histogram limits, or 
;              equivalently, the DE at zero slit opening
;       
;          The energy bandwith at a given slit opening Z is DE= Delta+ Z/P
;
;  MODIFICATION HISTORY:
; 	by M. Sanchez del Rio. ESRF. Grenoble, around 1995
; 	2013/03/18 srio@esrf.eu documented, extracted some parameters
; 
;-
Pro respower,input,col1,col2,nolost=nolost,nbins=nbins, hlimit = hlimit, title=title, $
   erange=erange, zrange=zrange,z2range=z2range, $
   resolvingPower=resolvingPower,deltaEmin=deltaE,pendent=pendent

;
on_error,2
if n_params() LT 3 then begin
  print,'RESPOWER: Usage: respower,'+"'"+'shadowdata'+"'"+',col_energy,col_x'
  return
endif
if  col1 ne 11 and col1 NE 19 then print, $
  'RESPOWER: Warning: first column is NOT energy or wavelength.'

if n_elements(erange) NE 2 then erange=[0,0] else erange=double(erange)
if n_elements(zrange) NE 2 then zrange=[0,0] else zrange=double(zrange)
;
; load shadow-idl structure and define  arrays and constants
;
str = readsh(input)
if not(keyword_set(str)) then return
if not(keyword_set(title)) then title=''
if not(keyword_set(hlimit)) then hlimit = 0.1d0

a=getshcol(str,[col1,col2],nolost=nolost)

if not(keyword_set(nbins)) then nbins=100L
calfwhm = 1
if not(keyword_set(gauss)) then gauss=2
degree=1

coeff=poly_fit(a(0,*),a(1,*),degree,yfit,error,sigma,amatrix)
print,'RESPOWER: the poly_fit parameters are y = c(0) + c(1)x + c(2)x^2 + ...  with: '
for i=0,degree do begin
  print,' c[',i,'] = ',coeff(i)
endfor
;
; substract the fit
;
diff = a
diff[1,*] = a[1,*]-yfit

sst  = stddev(diff[1,*])
mean1  = mean(diff[1,*])
print,'RESPOWER: Mean of residuals: ',mean1
print,'RESPOWER: StDev of residuals: ',sst

arr1=diff[0,*]
arr2=diff[1,*]
if n_elements(z2range) ne 2 then begin
  yran=fltarr(2)
  yran[0] = min(arr2)
  yran[1] = max(arr2)
endif else begin
  yran = z2range
endelse
binsize =(yran[1]-yran[0])/float(nbins)
hy = histogramw(arr2,binsize=binsize,min=yran[0],max=yran[1])
hx = fltarr(nbins)
for i=0,nbins-1 do hx[i] = yran[0]+ binsize/2 + binsize*i
hhx = fltarr (2*nbins)
hhy = fltarr (2*nbins)
for i=0,2*nbins-1,2  do begin
 hhx[i] = hx[fix(i/2)]-binsize/2
 hhx[i+1] = hx[fix(i/2)]+binsize/2
 hhy[i] = hy[fix(i/2)]
 hhy[i+1] = hy[fix(i/2)]
endfor
xrange = [0,1.1*max(hhy)]

; 
; get histo maximum and limits
; 
hmax = max(hy)
hmaxi =  where ( hy EQ hmax)
hmaxi = hmaxi(0)
tmp = where ( hy GT hmax*hlimit)
hlefti = tmp(0)
hrighti = tmp(n_elements(tmp)-1)
print,'Histo tolerance: ',hlimit
print,'Left: ',hx(hlefti),hy(hlefti)
print,'Peak: ',hx(hmaxi),hy(hmaxi)
print,'Right: ',hx(hrighti),hy(hrighti)

;
; --------- all plots --------
;
oldpos=!p.position
xx0 = !p.position[0]
yy0 = !p.position[1]
xx1 = !p.position[2]
yy1 = !p.position[3]
if xx0 EQ xx1 then begin
  xx0 = 0.0  &  xx1 = 1.0
endif
if yy0 EQ yy1 then begin
  xx0 = 0.0  &  yy1 = 1.0
endif
xx = abs(xx1-xx0)
yy = abs(yy1-yy0)

help,col1
CASE col1 OF
  11: xtitle='Energy [eV]'
  19: xtitle='Wavelength [A]'
  else: xtitle='col'+StrCompress(col1,/Rem)
ENDCASE
; 
; top graph
;
!p.position= [xx0,yy0,xx0,yy0] + [.1*xx,.60*yy,.7*xx,.95*yy]
plot,a[0,*],a[1,*],psym=3,xtitle=xtitle,ytitle='Z [cm]',xrange=erange,yrange=zrange
amin = min(a[0,*])
amax = max(a[0,*])
xfit=amin + findgen(100)/99.*(amax-amin)
oplot,xfit,coeff[1]*xfit+coeff[0]+hx(hmaxi)
;
; overplots fit shifted to the histogram center
;
oplot,xfit,xfit*0.
oplot,xfit,coeff[1]*xfit+coeff[0]+hx[hlefti],linestyle=3
oplot,xfit,coeff[1]*xfit+coeff[0]+hx[hrighti],linestyle=3
orig = -1.0*(coeff[0]+hx[hmaxi])/coeff[1]

deltax1 = abs( (hx[hlefti]-hx[hmaxi])/coeff[1])
IF (col1 EQ 19) THEN BEGIN
 oplot,[orig+deltax1,orig+deltax1],[-1000,1000]
ENDIF ELSE BEGIN
 oplot,[orig-deltax1,orig-deltax1],[-1000,1000]
ENDELSE

deltax2 = abs( (hx[hrighti]-hx[hmaxi])/coeff[1])
IF (col1 EQ 19) THEN BEGIN
  oplot,[orig-deltax2,orig-deltax2],[-1000,1000]
ENDIF ELSE BEGIN
  oplot,[orig+deltax2,orig+deltax2],[-1000,1000]
ENDELSE
deltaE = deltax1 + deltax2
;
; bottom graph
;
!p.position= [xx0,yy0,xx0,yy0] + [.1*xx,.1*yy,.7*xx,.50*yy]
plot,diff[0,*],diff[1,*],psym=3,/noerase,xtitle=xtitle, $
  ytitle='Z-Z!Dfit!N [cm]',xrange=erange,yrange=z2range
;
; overplots lines shifted to the histogram center, ans left and right sides
;
oplot,xfit,xfit*0.+hx(hmaxi)
oplot,xfit,xfit*0.+hx(hlefti),linestyle=3
oplot,xfit,xfit*0.+hx(hrighti),linestyle=3

;
; titles
;
csize = !p.charsize * 1.4
tsize = !p.charsize * 2.0
xyouts,/nor,xx0+xx*0.10,yy0+yy*0.98,title,siz=tsize
text='Command: respower,'+"'"+str.name+"',"+ $
 strcompress(col1,/rem)+','+strcompress(col2,/rem)+',nbins='+ $
 strcompress(nbins,/rem)
if keyword_set(hlimit) then text = text + ',hlimit='+strcompress(hlimit,/rem)
if keyword_set(nolost) then text = text + ',nolost='+strcompress(nolost,/re)
xyouts,/nor,xx0+xx*0.10,yy0+yy*0.96,text

;
; data
;
resolvingPower = orig/deltaE
pendent = coeff[1]
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.9,'Linear fit pendent P: '+strcompress(pendent,/rem)
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.87,'Histogram peak at: '+strcompress(hx[hmaxi],/rem)
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.84,'Histogram base line: '+strcompress(hlimit,/rem)
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.78,'DeltaE(DZ=0) = '+strcompress(deltaE,/rem)
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.75,'Origin E = '+strcompress(orig,/rem)
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.72,'Resolving Power = '+strcompress(resolvingPower,/rem)
intens,str,nolost=nolost,ii
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.69,'Intensity = '+strcompress(ii,/rem)
xyouts,/nor,xx0+xx*0.71,yy0+yy*0.66,'DE ~ DeltaE + DZ/P'

;
; histogram plot
;
!p.position= [xx0,yy0,xx0,yy0] + [.75*xx,.1*yy,.9*xx,.5*yy]
plot,hhy,hhx,yrange=yran,xrange = xrange,/noerase
plot,hhy*0.+hmax*hlimit,hhx,yrange=yran,xrange = xrange,/noerase

;
; reset !p.position
;
!p.position=oldpos

end

