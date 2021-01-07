import numpy

# https://www.rxoptics.de/design-parameters/
# view-source:https://www.rxoptics.de/design-parameters/
# https://www.rxoptics.de/wp-content/themes/rxoptics/js/crlcalc.js

def compute(e=10.0,
            desf=1000.0, #focal length (presf=0) in mm or image distance (presf=1) in mm
            dist=1000.0, #distance source-lens in mm
            R=0.05,
            w=150.0, h=150.0, F=1000.0, d=1000.0, th=30.0,
            presf=0, #0=focal length, 1=image distance,
            MATERIAL="Be", plot_cross_sections=False):


	# var tdelta,mu,lambda,rho,muoverrhos;
    lambda1 = 12.398 / e
    if MATERIAL == "Be":
        rho=1.85
        tdelta=5.4e-6 * lambda1**2 * rho * 4/9.012
        muoverrhos= [[1.0,604.1],[1.5,179.7],[2.0,74.69],[3.0,21.27],[4.0,8.685],[5.0,4.369],[6.0,2.527],[8.0,1.124],[10.0,0.6466],[15.0,0.307],[20.0,0.2251],[30.0,0.1792],[40.0,0.164],[50.0,0.1554],[60.0,0.1493],[80.0,0.1401],[100.0,0.1328],[150.0,0.119],[200.0,0.1089],[300.0,0.09463],[400.0,0.08471],[500.0,0.07739],[600.0,0.07155],[800.0,0.06286],[1000.0,0.05652],[1250.0,0.05054],[1500.0,0.04597],[2000.0,0.03938],[3000.0,0.03138],[4000.0,0.02664],[5000.0,0.02347],[6000.0,0.02121],[8000.0,0.01819],[10000.0,0.01627],[15000.0,0.01361],[20000.0,0.01227]]
    elif MATERIAL == "Al":
        rho=2.7
        tdelta=5.4*1e-6 * lambda1**2 * rho * 13/26.982
        muoverrhos=[[1.0,1185.0],[1.5,402.2],[1.5596,362.1],[1.5596,3957.0],[2.0,2263.0],[3.0,788.0],[4.0,360.5],[5.0,193.4],[6.0,115.3],[8.0,50.33],[10.0,26.23],[15.0,7.955],[20.0,3.441],[30.0,1.128],[40.0,0.5685],[50.0,0.3681],[60.0,0.2778],[80.0,0.2018],[100.0,0.1704],[150.0,0.1378],[200.0,0.1223],[300.0,0.1042],[400.0,0.09276],[500.0,0.08445],[600.0,0.07802],[800.0,0.06841],[1000.0,0.06146],[1250.0,0.05496],[1500.0,0.05006],[2000.0,0.04324],[3000.0,0.03541],[4000.0,0.03106],[5000.0,0.02836],[6000.0,0.02655],[8000.0,0.02437],[10000.0,0.02318],[15000.0,0.02195],[20000.0,0.02168]]
    else:
        rho=8.90
        tdelta=5.4*1e-6 * lambda1**2 * rho * 28/58.693
        muoverrhos=[[1.0,9855.0],[1.00404,9753.0],[1.0081,9654.0],[1.0081,10990.0],[1.5,4234.0],[2.0,2049.0],[3.0,709.4],[4.0,328.2],[5.0,179.3],[6.0,109.0],[8.0,49.52],[8.3328,44.28],[8.3328,329.4],[10.0,209.0],[15.0,70.81],[20.0,32.2],[30.0,10.34],[40.0,4.6],[50.0,2.474],[60.0,1.512],[80.0,0.7306],[100.0,0.444],[150.0,0.2208],[200.0,0.1582],[300.0,0.1154],[400.0,0.09765],[500.0,0.08698],[600.0,0.07944],[800.0,0.06891],[1000.0,0.0616],[1250.0,0.05494],[1500.0,0.05015],[2000.0,0.04387],[3000.0,0.03745],[4000.0,0.03444],[5000.0,0.03289],[6000.0,0.0321],[8000.0,0.03164],[10000.0,0.03185],[15000.0,0.0332],[20000.0,0.03476]]

    muoverrhos = numpy.array(muoverrhos)
    l= muoverrhos.shape[0]


    import xraylib
    # xrl_mu = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber(MATERIAL), e)


    if plot_cross_sections:
        XRL_MU = numpy.zeros(l)
        XRL_MU_E = numpy.zeros(l)

        for i in range(l):
            XRL_MU[i] = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber(MATERIAL), muoverrhos[i,0])
            XRL_MU_E[i] = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber(MATERIAL), muoverrhos[i,0])
            # print("     >>",muoverrhos[i,0],rho*muoverrhos[i,1], XRL_MU[i], XRL_MU_E[i], XRL_MU_E[i] / XRL_MU[i])

        from srxraylib.plot.gol import plot
        plot(muoverrhos[:,0], rho*muoverrhos[:,1],
             muoverrhos[:,0], XRL_MU,
             muoverrhos[:, 0], XRL_MU_E,
             xlog=1,ylog=1,title=MATERIAL,xrange=[1,1e3],xtitle="Photon energy [keV]",ytitle="mu [cm-1]",
             legend=["mu xro","mu xraylib","energy loss mu"])

    print(">>> len l", l)
    for i in range(l):
        if muoverrhos[i, 0] <= e:
            if i == 0:
                mu = rho * muoverrhos[0, 1]
            elif i == (l - 1):
                mu = rho * muoverrhos[l - 1, 1]
            else:
                # mu = rho * (muoverrhos[i,1] * \
                #             (e - muoverrhos[i-1,0]) / (muoverrhos[i, 0] - muoverrhos[i-1, 0]) - \
                #             muoverrhos[i-1, 1] * \
                #             (e - muoverrhos[i,0]) / (muoverrhos[i,0] - muoverrhos[i-1,0]) \
                #             )
                # mu=rho*(muoverrhos[i][1]*(e-muoverrhos[i-1][0])/(muoverrhos[i][0]-muoverrhos[i-1][0])-muoverrhos[i-1][1]*(e-muoverrhos[i][0])/(muoverrhos[i][0]-muoverrhos[i-1][0]));

                mu = rho * (
                        muoverrhos[i, 1] * (e - muoverrhos[i - 1, 0]) / (muoverrhos[i, 0] - muoverrhos[i - 1, 0]) - \
                        muoverrhos[i - 1, 1] * (e - muoverrhos[i, 0]) / (muoverrhos[i, 0] - muoverrhos[i - 1, 0])
                )




    xrl_mu = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber(MATERIAL), e)
    xrl_delta = 1.0 - (xraylib.Refractive_Index("Be", e, rho)).real
    print("e: %g, lambda: %g" % (e,lambda1))
    print(">>>", 5.4e-6 * lambda1**2 * rho * 4/9.012)
    print(">>>>>>>rxo_mu: %g , xrl_mu: %g" % (mu,  xrl_mu))
    print(">>>>>>>rxo_delta: %g , xrl_delta: %g" % (tdelta/2, xrl_delta))  #<<<<<<<<<<< factor 2!!!!


    # //compute desired focal length and real number of lenses
    dist=dist*1e+3 # //dist in mm

    if presf: # //compute focal length
        f=desf
    else:
        f=dist*desf/(dist+desf)

    N = R / (tdelta * F)
    # N = numpy.sqrt(R/(tdelta*F)) * numpy.arcsin( numpy.sqrt(R*F/tdelta) / f)#??????????????????

    print(">>>>f: %f" % (f) )

    if numpy.isnan(N):
        raise Exception("The focus lies within the CRL.")

    # //compute integer number of lenses and true focal length
    N=numpy.floor(N)
    print("N: %f" % (N) )

    exfa = numpy.sqrt(R*F/tdelta) * 1/(numpy.sin(N    *numpy.sqrt(tdelta*f/R))) # //in mm
    exfb = numpy.sqrt(R*F/tdelta) * 1/(numpy.sin((N+1)*numpy.sqrt(tdelta*f/R))) # //in mm

    exf = 0
    if (numpy.abs(exfa-f) > numpy.abs(exfb-f) and (N+1)<=numpy.sqrt(R/(tdelta*F))*numpy.pi/2):
        exf = exfb
        N = N + 1
    else:
        exf = exfa

    if exf >= dist:
        raise Exception("The source lies within the focal distance.")

    print("exf: ", exf)

        # for N+1 <= N_min
    # 	exf=exfb;
    # 	N=N+1;
    # } else {
    # 	exf=exfa;
    # }
    # if (exf>=dist) {
    # 	alert("The source lies within the focal distance.");
    # 	return false;
    # }


# for(var i=0; i<l && muoverrhos[i][0]<=e; i++) {}
# if (i==0) {
# 	mu=rho*muoverrhos[0][1];
# } else if (i==l) {
# 	mu=rho*muoverrhos[l-1][1];
# } else {
# mu=rho*(muoverrhos[i][1]*(e-muoverrhos[i-1][0])/(muoverrhos[i][0]-muoverrhos[i-1][0])-muoverrhos[i-1][1]*(e-muoverrhos[i][0])/(muoverrhos[i][0]-muoverrhos[i-1][0]));
# }
#
# //compute desired focal length and real number of lenses
# dist=dist*1e+3;//dist in mm
# if (document.getElementById("presf").checked) {//compute focal length
# 	f=desf;
# } else {
# 	f=dist*desf/(dist+desf);
# }
# var N=Math.sqrt(R/(tdelta*F))*Math.asin(Math.sqrt(R*F/tdelta)/f);
# if (isNaN(N)) {
# 	alert("The focus lies within the CRL.");
# 	return false;
# }
#
# //compute integer number of lenses and true focal length
# N=Math.floor(N);
# exfa=Math.sqrt(R*F/tdelta)*1/(Math.sin(N*Math.sqrt(tdelta*F/R)));//in mm
# exfb=Math.sqrt(R*F/tdelta)*1/(Math.sin((N+1)*Math.sqrt(tdelta*F/R)));//in mm
# var exf;
# if (Math.abs(exfa-f)>Math.abs(exfb-f) && (N+1)<=Math.sqrt(R/(tdelta*F))*Math.PI/2) {//test for: N+1 <= N_min
# 	exf=exfb;
# 	N=N+1;
# } else {
# 	exf=exfa;
# }
# if (exf>=dist) {
# 	alert("The source lies within the focal distance.");
# 	return false;
# }

#
# """
# function compute(){
# 	var e=parseFloat(document.getElementById("e").value);
# 	var desf=parseFloat(document.getElementById("desf").value);
# 	var dist=parseFloat(document.getElementById("dist").value);
# 	var R=parseFloat(document.getElementById("R").value);
# 	var w=parseFloat(document.getElementById("w").value);
# 	var h=parseFloat(document.getElementById("h").value);
# 	var F=parseFloat(document.getElementById("F").value);
# 	var d=parseFloat(document.getElementById("d").value);
# 	var th=parseFloat(document.getElementById("th").value);
# 	if (isNaN(e) || isNaN(desf) || isNaN(dist) || isNaN(w) || isNaN(h) || isNaN(F) || isNaN(d) || isNaN(th)) {
# 		alert("Please enter admissible values.");
# 		return false;
# 	}
#
# 	//compute mu
# 	var tdelta,mu,lambda,rho,muoverrhos;
# 	lambda=12.398/e;
# 	if (document.getElementById("Be").checked) {
# 		rho=1.85;
# 		tdelta=5.4*1e-6*Math.pow(lambda,2)*rho*4/9.012;
# 		muoverrhos=[[1.0,604.1],[1.5,179.7],[2.0,74.69],[3.0,21.27],[4.0,8.685],[5.0,4.369],[6.0,2.527],[8.0,1.124],[10.0,0.6466],[15.0,0.307],[20.0,0.2251],[30.0,0.1792],[40.0,0.164],[50.0,0.1554],[60.0,0.1493],[80.0,0.1401],[100.0,0.1328],[150.0,0.119],[200.0,0.1089],[300.0,0.09463],[400.0,0.08471],[500.0,0.07739],[600.0,0.07155],[800.0,0.06286],[1000.0,0.05652],[1250.0,0.05054],[1500.0,0.04597],[2000.0,0.03938],[3000.0,0.03138],[4000.0,0.02664],[5000.0,0.02347],[6000.0,0.02121],[8000.0,0.01819],[10000.0,0.01627],[15000.0,0.01361],[20000.0,0.01227]];
# 	} else if (document.getElementById("Al").checked) {
# 		rho=2.7;
# 		tdelta=5.4*1e-6*Math.pow(lambda,2)*rho*13/26.982;
# 		muoverrhos=[[1.0,1185.0],[1.5,402.2],[1.5596,362.1],[1.5596,3957.0],[2.0,2263.0],[3.0,788.0],[4.0,360.5],[5.0,193.4],[6.0,115.3],[8.0,50.33],[10.0,26.23],[15.0,7.955],[20.0,3.441],[30.0,1.128],[40.0,0.5685],[50.0,0.3681],[60.0,0.2778],[80.0,0.2018],[100.0,0.1704],[150.0,0.1378],[200.0,0.1223],[300.0,0.1042],[400.0,0.09276],[500.0,0.08445],[600.0,0.07802],[800.0,0.06841],[1000.0,0.06146],[1250.0,0.05496],[1500.0,0.05006],[2000.0,0.04324],[3000.0,0.03541],[4000.0,0.03106],[5000.0,0.02836],[6000.0,0.02655],[8000.0,0.02437],[10000.0,0.02318],[15000.0,0.02195],[20000.0,0.02168]];
# 	} else {
# 		rho=8.90;
# 		tdelta=5.4*1e-6*Math.pow(lambda,2)*rho*28/58.693;
# 		muoverrhos=[[1.0,9855.0],[1.00404,9753.0],[1.0081,9654.0],[1.0081,10990.0],[1.5,4234.0],[2.0,2049.0],[3.0,709.4],[4.0,328.2],[5.0,179.3],[6.0,109.0],[8.0,49.52],[8.3328,44.28],[8.3328,329.4],[10.0,209.0],[15.0,70.81],[20.0,32.2],[30.0,10.34],[40.0,4.6],[50.0,2.474],[60.0,1.512],[80.0,0.7306],[100.0,0.444],[150.0,0.2208],[200.0,0.1582],[300.0,0.1154],[400.0,0.09765],[500.0,0.08698],[600.0,0.07944],[800.0,0.06891],[1000.0,0.0616],[1250.0,0.05494],[1500.0,0.05015],[2000.0,0.04387],[3000.0,0.03745],[4000.0,0.03444],[5000.0,0.03289],[6000.0,0.0321],[8000.0,0.03164],[10000.0,0.03185],[15000.0,0.0332],[20000.0,0.03476]];
# 	}
# 	var l=muoverrhos.length;
# 	for(var i=0; i<l && muoverrhos[i][0]<=e; i++) {}
# 	if (i==0) {
# 		mu=rho*muoverrhos[0][1];
# 	} else if (i==l) {
# 		mu=rho*muoverrhos[l-1][1];
# 	} else {
# 	mu=rho*(muoverrhos[i][1]*(e-muoverrhos[i-1][0])/(muoverrhos[i][0]-muoverrhos[i-1][0])-muoverrhos[i-1][1]*(e-muoverrhos[i][0])/(muoverrhos[i][0]-muoverrhos[i-1][0]));
# 	}
#
# 	//compute desired focal length and real number of lenses
# 	dist=dist*1e+3;//dist in mm
# 	if (document.getElementById("presf").checked) {//compute focal length
# 		f=desf;
# 	} else {
# 		f=dist*desf/(dist+desf);
# 	}
# 	var N=Math.sqrt(R/(tdelta*F))*Math.asin(Math.sqrt(R*F/tdelta)/f);
# 	if (isNaN(N)) {
# 		alert("The focus lies within the CRL.");
# 		return false;
# 	}
#
# 	//compute integer number of lenses and true focal length
# 	N=Math.floor(N);
# 	exfa=Math.sqrt(R*F/tdelta)*1/(Math.sin(N*Math.sqrt(tdelta*F/R)));//in mm
# 	exfb=Math.sqrt(R*F/tdelta)*1/(Math.sin((N+1)*Math.sqrt(tdelta*F/R)));//in mm
# 	var exf;
# 	if (Math.abs(exfa-f)>Math.abs(exfb-f) && (N+1)<=Math.sqrt(R/(tdelta*F))*Math.PI/2) {//test for: N+1 <= N_min
# 		exf=exfb;
# 		N=N+1;
# 	} else {
# 		exf=exfa;
# 	}
# 	if (exf>=dist) {
# 		alert("The source lies within the focal distance.");
# 		return false;
# 	}
#
# 	var lcrl=N*F;//length of crl
# 	var ol=N*Math.sqrt(F*tdelta/R);//\omega*L
# 	var feff=exf*Math.cos(ol);//effective focal length in mm
# 	var pp=exf-feff-lcrl/2;//position of principal plane in mm
# 	var idist=exf*dist/(dist-exf);//image distance in mm
# 	var Neff=(N/2+0.25*Math.sin(2*ol)/Math.sin(ol/N));//effective number of lenses
# 	var gapp=2*Math.sqrt(R*1e+3*(th*1e+3-d));//geometric aperture in um
#
# 	//compute effective aperture
# 	var ap=0.1*mu*Neff*Math.pow(gapp/2e+3,2)/(2*R);
# 	var Deff;
# 	var twod=document.getElementById("D2").checked;
# 	if (twod) {
# 		Deff=gapp*Math.sqrt((1-Math.exp(-ap))/ap);//effective aperture in um for 2D lenses
# 	} else {
# 		Deff=gapp*0.5*Math.sqrt(Math.PI/ap)*erf(Math.sqrt(ap));//effective aperture in um for 1D lenses
# 	}
#
# 	//compute resolution and focal spot size
# 	var corr=2*Math.sqrt(2*Math.log(2));
# 	var res=corr*lambda*1e+2*exf/(Math.PI*Deff);//resolution limit in nm
# 	var fspotw=corr*idist*1e+6*Math.sqrt(Math.pow(w/(corr*dist*1e+3),2)+Math.pow(lambda*1e-4/(Math.PI*Deff),2));//focal spot size in nm
# 	if (twod) {
# 		var fspoth=corr*idist*1e+6*Math.sqrt(Math.pow(h/(corr*dist*1e+3),2)+Math.pow(lambda*1e-4/(Math.PI*Deff),2));//focal spot size (height) in nm
# 	}
#
# 	//compute transmission and gain
# 	var aD=0.1*mu*Neff*Math.pow(Deff/2e+3,2)/(2*R);
# 	var transm,transmgapp;
# 	if (twod) {
# 		transm=Math.exp(-mu*N*d*1e-4)*(1-Math.exp(-2*aD))/(2*aD);//transmission for 2D lenses (ill. of eff. ap.)
# 		transmgapp=Math.exp(-mu*N*d*1e-4)*(1-Math.exp(-2*ap))/(2*ap);//transmission for 2D lenses (ill. of geom. ap.)
# 	} else {
# 		transm=Math.exp(-mu*N*d*1e-4)*Math.sqrt(Math.PI)*0.5*erf(Math.sqrt(2*aD))/Math.sqrt(2*aD);//transmission for 1D lenses (ill. of eff. ap.)
# 		transmgapp=Math.exp(-mu*N*d*1e-4)*Math.sqrt(Math.PI)*0.5*erf(Math.sqrt(2*ap))/Math.sqrt(2*ap);//transmission for 1D lenses (ill. of geom. ap.)
# 	}
# 	var gain;
# 	if (twod) {
# 		gain=transmgapp*Math.pow(1e+3*gapp,2)/(fspotw*fspoth);//gain for 2D lenses
# 	} else {
# 		gain=transmgapp*1e+3*gapp/fspotw;//gain for 1D lenses
# 	}
#
# 	document.getElementById("delta").innerHTML=(1e+6*tdelta/2).toFixed(3)+" &times; 10<sup>-6</sup>";
# 	document.getElementById("mu").innerHTML=mu.toFixed(3)+" cm<sup>-1</sup>";
# 	document.getElementById("N").innerHTML=N;
# 	document.getElementById("exf").innerHTML=exf.toFixed(1)+" mm";
# 	document.getElementById("gapp").innerHTML=Math.round(gapp)+" &mu;m";
# 	document.getElementById("Deff").innerHTML=Math.round(Deff)+" &mu;m";
# 	document.getElementById("eidist").innerHTML=(idist-pp-N*F/2).toFixed(1)+" mm";
# 	document.getElementById("lcrl").innerHTML=lcrl.toFixed(1)+" mm";
# 	document.getElementById("pp").innerHTML=pp.toFixed(1)+" mm";
# 	document.getElementById("res").innerHTML=Math.round(res)+" nm";
# 	document.getElementById("idist").innerHTML=idist.toFixed(1)+" mm";
# 	if (twod) {
# 		document.getElementById("fspot").innerHTML=Math.round(fspotw)+" nm &times; "+Math.round(fspoth)+" nm";
# 	} else {
# 		document.getElementById("fspot").innerHTML=Math.round(fspotw)+" nm";
# 	}
# 	document.getElementById("transm").innerHTML=transmgapp.toFixed(3);
# 	document.getElementById("gain").innerHTML=Math.round(gain);
# }
#
# function erf(x) {//A&S formula 7.1.26
# 	var a1=0.254829592;
# 	var a2=-0.284496736;
# 	var a3=1.421413741;
# 	var a4=-1.453152027;
# 	var a5=1.061405429;
# 	var p=0.3275911;
# 	var t=1.0/(1.0+p*x);
# 	return y=1.0-(((((a5*t+a4)*t)+a3)*t+a2)*t+a1)*t*Math.exp(-x*x);
# }
#
# function choose() {
# 	if (document.getElementById("presf").checked) {
# 		document.getElementById("choice").innerHTML="Desired focal length:";
# 	} else {
# 		document.getElementById("choice").innerHTML="Desired image distance:";
# 	}
#
# 	if (document.getElementById("D2").checked) {
# 		document.getElementById("htext").style.visibility="visible";
# 		document.getElementById("hinput").style.visibility="visible";
# 		document.getElementById("ttwidth").innerHTML="Please specify the width (FWHM) of the x-ray source.";
# 		document.getElementById("ttfspot").innerHTML="The lateral extent (FWHM, width &times; height) of the secondary source (geometric plus diffraction contribution).";
# 	} else {
# 		document.getElementById("htext").style.visibility="hidden";
# 		document.getElementById("hinput").style.visibility="hidden";
# 		document.getElementById("ttwidth").innerHTML="Please specify the lateral extent (FWHM) of the x-ray source in the focused direction.";
# 		document.getElementById("ttfspot").innerHTML="The lateral extent (FWHM) of the secondary source (geometric plus diffraction contribution) in the focused direction.";
# 	}
# }
#
# document.getElementById("nojs").innerHTML="";
# document.getElementById("compute").addEventListener("click",compute);
# document.getElementById("presf").addEventListener("click",choose);
# document.getElementById("presi").addEventListener("click",choose);
# document.getElementById("D2").addEventListener("click",choose);
# document.getElementById("D1").addEventListener("click",choose);
#
# """


if __name__ == "__main__":

    compute()