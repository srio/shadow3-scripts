import at
import math
import numpy as np

lattice_file = '/operation/beamdyn/optics/sr/theory/betamodel.mat'
r = at.load_lattice(lattice_file)

def get_parameters_at_s_locations(r0=r, s_locs=None, filename='ebs_optics.spec', verbose=False):
    r"""Computes optics geometry and beam sigmas at given s locations

    See :py:func:`get_parameters_at_s_locations`

    Parameters:
        r0:         pyAT Lattice 
        s_locs:     Minimum distance between 2 grid points
        filename:   file where to save data, default: ebs_optics.spec
        verbose:    print some text

    Returns:
        saves data in spec format into the file filename
        data:   (len(s_locs), 27) array: L s-Circumference[m]  s[m]  x  y  angle  alphaX  alphaY  betaX  betaY  gammaX  gammaY  etaX  etaY  etaXp  etaYp  <xx>  <yy>  <xxp>  <yyp>  <xpxp>  <ypyp>  SigX[um]  SigY[um]  SigXp[urad]  SigYp[urad]  epsX  epsY

    """

    # get circumference
    Circumference = r0.get_s_pos(at.End)

    # define break locations
    if type(s_locs) != np.ndarray:
        npoints = 1 + int(Circumference[0]*1)
        s_locs = np.linspace(0.0, Circumference, npoints)[:, 0] # all locations along lattice
        
    npoints = s_locs.shape[0]
    
    # insert markers at break locations
    r = r0.sbreak(break_s=list(s_locs))

    # indexes of s locations
    s_ind = r.get_uint32_index('sbreak')

    # get lattice parameters with radiation
    if verbose:
        print('get lattice parameters')
    r.enable_6d()
    p0 = r.envelope_parameters()

    epsilonX = p0.emittances[0]; 
    epsilonY = 10*1e-12; # tuned to this value during operation
    delta = p0.sigma_e; 

    # gert optics
    if verbose:
        print('get orbit, dispersion, beta functions')
    _, _, l = r.linopt6(refpts=s_ind)
    
    # get geometry
    if verbose:
        print('get geomtery')
    geom, _ = r.get_geometry(refpts=s_ind)
    
    data = []

    # write to file
    with open(filename,'w') as f:
        f.write(f'#F {filename}\n')
        f.write(f'#ULATTICEFILE {lattice_file} \n')
        f.write(f'#UEPSILONX {epsilonX:15g} \n')
        f.write(f'#UEPSILONY {epsilonY:15g} \n')
        f.write(f'#UENERGYSPREAD {delta:15g} \n')
        f.write('\n#S 1 data from pyAT using script get_ebs_optics.py \n')
        f.write('#N 27\n')
        f.write('#L s-sW[m]  s[m]  x  y  angle  alphaX  alphaY  betaX  betaY  gammaX  gammaY  etaX  etaY  etaXp  etaYp  <xx>  <yy>  <xxp>  <yyp>  <xpxp>  <ypyp>  SigX[um]  SigY[um]  SigXp[urad]  SigYp[urad]  epsX  epsY\n')
        
        for i in range(0, npoints):
            s = s_locs[i]
            s0 = s-Circumference

            alpha = l[i].alpha
            alphaX = alpha[0]
            alphaY = alpha[1]

            beta = l[i].beta
            betaX = beta[0] 
            betaY = beta[1]

            gammaX = (1.0 + alpha[0]*alpha[0])/beta[0]
            gammaY = (1.0 + alpha[1]*alpha[1])/beta[1]

            eta = l[i].dispersion
            etaX = eta[0]
            etaXp = eta[1]
            etaY = eta[2]
            etaYp = eta[3]

            xx = betaX*epsilonX + (etaX * delta)**2
            yy = betaY*epsilonY + (etaY * delta)**2
            xxp = -alphaX * epsilonX + etaX * etaXp * delta**2
            yyp = -alphaY * epsilonY + etaY * etaYp * delta**2
            xpxp = gammaX * epsilonX + (etaXp * delta)**2
            ypyp = gammaY * epsilonY + (etaYp * delta)**2

            lab_x = geom[i].x
            lab_y = geom[i].y
            angle = geom[i].angle

            tmp = [s0[0],
                   s,
                   lab_x,
                   lab_y,
                   angle,
                   alphaX,
                   alphaY,
                   betaX,
                   betaY,
                   gammaX,gammaY,
                   etaX,etaY,etaXp,etaYp,
                   xx,yy,xxp,yyp,xpxp,ypyp,
                   1e6*math.sqrt(xx),1e6*math.sqrt(yy),1e6*math.sqrt(xpxp),
                   1e6*math.sqrt(ypyp),math.sqrt(xx*xpxp),math.sqrt(yy*ypyp)]
            
            data.append(tmp)

            tmp_strings = [f'{t:12.6f}' for t in tmp[0:-13]] + [f'{t:2.12g}' for t in tmp[-12::]]

            f.write('\t'.join(tmp_strings)+'\n')

    return np.array(data) 

if __name__=='__main__':

    # get parameters at some s locations
    data = get_parameters_at_s_locations(s_locs=np.array([0.0, 10.3, 26.374, 567.84]))
    print(data)
    
    # DISPLAY PARAMETERS AT center of Straigth section
    id = 'ID08'
    filename=f'{id}_params.txt'
    IDind = r.get_uint32_index(id)
    s = r.get_s_pos(IDind)
    get_parameters_at_s_locations(s_locs=s, filename=filename)
    with open(filename,'r') as f:
        print(f.read())
    