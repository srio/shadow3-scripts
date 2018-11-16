rx=1.998731
drx=0.546210
ry=1.998731
dry=0.546210
rescale_x=1.000000
rescale_y=1.000000
s_id="0_hib3-3302"
i_mode=49

import pickle
from srwlib import *
from comsyl.waveoptics.SRWAdapter import SRWAdapter
from comsyl.waveoptics.Wavefront import NumpyWavefront, SRWWavefront

wavefront = NumpyWavefront.load("./tmp%s_in.npz"%s_id)
adapter = SRWAdapter()
wfr = adapter.SRWWavefrontFromWavefront(wavefront,
                                        rx,
                                        drx,
                                        ry,
                                        dry,rescale_x,rescale_y)
#print("Doing propagation in external call")
srw_beamline = pickle.load(open("./tmp%s_beamline.p"%s_id,"rb"))
srwl.PropagElecField(wfr, srw_beamline)

tmp = SRWWavefront(wfr).toNumpyWavefront()
wfr = None
tmp.save("./tmp%s_out" % s_id)
