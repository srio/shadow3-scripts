#!/usr/bin/env python
# encoding: utf-8
"""
@file test_zern.py
@brief Test libtim/zern.py library

@author Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
@copyright Creative Commons Attribution-Share Alike license versions 3.0 or higher, see http://creativecommons.org/licenses/by-sa/3.0/
@date 20120927

Testcases for zern.py library.
"""

from zern import *
import unittest
# import libtim.im
# import pyfits
from timeit import Timer
import os
from functools import reduce
# from im import inter_imshow


SHOWPLOTS=True


if SHOWPLOTS:
	import pylab as plt
	from matplotlib import cm

#
# copied from im.py (srio)
#

def raw_input(x):
	input(x)

def inter_imshow(data, desc="", doshow=True, dowait=False, log=False, rollaxes=False, cmap='RdYlBu', figid=None, **kwargs):
	"""Show data using matplotlib.imshow if **doshow** is true.

	Additionally, print **desc** just before plotting so users know what they see. If **dowait** is True (default), wait for input before continuing.

	This function is used to show intermediate results of analysis programs conditionally, i.e. only when a certain flag is set.

	@param [in] data 2D data to plot
	@param [in] desc Text to print before plot and plot title
	@param [in] doshow Show only if this is True
	@param [in] dowait If set, wait before continuing
	@param [in] log Take logarithm of data before plotting
	@param [in] rollaxes Roll axes for plot such that (0,0) is the center
	@param [in] cmap Colormap to use (RdYlBu or YlOrBr are nice)
	@param [in] figid Figure id to use in plt.figure(figid), can be used to re-use plot windows.
	@param [in] **kwargs Additional arguments passed to imshow()
	"""

	if (not doshow):
		return

	import pylab as plt
	print ("inter_imshow(): " + desc)

	# Pre-format data
	data_arr = np.asanyarray(data)
	if (log):
		data_arr = np.log10(data_arr)

	# Check if we want to roll the axis
	extent = None
	if (rollaxes):
		sh = data_arr.shape
		extent = (-sh[1]/2., sh[1]/2., -sh[0]/2., sh[0]/2.)

	fig = plt.figure(figid)
	fig.clf()
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
	ax.set_title(desc)

	img = ax.imshow(data_arr, extent=extent, cmap=cm.get_cmap(cmap), **kwargs)
	fig.colorbar(img, aspect=30, pad=0.05)
	plt.draw()

	plt.show()

	# If we want to wait, ask user for input, discard it and continue
	if (dowait):
		raw_input("<?>")

#
#
#
class PlotZernikes(unittest.TestCase):
	def test0a_masking(self):
		"""Plot masking"""
		if (not SHOWPLOTS):
			return
		rad = 128
		zn_m = calc_zernike([1,2,3,4,5,6], rad, mask=True)
		zn_unm = calc_zernike([1,2,3,4,5,6], rad, mask=False)
		mask = mk_rad_mask(2*rad) <= 1
		inter_imshow(zn_m, desc="Masked Zernike")
		inter_imshow(zn_unm, desc="Unmasked Zernike")
	def test1_plot_basis(self):
		"""Plot the first few Zernike basis functions, until the user quits"""
		if (not SHOWPLOTS):
			return
		rad = 257
		# for idx in range(100):
		for idx in range(5):
			vec = [0]*(idx+1)
			vec[idx] = 1
			mode = calc_zernike(vec, rad, mask=True)
			nn, mm = noll_to_zern(idx+1)
			# inter_imshow(mode, desc="Zernike mode j=%d (n=%d, m=%d) " % (n+1, nn, mm))
			inter_imshow(mode, desc="Zernike mode j=%d (n=%d, m=%d) " % (idx+1, nn, mm))
class StoreZernikes(unittest.TestCase):
	def setUp(self):
		"""Generate source image, darkfield, flatfield and simulated data"""
		self.outdir = './test_zern_out/'
		try:
			os.makedirs(self.outdir)
		except OSError:
			pass
		rad = 257
		self.rad = rad
		self.nmodes = 30
		self.basis = []
		self.basis_data = calc_zern_basis(self.nmodes, self.rad)
		self.basis = self.basis_data['modes']*self.basis_data['mask']

	# def test1_file_write(self):
	# 	"""Test disk writing & plotting (might take a while...)"""
	# 	if (not SHOWPLOTS): return
	# 	# Store as one big FITS file
	# 	outf = os.path.join(self.outdir, 'zernike-basis.fits')
	# 	pyfits.writeto(outf, np.r_[self.basis], clobber=True)
	#
	# 	# Plot each mode
	# 	for n, basis in enumerate(self.basis):
	# 		nn, mm = noll_to_zern(n+1)
	# 		libtim.im.store_2ddata(basis, 'zernike-basis_%02d' % (n+1),
	# 			dir=self.outdir, pltitle='Zernike mode j=%d, (n=%d, m=%d) ' % (n+1, nn, mm), xlab='X', ylab='Y', ident=False)

class TestZernikes(unittest.TestCase):
	def setUp(self):
		"""Generate source image, darkfield, flatfield and simulated data"""
		rad = 257
		rad = 127
		self.rad = rad
		self.nmodes = 25
		self.vec = np.random.random(self.nmodes)
		self.basis = []
		self.basis_data = calc_zern_basis(self.nmodes, self.rad)
		self.basis = self.basis_data['modes']
		self.wf = reduce(lambda x,y: x+y[1]*self.basis[y[0]], enumerate(self.vec), 0)
		self.wf_msk = reduce(lambda x,y: x+y[1]*self.basis[y[0]]*self.basis_data['mask'], enumerate(self.vec), 0)

	# Shallow data tests
	def test0a_basis_data(self):
		"""Length of basis should be the same as amplitude vector length"""
		self.assertEqual(len(self.basis), len(self.vec))

	def test0b_data(self):
		"""Shape of basis and test wf should be equal"""
		for mode in self.basis:
			self.assertEqual(self.wf.shape, mode.shape)

	def test0c_indices(self):
		"""Test Noll to Zernike index conversion, see
		<https://oeis.org/A176988>"""
		zern = [(0,0), (1,1), (1,-1), (2,0), (2,-2), (2,2), (3,-1), (3,1), (3,-3), (3,3)]
		for j, nm in enumerate(zern):
			self.assertEqual(noll_to_zern(j+1), nm)

	# Shallow function test
	def test1a_zero_wf(self):
		"""Zero-wavefront should return zero amplitude vector"""
		fitdata = fit_zernike(np.zeros((64, 64)), nmodes=self.nmodes)
		fitvec = fitdata[0]
		self.assertTrue(np.allclose(fitvec, 0.0))

	def test1b_shape_detection(self):
		"""Test calc_zern_basis basis shape detection"""
		basis_data = calc_zern_basis(nmodes=10, rad=64)
		basis = basis_data['modes']
		self.assertEqual(len(basis), 10)
		self.assertEqual(basis[0].shape, (128, 128))
		self.assertEqual(basis[0].shape, basis[-1].shape)

		basis_data = calc_zern_basis(nmodes=5, rad=128)
		basis = basis_data['modes']
		self.assertEqual(len(basis), 5)
		self.assertEqual(basis[0].shape, (256, 256))
		self.assertEqual(basis[0].shape, basis[-1].shape)

		basis_data = calc_zern_basis(nmodes=5, rad=64)
		basis = basis_data['modes']
		self.assertEqual(len(basis), 5)
		self.assertEqual(basis[0].shape, (128, 128))
		self.assertEqual(basis[0].shape, basis[-1].shape)

	def test1d_masking(self):
		"""Test if masking works"""
		rad = 128
		zn_m = calc_zernike([1,2,3,4,5,6], rad, mask=True)
		zn_unm = calc_zernike([1,2,3,4,5,6], rad, mask=False)
		mask = mk_rad_mask(2*rad) <= 1

		# Should be equal inside mask
		self.assertAlmostEqual(np.sum(zn_unm[mask] - zn_m[mask]), 0.0)
		self.assertTrue(np.allclose(zn_unm[mask], zn_m[mask]))

		# Should be unequal outside mask
		self.assertFalse(np.allclose(zn_unm[mask==False], zn_m[mask==False]))

		# Mean outside the mask should be larger than the mean inside the mask, in general
		self.assertGreater(np.mean(np.abs(zn_unm[mask==False])),
			0.5*np.mean(np.abs(zn_unm[mask])))

	# Deep function tests
	# calc_zern_basis(nmodes, rad, mask=True):
	# fit_zernike(wavefront, nmodes=10, center=(-0.5, -0.5), rad=-0.5):
	# calc_zernike(zern_vec, rad=-1):

	def test2b_zern_calc(self):
		"""Compare calc_zernike output with pre-computed basis"""
		vec = [0]*self.nmodes
		for i in range(self.nmodes):
			vec[i] = 1
			testzern = calc_zernike(vec, self.rad, mask=False)
			self.assertTrue(np.allclose(self.basis[i], testzern))
			vec[i] = 0

	def test2c_variance(self):
		"""Test whether all Zernike modes have variance unity"""
		self.mask = mk_rad_mask(2*self.rad) <= 1

		for idx, m in enumerate(self.basis):
			if (idx == 0):
				continue
			# The error in the variance should scale with the number of
			# pixels, more pixels means less error because of better sampling.
			# Because of this we take 1/npixels as an error margin. The factor
			# 1.1 is added for numerical roundoff and other computer errors.
			# We use the mask to test only the relevant part of the Zernike
			# modes
			self.assertAlmostEqual(np.var(m[self.mask]), 1.0, delta=1.1/(self.rad**2.)**0.5)

	def test2d_equal_mode(self):
		"""Test equal-mode Zernike reconstruction with non-masked input"""
		fitvec, fitrec, fitdiff = fit_zernike(self.wf, nmodes=self.nmodes)

		self.assertAlmostEqual(np.sum(self.vec - fitvec), 0.0)
		self.assertTrue(np.allclose(self.vec, fitvec))

	def test2d2_equal_mode_masked(self):
		"""Test equal-mode Zernike reconstruction with masked input"""
		fitvec, fitrec, fitdiff = fit_zernike(self.wf_msk, nmodes=self.nmodes)

		self.assertAlmostEqual(np.sum(self.vec - fitvec), 0.0)
		self.assertTrue(np.allclose(self.vec, fitvec))

	def test2e_unequal_mode(self):
		"""Test unequal-mode Zernike reconstruction with non-masked input. This will probably not be exactly equal because we try to fit a surface generated with 20 modes with only 10 modes."""
		fitvec, fitrec, fitdiff = fit_zernike(self.wf, nmodes=10)

		self.assertAlmostEqual(np.mean(self.vec[:10] / fitvec), 1.0, delta=0.1)

	def test2e2_unequal_mode(self):
		"""Test unequal-mode Zernike reconstruction with masked input. This will probably not be exactly equal because we try to fit a surface generated with 20 modes with only 10 modes."""
		fitvec, fitrec, fitdiff = fit_zernike(self.wf_msk, nmodes=10)

		self.assertAlmostEqual(np.mean(self.vec[:10] / fitvec), 1.0, delta=0.1)

	def test2f_fit_startmode(self):
		"""Test startmode parameter in fit_zernike"""
		# startmode == 0 should raise an error, as this is not a valid Noll
		# index
		with self.assertRaises(ValueError):
			fitdata = fit_zernike(self.wf, nmodes=10, startmode=0)

		# Setting startmode higher should block out the first few modes
		for s in range(10):
			fitdata = fit_zernike(self.wf, nmodes=10, startmode=s+1)
			fitvec = fitdata[0]
			self.assertEqual(tuple(fitvec[:s]), (0,)*s)

	def test2g_weighed_fit(self):
		"""Test weighed Zernike reconstruction"""
		this_wf = self.wf_msk.copy()
		this_weight = self.wf_msk.copy()*0+1

		# Set region in center to random values, and set weight to 0
		this_wf[int(this_wf.shape[0]*2/3):] = 2*this_wf.max()
		this_weight[int(this_wf.shape[0]*2/3):] = 0

		# Now fit and plot
		fitvecw, fitrecw, fitdiffw = fit_zernike(this_wf, fitweight=this_weight, nmodes=self.nmodes)
		fitvec, fitrec, fitdiff = fit_zernike(this_wf, nmodes=self.nmodes)

		# When fitting the corrupted wavefront, the residual with the original
		# wavefront should be much bigger than with the corrupted wavefront
		fitw_orig = np.abs(fitrecw - self.wf_msk).mean()
		fit_orig = np.abs(fitrec - self.wf_msk).mean()
		fitw_corr = np.abs(fitrecw - this_wf).mean()
		fit_corr = np.abs(fitrec - this_wf).mean()

		self.assertGreater(fit_orig/1.2, fit_corr)
		self.assertGreater(fitw_corr/10., fitw_orig)
		self.assertGreater(fitrec.max()/1.5, fitrecw.max())

		if (SHOWPLOTS):
			inter_imshow(self.wf_msk, desc="Input data")
			inter_imshow(this_wf, desc="Corrupted input data")
			inter_imshow(this_weight, desc="Input weight")
			inter_imshow(fitrecw, desc="Weighed Reconstruction")
			inter_imshow(fitrec, desc="Normal Reconstruction")

class TestZernikeSpeed(unittest.TestCase):
	def setUp(self):
		self.calc_iter = 3
		self.fit_iter = 5
		self.nmodes = 25
		self.rad = 257

	def test3a_timing_calc(self):
		"""Test Zernike calculation timing and cache functioning"""

		t1 = Timer("""
a=calc_zernike(vec, rad, z_cache)
		""", """
from zern import calc_zern_basis, fit_zernike, calc_zernike
import numpy as np
rad = %d
nmodes = %d
vec = np.random.random(nmodes)
z_cache = {}
		""" % (self.rad, self.nmodes) )
		t2 = Timer("""
a=calc_zernike(vec, rad, {})
		""", """
from zern import calc_zern_basis, fit_zernike, calc_zernike
import numpy as np
rad = %d
nmodes = %d
vec = np.random.random(nmodes)
		""" % (self.rad, self.nmodes) )
		t_cache = t1.timeit(self.calc_iter)/self.calc_iter
		t_nocache = t2.timeit(self.calc_iter)/self.calc_iter
		# Caching should be at least twice as fast as no caching
		# Note that here we do not initialize the cache in the setup, it is
		# set to an empty dict which is filled on first run. This test should
		# test that this automatic filling works properly
		self.assertGreater(t_nocache/2.0, t_cache)

	def test3b_timing_calc(self):
		"""Test Zernike calculation performance with and without cache, print results"""

		t1 = Timer("""
a=calc_zernike(vec, rad, z_cache)
		""", """
from zern import calc_zern_basis, fit_zernike, calc_zernike
import numpy as np
rad = %d
nmodes = %d
vec = np.random.random(nmodes)
z_cache = calc_zern_basis(len(vec), rad)
		""" % (self.rad, self.nmodes) )

		t2 = Timer("""
a=calc_zernike(vec, rad, {})
		""", """
from zern import calc_zern_basis, fit_zernike, calc_zernike
import numpy as np
rad = %d
nmodes = %d
vec = np.random.random(nmodes)
		""" % (self.rad, self.nmodes) )

		t_cached = min(t1.repeat(2, self.calc_iter))/self.calc_iter
		t_nocache = min(t2.repeat(2, self.calc_iter))/self.calc_iter
		print ("test3b_timing_calc(): rad=257, nmodes=25 cache: %.3g s/it no cache: %.3g s/it" % (t_cached, t_nocache))

	def test3c_timing_fit(self):
		"""Test Zernike fitting performance"""

		t1 = Timer("""
a=fit_zernike(wf, z_cache, nmodes=nmodes)
		""", """
from zern import calc_zern_basis, fit_zernike, calc_zernike
import numpy as np
rad = %d
nmodes = %d
vec = np.random.random(nmodes)
z_cache = calc_zern_basis(len(vec), rad)
wf = np.random.random((rad, rad))
		""" % (self.rad, self.nmodes) )

		t_cached = min(t1.repeat(2, self.fit_iter))/self.fit_iter
		# Caching should be at least twice as fast as no caching
		print ("test3c_timing_fit(): rad=257, nmodes=25 %.3g sec/it" % (t_cached))

if __name__ == "__main__":
	import sys
	sys.exit(unittest.main())
