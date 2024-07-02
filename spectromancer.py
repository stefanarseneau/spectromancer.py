import argparse
import numpy as np
from tqdm import tqdm
import glob
import os

from astropy.io import fits
from astropy.table import Table
import corv

class Spectrum:
    def __init__(self, wave, flux, ivar, mask = None):
        self.wave = wave
        self.flux = flux
        self.ivar = ivar
        self.mask = mask

    def fit_rv(self, model):
        return corv.fit.fit_corv(self.wave, self.flux, self.ivar, model)
    
    def plot(self, model, params):
        corv.utils.lineplot(self.wave, self.flux, self.ivar, model, params)

    def spectrum_coadd(self, other):
        # interpolate the other spectrum onto a consistent wavelength grid
        flux_rebase = np.interp(self.wave, other.flux, other.ivar)
        ivar_rebase = np.interp(self.wave, other.wave, other.ivar)

        # create holder objects for fluxes and ivars
        fls = np.array([self.flux, flux_rebase])
        ivars = np.array([self.ivar, ivar_rebase])

        # perform the coadd
        mask = (ivars == 0) # create a boolean mask of ivar = 0
        ivars[mask] = 1 # set those points to 1 to avoid divide by zero
        variances = ivars**-1 # compute variances of both spectra
        variance = np.sum(variances, axis = 0) / len(variances)**2 # central limit theorem
        ivar = variance**-1 # put this back into ivar
        flux = np.median(fls, axis=0) # compute the median of each point
        smask = (mask).all(axis=0) # reduce the dimensionality of the mask
        ivar[smask] = 1e-6 # reset the points to something near zero
        
        # return a new object that is the co-added spectra
        return Spectrum(self.wave, flux, ivar, mask)

    # the + operation performs a coadd of two spectra
    def __add__(self, other):
        return self.spectrum_coadd(other)

    def __radd__(self, other):
        return self + other

    
class Target:
    def __init__(self, name, path):
        self.name = name # save the name of the target
        self.exp_paths = glob.glob(path + '/*.fits') # find all the fits files in the directory
        
        self.exp_spectra = []
        for path in self.exp_paths:
            wl, fl, ivar = self.read_file(path)
            self.exp_spectra.append(Spectrum(wl, fl, ivar)) # read in the spectra for each coadd

        # perform a coadd by summing the spectra
        self.spectrum = self.exp_spectra[0] # set the spectrum to the first exposure
        if len(self.exp_spectra) > 1: # if there is only one exposure, end
            for i in range(len(self.exp_spectra)-1):
                # now coadd each exposure
                self.spectrum += self.exp_spectra[i+1]

    def read_file(self, path, use_mask=False):
        coadd = fits.open(path)[1]
        cl = coadd.data['wave'] > 0
        wave = coadd.data['wave'][cl]
        fl = coadd.data['flux'][cl]
        ivar = coadd.data['ivar'][cl]

        if use_mask:
            mask = coadd.data['mask'][cl]
            return wave, fl, ivar, mask
        
        return wave, fl, ivar


    def fit_rv(self, model, plot=False):
        rv, e_rv, redchi, param_res = self.spectrum.fit_rv(model)
        if plot:
            self.spectrum.plot(model, param_res.params)
        return rv, e_rv, redchi, param_res
        

class Observation:
    """
    assume that the file structure for n coadds is is something like:
    path/
        obsfile.csv
        target1/
            coadd1
            ...
            coaddn
        target2/
            coadd1
            ...
            coaddn
    """
    def __init__(self, path):
        self.table = Table.read(path+'obsfile.csv') # read the observation file
        self.paths = [path + name for name in self.table['name']] # build paths to each target
        self.targets = [Target(name, path) for name, path in zip(self.table['name'], self.paths)] # create the target object

    def fit_rvs(self, model, save_column = False, verbose=False):
        if verbose:
            print('Measuring RVs')
        rvs = np.array([target.fit_rv(model, verbose) for target in tqdm(self.targets)]) # measure the rvs for each target

        if save_column:
            # if save column, add it to the observation table
            self.table['rv'] = rvs[:,0] 
            self.table['e_rv'] = rvs[:,1]
            self.table['chisqr'] = rvs[:,3]
        return rvs

    def write(self, outfile, **kwargs):
        self.table.write(outfile, **kwargs)

if __name__ == '__main__':
    # read in the arguments from the cli
    parser = argparse.ArgumentParser(prog='spectromancer.py',
                    description='spectromancer.py is a wizard for managing MAGE spectrum observations of white dwarfs')
    parser.add_argument("path", help='path to the observation folder')
    parser.add_argument('--measure-rvs', default=False, action='store_true', help='append RVs to the observation table?')
    parser.add_argument('-o', '--outfile', nargs='?', default=None, help='file to which to write the observation table')
    parser.add_argument('-v', '--verbose', action='store_true', help='display plots?')
    args = parser.parse_args()

    model = corv.models.make_balmer_model(nvoigt=2)
    observation = Observation(args.path)

    if args.measure_rvs:
        observation.fit_rvs(model, verbose=args.verbose)

    if args.outfile is not None:
        observation.write(args.outfile)
    