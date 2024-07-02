import argparse

import corv
from .spectromancer import *

# read in the arguments from the cli
parser = argparse.ArgumentParser(prog='spectromancer.py',
                description='spectromancer.py is a wizard for managing MAGE spectrum observations of white dwarfs')
parser.add_argument("path", help='path to the observation folder')
parser.add_argument('--measure-rvs', default=False, action='store_true', help='append RVs to the observation table?')
parser.add_argument('-o', '--outfile', nargs='?', default=None, help='file to which to write the observation table')
parser.add_argument('-v', '--verbose', action='store_true', help='display plots?')
args = parser.parse_args()

model = corv.models.make_warwick_da_model(names=['a','b','g','d'])
observation = Observation(args.path)

if args.measure_rvs:
    observation.fit_rvs(model, save_column=True, verbose=args.verbose)

if args.outfile is not None:
    observation.write(args.outfile)