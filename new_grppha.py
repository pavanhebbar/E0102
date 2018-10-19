#!/Users/pavanrh/anaconda3/bin/python
"""Python program for optimal binning of spectra."""

import os
import argparse
import numpy as np
from astropy.io import fits


def make_parser():
    """Make parser."""

    desc = """Prog. to group pha file for Xspec over a energy range.

              User must specify input and ouput file.
              User can specifiy minimum binning based on spectral resolution
              and minimum number of counts per bin."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-i', '--input',
                        help='input pi file',
                        dest='infile')
    parser.add_argument('-r', '--rmf',
                        help='rmf file',
                        dest='rmf')
    parser.add_argument('-o', '--output',
                        help='output (grouped) pi file',
                        dest='output_pi')
    parser.add_argument('-e', '--energy_bin',
                        help='minimum energy change per bin in eV or fraction [%default]',
                        dest='bin_min_energy',
                        type=float,
                        default=100.0)
    parser.add_argument('-f', '--fractional_energy_binning',
                        help='switch to set fractional energy binning [%default]',
                        dest='fraction_switch',
                        default=False,
                        action='store_true')
    parser.add_argument('-l', '--lower_energy',
                        help='minimum energy of first bin in eV [%default]',
                        dest='lower_energy',
                        type=float,
                        default=300.0)
    parser.add_argument('-u', '--upper_energy',
                        help='maximum energy of last bin in eV [%default]',
                        dest='upper_energy',
                        type=float,
                        default=10000.0)
    parser.add_argument('-c', '--count_bin',
                        help='minimum counts per bin [%default]',
                        dest='bin_min_counts',
                        type=float,
                        default=15)

    return parser


def proc_args(args):
    """Check validity, maniuplate inputs and display errors."""
    if args.infile is None:
        print("Input file is not specified.")
        raise SystemExit

    if not args.infile[-3:] == ".pi":
        print("Not a valid input. Valid input file ends with .pi")
        raise SystemExit

    if not os.path.isfile(args.infile):
        print("Input file doesn't exist")
        raise SystemExit

    if args.rmf is None:
        try:
            args.rmf = fits.open(args.infile)[1].header['RESPFILE']
        except KeyError:
            print("Keyword 'RESPFILE' not found in input file header")
            raise SystemExit

    if not os.path.isfile(args.rmf):
        print("Rmf file doesn't exist")
        raise SystemExit

    if args.output_pi is None:
        args.output_pi = args.infile[:-3] + '_grouped' + '.pi'

    if os.path.isfile(args.output_pi):
        query = input("Output file already exists. Do you want to replace it?")
        if query == 'no':
            raise SystemExit

    return args


def group_pha(indata, ebounds, e_limits, minsize, mincount, frac_bin):
    """Group the pha file.

    Inputs:
    indata: input fits data
    ebounds: EBOUNDS block of the rmf file
    e_limits: [lower energy limit, upper energy limit]
    minsize: minimum bin size
    mincount: minimum counts in each group
    frac_bin: bool variable for fractional energy bin
    """
    grouping = np.zeros(len(indata)) - 1
    quality = np.zeros(len(indata)) + 5
    grouping[0] = 1

    channels = ebounds.data['channel']
    emin = ebounds.data['e_min']
    emax = ebounds.data['e_max']
    energy = (emin*emax)**0.5
    ewidth = emax - emin
    min_ewidth = np.zeros(len(indata)) + minsize

    if frac_bin:
        min_ewidth = min_ewidth*energy

    group_counts = 0
    group_width = 0

    for rownum, row in enumerate(indata):
        pi_channel = row[0]
        pi_counts = row[2]
        match_index = np.where(channels == pi_channel)[0][0]
        if emin[match_index] > e_limits[0] and emax[match_index] < e_limits[1]:
            quality[rownum] = 0
            if group_width == 0:
                grouping[rownum] = 1

            if group_counts > mincount and group_width > min_ewidth[match_index]:
                group_counts = pi_counts
                group_width = ewidth[match_index]
                grouping[rownum] = 1
            else:
                group_counts += pi_counts
                group_width += ewidth[match_index]

    if group_counts < mincount:
        lastgroup = np.where(grouping == 1)[0][-1]
        grouping[lastgroup] = -1

    return quality, grouping


def process_fits(args):
    """Process the parser and fits files and feed into group_pha."""
    infits = fits.open(args.infile)
    indata = infits[1].data
    ebounds = fits.open(args.rmf)['ebounds']
    quality, grouping = group_pha(indata, ebounds, [args.lower_energy,
                                                    args.upper_energy],
                                  args.bin_min_energy, args.bin_min_counts,
                                  args.fraction_switch)
    group_col = fits.Column(name='GROUPING', array=grouping, format='I')
    qual_col = fits.Column(name='QUALITY', array=quality, format='I')
    outdata = fits.BinTableHDU.from_columns(indata.columns + group_col
                                            + qual_col)
    infits[1].data = outdata.data
    infits.writeto(args.output_pi, overwrite=True)


def main():
    """Execute all other functions."""
    parser = make_parser()
    args = parser.parse_args()
    args = proc_args(args)
    print(args)
    process_fits(args)


if __name__ == '__main__':
    main()
