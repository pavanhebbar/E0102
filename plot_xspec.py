"""Program to plot figures in python from saved models of Xspec.

Creates a temperory tcl file adding few more lines into the saved .xcm file
Then writes the data and the model seperately into different files
Uses python to plot from these figures.
"""

import os
import copy
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator, FixedFormatter, FormatStrFormatter


def makeparser():
    """Construct parser."""
    parser = argparse.ArgumentParser(
        description="""Program to plot saved xspec models in python.

        Creates a temperory tcl file adding few more lines into the saved .xcm file
        Then writes the data and the model seperately into different files
        Uses python to plot from these figures.
        """,
        epilog="""Version - v1.0.

        Last Update - 7 May 2019
        """)
    parser.add_argument('model_in', help='Saved model to be plot. If the input'
                        'is a foldername, all the .xcm files shall pe plotted',
                        type=str)
    parser.add_argument('outroot', nargs='?', default=None, help='Rootname for the output'
                        'files. If no argument is given the he basename of the'
                        'input model will be used.', type=str)
    parser.add_argument('--models', nargs='?', default=None, type=str,
                        help='Names of the additive models used')
    parser.add_argument('--data_style', '--dp', default='k.', type=str,
                        help='Default point style for plotting data points')
    parser.add_argument('--model_style', '--ms', default='-', type=str,
                        help='Default line style for final model')
    parser.add_argument('--comp_style', '--cs', default='--', type=str,
                        help='Default line style for component models')
    parser.add_argument('--legend', default=False, type=bool,
                        help='Whether or not legend must be shown')
    parser.add_argument('--label', default='automatic', type=str,
                        help='Whether to the label names must be "automatic"'
                        'or "ask" for labels')
    return parser


def process_input(args):
    """Manipulate the input arguments and set default values.

    Refer the help statements in the argparser to know what happens to each
    input.

    Inputs:
    args - Arguments passed by the user after being parsed

    Outputs:
    args - Modified argepassed with default values
    """
    # Check if input model exists
    if not os.path.exists(args.model_in):
        print("Input model file doesnot exist")
        raise SystemError

    # Default output name if nothing is specified
    if args.outroot is None:
        args.outroot = os.path.basename(args.model_in)
        if not args.outroot.find('.') == -1:
            args.outroot = args.outroot[:args.outroot.index('.')]

    # Find the names of the additive models used
    if args.models is None:
        infile = open(args.model_in, 'r')
        for line in infile:
            if not line.find('model') == -1:
                models = line.split()
                models.remove('model')
                if not '+' in models:
                    models[0] = models[0][models[0].index('*')+1:]
                    args.fullmo = models[0]
                else:
                    args.fullmo = ''.join(models)
                    args.fullmo = args.fullmo[args.fullmo.index('('):
                                              args.fullmo.index(')')+1]
                    for i, model in enumerate(models):
                        if model == '+':
                            models.remove('+')
                        if '(' in models[i]:
                            models[i] = models[i][models[i].index('(')+1:]
                        if ')' in models[i]:
                            models[i] = models[i][:models[i].index(')')]
        args.models = copy.copy(models)
    args.nmodels = len(args.models)
    if args.nmodels > 3:
        args.fullmo = "Complete model"
    return args


def create_tcl(args):
    """Creates output tcl file to write the data from saved model in infile.

    Output tcl file is deleted in the end
    Data and model files are written sepertely to take more points into
    consideration for model
    Data inculde instrumental effects i.e folded model

    Inputs:
    args - Modified parsed arguments
    """
    in_file = open(args.model_in, 'r')
    outfile = open('temp_'+ args.model_in, 'w+')
    for line in in_file:
        outfile.write(line)
    # outfile.write('ignore 0.0-0.3, 10.0-**\n')
    # outfile.write('notice 0.3-10.0 \n')
    # outfile.write('query yes \n')
    # outfile.write('statistic test cvm \n')
    # outfile.write('abund wilm \n')
    outfile.write('fit \n')
    outfile.write('cpd /xw \n')
    outfile.write('setp en \n')
    outfile.write('setp add \n')
    # outfile.write('setp rebin 2 10 1 \n')
    outfile.write('pl ldata \n')
    os.system('rm -rf ' + args.outroot + '_data.*')
    os.system('rm -rf ' + args.outroot + '_mo.*')
    os.system('rm -rf ' + args.outroot + '_resid.*')
    outfile.write('setp comm we ' + args.outroot + '_data \n')
    outfile.write('plot \n')
    outfile.write('setp delete 1\n')
    outfile.write('plot model \n')
    outfile.write('setp comm we ' + args.outroot + '_mo \n')
    outfile.write('plot \n')
    outfile.write('setp delete 1\n')
    outfile.write('plot delchi \n')
    outfile.write('setp comm we ' + args.outroot + '_resid \n')
    outfile.write('plot \n')
    outfile.write('quit \n')
    outfile.write('y \n')
    outfile.close()


def plot_figure(args):
    """Plot the .qdp files generated previously."""
    data = np.loadtxt(args.outroot + '_data.qdp', skiprows=3)
    #model = np.loadtxt(args.outroot + '_mo.qdp', skiprows=3)
    residue = np.loadtxt(args.outroot + '_resid.qdp', skiprows=3)
    # f, (ax1, ax2) = plt.subplots(2, sharex=True)
    minorloc = FixedLocator([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                             2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    minorform = FixedFormatter(['', '', '', '0.5', '', '', '', '',
                                '2.0', '', '', '5.0', '', '', '', ''])
    majorloc = FixedLocator([0.1, 1.0, 10.0])
    majorform = FormatStrFormatter('%1.1f')
    fig = plt.figure(figsize=(5, 4))
    grids = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax1 = fig.add_subplot(grids[0])
    ax2 = fig.add_subplot(grids[1], sharex=ax1)
    plt.xlabel('Energy (keV)')
    ax1.set_ylabel(r'Normalized counts s$^{-1}$ keV$^{-1}$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    plt.xlim(0.5, 10.0)
    ax1.set_ylim(np.min(data[:, 2])*0.5, np.max(data[:, 2])*2.0)
    ax1.errorbar(data[:, 0], data[:, 2], xerr=data[:, 1], yerr=data[:, 3],
                 fmt=args.data_style)
    xdata = np.append(data[:, 0] - data[:, 1], data[-1, 0] + data[-1, 1])
    ydata = np.append(data[:, 4], data[-1, 4])
    ax1.step(xdata, ydata, where='post', label='Model')
    ax1.tick_params(axis='both', which='both', bottom=True, top=True,
                    left=True, right=True,
                    direction='in', labelbottom=False)
    ax1.tick_params(axis='both', which='major', length=8)
    ax1.tick_params(axis='both', which='minor', length=5)
    ax2.xaxis.set_major_locator(majorloc)
    ax2.xaxis.set_major_formatter(majorform)
    ax2.xaxis.set_minor_locator(minorloc)
    ax2.xaxis.set_minor_formatter(minorform)

    #ax1.plot(model[:, 0], model[:, 2], args.model_style, label=args.fullmo)
    comp_styles = ['g--', 'r-.', 'm:', 'm:', 'm:']
    # if len(data[0, :]) > 5:
    #     for i in range(args.nmodels):
    #         if args.label == 'ask':
    #             args.models[i] = input("Input label for model" +
    #                                    args.models[i])
    #         ax1.plot(data[:, 0], data[:, 5+i], comp_styles[i],
    #                  label=args.models[i])
    # ax1.plot(data[:, 0], data[:, -1], comp_styles[0],
    #           label=r'PL, $\Gamma = 2.0$')
    ax1.plot(data[:, 0], data[:, -1], comp_styles[1],
             label=r'NS, $B = 10^{12}$ G')
    ax1.plot(data[:, 0], data[:, 4] - data[:, -1], #  - data[:, -2],
             comp_styles[2], label='Background')
    ax1.legend()
    ax2.tick_params(axis='both', which='both', bottom=True, top=True, left=True,
                    right=True, direction='in', length=8)
    ax2.tick_params(axis='both', which='major', length=8)
    ax2.tick_params(axis='both', which='minor', length=5)
    ax2.set_ylabel(r'$\Delta \chi$')
    ax2.set_yscale('linear')
    ax2.set_ylim(np.min(residue[:, 2])*1.1, np.max(residue[:, 2])*1.1)
    ax2.errorbar(residue[:, 0], residue[:, 2], xerr=residue[:, 1],
                 yerr=residue[:, 3],
                 fmt=args.data_style, label='data')
    ax2.plot(np.linspace(0.3, 10.0, 10), np.zeros(10))
    ax2.set_yticks(np.arange(int(min(residue[:, 2])), int(max(residue[:, 2]))
                             + 0.5, 1.0, dtype=int))
    if args.legend:
        ax1.legend()
    fig.subplots_adjust(hspace=0)
    plt.setp(ax1.get_xticklabels(), visible=False)
    # plt.tight_layout()
    plt.savefig(args.outroot+'.png')
    plt.show()


def main():
    """Call all functions in order."""
    parser = makeparser()
    args = parser.parse_args()
    args = process_input(args)
    create_tcl(args)
    os.system('heasoft')
    os.system('xspec - ' + 'temp_' + args.model_in)
    plot_figure(args)


if __name__ == '__main__':
    main()
