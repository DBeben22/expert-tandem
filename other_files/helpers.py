import math as math
from xml.etree.ElementTree import fromstring

import numpy as np
import matplotlib.pyplot as plt



"""
code from:
https://github.com/PASTAQ-MS/PASTAQ
"""
# plot mass spectra considering that only non-eros values are included in the spectra
def plot_msSpectra(mz, intensity, norm_mz_diff = 0.0035, diffFactor = 1.3, scanIdx = 0, mzRange = (0, 0)):
    """
    Plotting mass spectra, expecting spectra where 0 intensity values were omitted.

    Args:
        mz: A list of mz values of the mass spectra.
        intensity: A list of intensity (non-zero) values of the mass spectra.
        norm_mz_diff: Difference between two adjacent mz at measurement point corresponding to the original sampling frequency of the mass spectra.
        diffFactor: Tolerance factor allowing to vary sampling frequency of mass spectra. It is typically set to 30% (factor 1.3).

    Returns:
        Figure object as fig.
    """
    spectra = {'mz': mz, 'intensity': intensity}
    newSpectra = {'mz': [], 'intensity': []}
    idxSpectra = 0
    for i in range(1, len(spectra['mz'])):
        diff = spectra['mz'][i ] -spectra['mz'][ i -1]
        # takes the m/z difference between 2 adjacent values
        # and if that difference exceeds the set norm_mz_diff * diffFactor
        if diff > norm_mz_diff *diffFactor:
            if diff < norm_mz_diff *2:
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][ i -1] + diff /2)
                newSpectra['intensity'].insert(idxSpectra, 0)
                idxSpectra += 1
                # print(norm_mz_diff *diffFactor, diff, norm_mz_diff *diffFactor + diff, norm_mz_diff)
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i])
                newSpectra['intensity'].insert(idxSpectra, spectra['intensity'][i])
                idxSpectra += 1
            else:
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][ i -1] + norm_mz_diff)
                newSpectra['intensity'].insert(idxSpectra, 0)
                idxSpectra += 1
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i] - norm_mz_diff)
                newSpectra['intensity'].insert(idxSpectra, 0)
                idxSpectra += 1
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i])
                newSpectra['intensity'].insert(idxSpectra, spectra['intensity'][i])
                idxSpectra += 1
        else:
            newSpectra['mz'].insert(idxSpectra, spectra['mz'][i])
            newSpectra['intensity'].insert(idxSpectra, spectra['intensity'][i])
            if (diff /diffFactor) > norm_mz_diff:
                norm_mz_diff = diff
            idxSpectra += 1


    if mzRange == (0,0):
        mzRange = (min(newSpectra['mz']), max(newSpectra['mz']))

    fig = plt.figure(figsize=(25, 6), facecolor='white')  # Set the figure size and white background
    plt.plot(newSpectra['mz'], newSpectra['intensity'], color = 'red', marker='', linestyle='-')  # Plot mz vs. intensity
    plt.xlabel('m/z')  # Set the x-axis label
    plt.ylabel('Intensity')  # Set the y-axis label
    plt.xlim(mzRange)
    plt.title('Mass Spectrum of scan {}' .format(scanIdx))  # Set the title
    plt.grid(False)  # Show grid
    # return fig