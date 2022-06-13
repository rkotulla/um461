#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits as pyfits
import scipy.ndimage

if __name__ == "__main__":

    catalog_fn = sys.argv[1]
    reference_filename = sys.argv[2]
    column_number = int(sys.argv[3])
    output_filename = sys.argv[4]

    try:
        # read catalog
        catalog = numpy.loadtxt(catalog_fn)
        # extract the data from the right catalog column
        data = catalog[:, column_number-1]
        pixel_index = catalog[:, 0]

    except:
        data = []
        # open the file
        pixel_index = []
        with open(catalog_fn) as file:
            # read all the lines, and ave them as an array of lines
            lines = file.readlines()

            # now go over every single line
            for line in lines:
                # if it's a comment, ignore it
                if (line.startswith("#")):
                    continue

                # split lines into columns, using whitespace as a separator
                items = line.split()

                # convert just the single column item into a number
                value = float(items[column_number-1])
                _pi = int(items[0])

                # and add it to our data array
                data.append(value)
                pixel_index.append(_pi)

        # once that's done for all lines, convert the array into a proper numpy array
        data = numpy.array(data)
        pixel_index = numpy.array(pixel_index)

    print(data.shape)

    # sort data by pixelindex
    pi_sort = numpy.argsort(pixel_index)
    data = data[pi_sort]

    reference = pyfits.open(reference_filename)
    hdr = reference[0].header

    data2d = data.reshape((hdr['NAXIS1'], hdr['NAXIS2']))
    print(data2d.shape)

    # put the newly reconstructed image data back into the original FITS file
    reference[0].data = data2d.T

    print("writing image back to %s" % (output_filename))
    reference.writeto(output_filename, overwrite=True)


