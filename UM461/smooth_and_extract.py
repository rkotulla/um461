#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits as pyfits
import scipy.ndimage

if __name__ == "__main__":

    _x = sys.argv[1]
    _y = sys.argv[2]
    # print(_x)
    cutout_x = [int(i) for i in _x.split(":")]
    cutout_y = [int(i) for i in _y.split(":")]
    # print(cutout_x, cutout_y)

    binning = int(sys.argv[3])
    filtering_length = int(sys.argv[4])

    image_filenames = sys.argv[5:-1]

    multiband_catalog = None

    for fn in image_filenames:
        print(fn)

        # open filename
        image_file = pyfits.open(fn)
        # image_file.info()

        # extract image data
        img_data = image_file[0].data.T


        # print how big the image is (in x,y)
        print(img_data.shape)

        # now let's do the median filtering
        print("Filtering data, that'll take a while, be patient, young padawan...")
        filtered_image = scipy.ndimage.median_filter(
            input=img_data,
            size=filtering_length,
            mode='nearest'
        )
        print(filtered_image.shape)

        print("saving files")

        dump_file = pyfits.PrimaryHDU(data=filtered_image.T, header=image_file[0].header)
        dump_file.writeto(fn[:-5]+"_filtered.fits", overwrite=True)

        highpass_data = img_data - filtered_image
        dump2_file = pyfits.PrimaryHDU(data=highpass_data.T, header=image_file[0].header)
        dump2_file.writeto(fn[:-5]+"_highpass.fits", overwrite=True)

        # now we have the media-fitlered data, continue to bin it up
        # but first, create cutout

        x = cutout_x
        y = cutout_y
        cutout_image = filtered_image[x[0]-1:x[1]-1, y[0]-1:y[1]-1]
        print(cutout_image.shape)

        print("binning")
        binned_size_x = (x[1] - x[0]) // binning
        binned_size_y = (y[1] - y[0]) // binning

        rb1 = numpy.reshape(cutout_image, (binned_size_x, binning, binned_size_y, binning))
        rb2 = numpy.nansum(rb1, axis=-1)
        rebinned = numpy.nansum(rb2, axis=1)
        print(rebinned.shape)

        binned_file = pyfits.PrimaryHDU(data=rebinned.T, header=image_file[0].header)
        binned_hdr = binned_file.header
        # fix the WCS to account for the cutout
        binned_hdr['CRPIX1'] -= x[0]
        binned_hdr['CRPIX2'] -= y[0]
        # and now account for the binning
        binned_hdr['CRPIX1'] /= binning
        binned_hdr['CRPIX2'] /= binning
        binned_hdr['CD1_1'] *= binning
        binned_hdr['CD1_2'] *= binning
        binned_hdr['CD2_1'] *= binning
        binned_hdr['CD2_2'] *= binning
        binned_file.writeto(fn[:-5]+"_binned.fits", overwrite=True)

        hdr = image_file[0].header
        mag_zeropoint = -2.5*numpy.log10(hdr['PHOTFLAM']) - 5*numpy.log10(hdr['PHOTPLAM']) - 2.408
        print("Using zeropoints of %f for %s" % (mag_zeropoint, fn))
        mag_image = -2.5 * numpy.log10(rebinned) + mag_zeropoint

        phot_error = numpy.sqrt( rebinned * hdr['EXPTIME'] )
        mag_error = phot_error / (rebinned * hdr['EXPTIME'])

        mag_and_errors_combined = numpy.array( [
            mag_image.ravel(), mag_error.ravel()
        ]).T
        print(mag_and_errors_combined.shape)

        ids_for_gazelle = numpy.arange(mag_and_errors_combined.shape[0]).reshape((-1,1))

        if (multiband_catalog is None):
            # on the first go-around, just use the IDs
            multiband_catalog = ids_for_gazelle

        # and add all the photometry and errors
        multiband_catalog = numpy.append(multiband_catalog, mag_and_errors_combined, axis=1)

        print(multiband_catalog.shape)
        print("done, yay (well, with this file at least)")
        #break

    print("all files completely done, writing output catalog")
    output_catalog = sys.argv[-1]
    numpy.savetxt(output_catalog, multiband_catalog)


