"""Module for auto-processing Varian NMR fid files for the IDREAM project.
"""

# General Python imports.
import glob
import os
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import multiprocessing as mp
import re
from tqdm import tqdm


# NMR imports.
import nmrglue as ng


def trim_file(path):
    return os.sep.join(os.path.normpath(path).split(os.sep)[6:])


def mp_proc_pipe(in_file, output_dir):
    return process_pipe_file(in_file, output_dir)


def pool_nmr_proc(file_list, output_dir, processes=mp.cpu_count() - 1):
    pool = mp.Pool(processes=processes)
    results = [pool.apply_async(mp_proc_pipe, args=(v, output_dir)) for v in file_list]
    results = [p.get() for p in results]
    return results


def batch_nmrpipe_write(file_list, output_dir):
    pipe_list = list()
    for x in tqdm(file_list):
        pipe_list.append(write_varian_as_pipe(x, output_dir))
    return pipe_list


def read_varian_as_nmrpipe(fid_file):
    """Reads a Varian .fid file and converts it into an NMR pipe file.
    """
    dic, data = ng.varian.read(fid_file)

    udic = ng.varian.guess_udic(dic, data)

    C = ng.convert.converter()

    C.from_varian(dic, data, udic)

    dic, data = C.to_pipe()

    return dic, data


def write_varian_as_pipe(fid_file, output_folder):
    """Reads a Varian .fid file and writes it as an NMR pipe file in the
    output_folder directory.
    """
    # Get the basename of the fid_file.
    # base_name = os.path.basename(fid_file)
    base_name = os.sep.join(os.path.normpath(fid_file).split(os.sep)[5:])

    dic, data = ng.varian.read(fid_file)
    udic = ng.varian.guess_udic(dic, data)
    convert = ng.convert.converter()
    convert.from_varian(dic, data, udic)
    output_path = os.path.join(output_folder, f"pipe-{base_name}")
    ng.pipe.write(output_path, *convert.to_pipe(), overwrite=True)
    return output_path


def process_pipe_file(fid_file, output_folder):
    """"""
    # read in the file
    dic, data = ng.pipe.read(fid_file)

    # Get the base file name without an extension.
    base = os.sep.join(os.path.normpath(fid_file).split(os.sep)[5:])
    base = os.path.splitext(base)[0]

    # Apodization.
    # dic, data = ng.pipe_proc.sp(dic, data)

    # Zero fill.
    dic, data = ng.pipe_proc.zf(dic, data, auto=True)

    # Complex Fourier transform.
    dic, data = ng.pipe_proc.ft(dic, data, auto=True)

    # Phase shift.
    # dic, data = ng.pipe_proc.ps(dic, data, p0=-17.7, p1=-36.0)
    data = ng.process.proc_autophase.autops(data, 'acme')

    # Delete imaginaries.
    dic, data = ng.pipe_proc.di(dic, data)

    # write out processed data
    output_path = os.path.join(output_folder, f"processed-{base}.ft")
    ng.pipe.write(output_path, dic, data, overwrite=True)
    return output_path


def process_fid(dic, data):
    """Processes an NMR pipe file.

    This runs a fourier transform, autophases the datra, then removes
    the imaginary portion of the spectra.
    """
    uc = ng.pipe.make_uc(dic, data)

    data = ng.process.proc_base.fft(data)
    data = ng.process.proc_autophase.autops(data, 'acme')


    dic, data = ng.pipe_proc.di(dic, data)

    return dic, data, uc


def find_nmr_peaks(data, intensity=500000):
    """
    """
    peaks = ng.analysis.peakpick.pick(data, intensity)

    return peaks


def integrate_peak_by_linewidth(peaks, dic, data, uc):
    """Integrates a spectra based on the linewidth of its peaks.
    """
    integrals = dict()
    # TODO: Create a unit conversion object.
    for peak in peaks:

        peak_center = peak[0]
        peak_id = peak[1]
        peak_width = peak[2]

        int_start = (peak_center - peak_width) / 500
        int_end = (peak_center + peak_width) / 500

        integrals[peak_id] = ng.analysis.integration.integrate(
            data, uc, [int_start, int_end])

    return integrals


def main():
    """Main function.
    """
    pass


if __name__ == "__main__":
    main()
