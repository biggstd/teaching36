"""Module for auto-processing Varian NMR fid files for the IDREAM project.
"""

# General Python imports.
import glob

# NMR imports.
import nmrglue as ng


def read_varian_as_nmrpipe(fid_file):
    """Reads a Varian .fid file and converts it into an NMR pipe file.
    """
    dic, data = ng.varian.read(fid_file)

    udic = ng.varian.guess_udic(dic, data)

    C = ng.convert.converter()

    C.from_varian(dic, data, udic)

    dic, data = C.to_pipe()

    return dic, data


def process_fid(dic, data):
    """Processes an NMR pipe file.

    This runs a fourier transform, autophases the datra, then removes
    the imaginary portion of the spectra.
    """
    data = ng.process.proc_base.fft(data)
    data = ng.process.proc_autophase.autops(data, 'peak_minima')
    dic, data = ng.pipe_proc.di(dic, data)

    return dic, data


def find_nmr_peaks(data, intensity=500000):
    """
    """
    peaks = ng.analysis.peakpick.pick(data, intensity)

    return peaks


def integrate_peak_by_linewidth(peaks, dic, data):
    """Integrates a spectra based on the linewidth of its peaks.
    """
    integrals = dict()

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
