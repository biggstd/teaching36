# Bulk Processing of NMR Data

## Questions

1. What data is not included in these .fid files? That is, what other metadata is needed?


## Overview of Processing


1. Read in the varian .fid files.
2. Read / guess the metadata based on the .fid files.
3. Convert to an nmr_pipe file.
4. Run an NMR fourier transform.
5. Autophase the data.
6. Remove the imaginary points.
7. Pick peak by intensities.
   a. Peak Center.
   b. Peak ID.
   c. Peak width / line width. How is this calculated?
8. Integrate peaks by using the line width.


