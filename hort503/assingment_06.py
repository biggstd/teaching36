import numpy as np

a = np.arange(15).reshape(3, 5)
a

a.shape

a.ndim

a.dtype.name

a.itemsize

a.size

type(a)

b = np.array([6, 7, 8])

b

type(b)

# Important to remember how arguments are passed in python!

# a = np.array(1, 2, 3, 4) # WRONG!!!

b = np.array([1, 2, 3, 4]) # CORRECT!!

# %%

# Import pyplot plt for convenience.
import matplotlib.pyplot as plt
# Import matplotlib into the notebooks namespace.
import matplotlib
# Use jupyter notebooks inline function to simplify plot display.
%matplotlib inline

# %%

plt.plot([1,2,3,4])
plt.ylabel('Demo Plot.')
plt.show()

# %%
