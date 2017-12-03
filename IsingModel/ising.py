"""
==========================
Ising Class Implementation
==========================

By: Tyler Biggs
"""

# Import basic packages.
import numpy as np
import matplotlib
import matplotlib.animation as manimation
from tqdm import tqdm


class Ising:
    """
    This class models a 2-Dimensional Ising array.
    """

    def __init__(self, array_size, initial_state='random'):
        """
        Initialization of the Ising array.
        """

        # Create class instance variables.
        self.array_size = array_size
        self.ising_array = self.generate_initial_array(initial_state)

    def generate_initial_array(self, initial_state):
        """
        Generates an array. Optional methods of initial array
        generation can be added under this function.
        """

        # Check if the initial state is random.
        if initial_state == 'random':
            new_ising_array = np.random.randint(
                0, 2, (self.array_size, self.array_size))
            # Set all instances of 0 to -1.
            new_ising_array[new_ising_array == 0] = -1

        # Otherwise return an array of all ones.
        else:
            new_ising_array = np.ones((self.array_size, self.array_size))

        return new_ising_array

    @staticmethod
    def calc_Boltz_Dist_Prob(energy, temperature, k=1.0):
        """
        Returns the probability of a given energy. This function
        defaults to a unit-less output.

        Implemented as a static method so that this function can
        be called without running the initialization above.
        """
        return np.exp(-energy / (k * temperature))

    def calc_spin_flip_energy(self, i, j):
        """
        Calculate the energy, using a simplified Hamiltonian, of
        the current class instance array if the spin at index i, j
        is flipped.
        """
        new_energy = -1 * self.ising_array[i, j] * (
            self.ising_array[i + 1, j] +
            self.ising_array[i - 1, j] +
            self.ising_array[i, j - 1] +
            self.ising_array[i, j + 1]
        )

        return new_energy

    def apply_periodic_boundry(self, n):
        """
        Apply periodic boundary conditions to a set of input
        index values.
        """
        # If the given index is on the right-most edge, return 0,
        # which is the beginning array index.
        if n + 1 >= self.array_size:
            return 0

        # If the given index is on the left-most edge, return the
        # size - 1, which is the right most edge index of the array.
        if n < 0:
            return self.array_size - 1

        # Otherwise return the given index unchanged.
        else:
            return n

    def simulate_epoch(self, temp):
        """
        Simulate a single epoch.
        """

        # Generate two random indices.
        r1, r2 = np.random.randint(0, self.array_size, 2)

        # Apply a given boundary condition.
        r1 = self.apply_periodic_boundry(r1)
        r2 = self.apply_periodic_boundry(r2)

        # Calculate the energy of the system if the spin at
        # r1, r2 is flipped.
        spin_flip_energy = self.calc_spin_flip_energy(r1, r2)

        # If the flip results in a lower energy, simply flip the spin.
        if spin_flip_energy < 0:
            self.ising_array[r1, r2] *= -1

        # Otherwise, generate a new random value between 0 and 1.
        # If the calculated probability is greater than the random
        # number, flip the spin at (r1, r2).
        elif (self.calc_Boltz_Dist_Prob(spin_flip_energy, temp) >
              np.random.rand()):
            self.ising_array[r1, r2] *= -1

    def simulate(self, epochs, temperature, epoch_size=500, video=True):
        """
        Run a simulation for a given number of epochs and temperature.
        """

        FFMpegWriter = manimation.writers['ffmpeg']
        writer = FFMpegWriter(fps=15)

        fig = matplotlib.pyplot.figure()

        with writer.saving(fig, "ising.mp4", 220):

            # Iterate through the range of epochs.
            for epoch in tqdm(range(epochs)):

                # Run an epoch.
                self.simulate_epoch(temperature)

                # If this is an interval of epoch_size, record a frame.
                if epoch % epoch_size == 0:

                    if video:
                        matplotlib.pyplot.axis('off')
                        frame = matplotlib.pyplot.imshow(
                            self.ising_array,
                            cmap='copper',
                            interpolation='nearest')
                        writer.grab_frame()
                        frame.remove()

        matplotlib.pyplot.close('all')
