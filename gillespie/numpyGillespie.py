#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd


class Gillespie:

    def __init__(self, species, rates, species_changes, permutations,
                 max_sim_rxn):
        """
        Assigns the input values to appropriate numpy arrays.

        :param species:
            A list of the species to be tracked. Should be integers.

        :param rates:
            A list of the rates to be simulated.

        :param species_changes:
            Changes to the species, indexed the same as the reaction
            rate list.

        :param permutations:
            A list of functions that calculate the permutations
            available for each of the reactions tracked. Should be
            indexed the same as `species_change` and `rates`.

        :param max_sim_time:
            The maxiumum simulation time the program will be allowed
            to run. This time depends on the scale chosen in the rates.

        :param max_sim_rxn:
            The maximum number of reactions that the program will be
            allowed to simulate.

        """

        # Set the maximum number of allowed simulations.
        self.max_sim_rxn = max_sim_rxn

        # Define numpy arrays of the species and rates inputs.
        self.species = np.array(species)
        self.rates = np.array(rates)

        # Create an array for the time output, the size is determined
        # by the maximum number of simulated events, as set by
        # `self.max_sim_rxn`.
        self.time_out = np.empty(shape=self.max_sim_rxn, dtype=np.float)

        # Create an array for the species count output. The same length
        # dimension is used as above, but species to track as well.
        self.species_out = \
            np.empty(shape=(self.max_sim_rxn, len(self.species)))

        # Append / set the initialization values to the arrays created.
        self.time_out[0] = 0
        self.species_out[0] = self.species

        # Set the reaction count index to one (to allow for the starting)
        # values to be recoreded.
        self.rxn_count = 1

        # Assign the inputs to the class instance.
        self.species_changes = species_changes
        self.permutations = permutations

    def calc_av(self, curr_species):
        """
        Calculate the Av value.

        The current permutations multiplied be the rates.

        Will use `np.ma.fromiter(iterable, dtype, count=-1)`. Which
        creates a new 1-dimensional array from an iterable object.
        """

        # Create an generator that is iterable from the lambda
        # functions found within `self.permutations`.
        iterator = [fn(curr_species) * self.rates[i]
                    for i, fn in enumerate(self.permutations)]

        # Return a numpy array of floats from the iterator above.`
        return np.fromiter(iterator, np.float)

    def calc_tau(self, Av_sum, random_value):
        """
        Calculate the Tau value, which is the probable length of time
        before any given simulated reaction occurs.

        See the Gillespie paper for a discussion of this.
        """
        return (1. / Av_sum) * np.log(1.0 / random_value)

    def calc_mu(self, Av_vals, Av_sum, random_value):
        """
        Calculate the mu value.

        :param Av_vals:
            The possible reaction permutations * their rates. Given
            as a numpy array.

        :param Av_sum:
            The sum of the `Av_vals` array. This sum is used in
            multiple places, so it is not calculated within this
            function.

        :param random_value:
            A randomly generated value between zero and one.

        :returns mu:
            The index corresponding to the reaction that should be
            simulated.

        Essentially we generate blocks of value on a number line
        from zero to one. A random number determines where on this
        line a reaction "occurs".

            [================================================]
            0.0                  *                          1.0
                                  A random point.

        We cast the possible reactions to this scale by multiplying
        the random value by the sum of Av values. Such casting is
        done by chunks.

            [================================================]
            [=Chunk 1=][======Chunk 2======][=====Chunk 3====]

        The sums of these chucnks are examined iteratively, and
        when the sum is found to be greated than the randomly cast
        point defined above, the corresponding reaction is simulated.

        ..warning::
            A different random value should be used for `calc_mu()` and
            `calc_tau()`.
        """

        # Initialize the integer mu as 0.
        mu = 0

        # Cast the Av_sum to the random_value scale.
        cast = Av_sum * random_value

        # While the current sum is below the randomly cast value,
        # increase the Av_vals sum range by one and check again.
        # The + 1 below to the mu index moves the sum to the
        # correct position, and maintains the integrity of the mu
        # index.
        while sum(Av_vals[:mu + 1]) < cast:
            mu += 1

        # When the sum is exceeded, return the mu index.
        return mu

    def simulate(self):
        """
        Runs a stochastic simulation based on the input provided.

        See the Gillespie paper for a discussion of this.

        :returns:
            A tuple of `(time, species)`.
        """

        # Begin the simulation.
        # Check if the maximum number of allowed reaction has been reached.
        while self.rxn_count < self.max_sim_rxn:

            # Add the current species counts to the output list.
            self.species_out[self.rxn_count] = self.species
            # print("current_species", self.species)

            # Add the current time to the output list.
            self.time_out[self.rxn_count] = self.time

            # Calculate the Av values.
            Av = self.calc_av(self.species)
            # print("Av", Av)

            # Sum the Av values.
            Av_0 = sum(Av)
            # print("Av_0", Av_0)

            # Generate two random numbers.
            r_one, r_two = np.random.random(), np.random.random()

            # Calculate tau, the time until a reaction occurs. Then
            # increment the current time by that value.
            tau = self.calc_tau(Av_0, r_one)
            self.time += tau

            # Determine mu, the index of which reaction occurs.
            mu = self.calc_mu(Av, Av_0, r_two)

            # Change the species by the function indexed by mu. The current
            # index is given by the current reaction count.
            self.species_out[self.rxn_count] = \
                self.species_changes[mu](current_species)

            # Increment the reaction counter by one.
            self.rxn_count += 1

            # Calculate the new population of species based on the mu index.
            self.species = self.species_changes[mu](self.species)

        # When the loop is over (the maximum number of reactions to be
        # simulated has been reached) return a tuple of the time and
        # species values.
        return (self.time_out, self.species_out)


class CompleteGillespie(Gillespie):
    """
    Assigns the input values to appropriate numpy arrays. This
    version of the Gillespie algorithm implementation returns
    much more information.

    The output from CompleteGillespie.simulate() is a dictionary
    of arrays.
    """

    def __init__(self, *args, **kwargs):

        """
        This `CompleteGillespie` initialization creates new output
        arrays to return more data concerning the reactions simulated.
        
        :param species:
            A list of the species to be tracked. Should be integers.

        :param rates:
            A list of the rates to be simulated.

        :param species_changes:
            Changes to the species, indexed the same as the reaction
            rate list.

        :param permutations:
            A list of functions that calculate the permutations
            available for each of the reactions tracked. Should be
            indexed the same as `species_change` and `rates`.

        :param max_sim_time:
            The maxiumum simulation time the program will be allowed
            to run. This time depends on the scale chosen in the rates.

        :param max_sim_rxn:
            The maximum number of reactions that the program will be
            allowed to simulate.
        """
        
        # Call the parent classes __init__() function, pass any
        # appropriate *args or **kwargs.
        super().__init__(*args, **kwargs)

        # Create the advanced output arays.
        self.av_out = np.empty(shape=(self.max_sim_rxn, len(self.rates)))
        self.mu_out = np.empty(self.max_sim_rxn)

    def simulate(self):
        """
        Runs a stochastic simulation based on the input provided.

        See the Gillespie paper for a discussion of this.

        :returns:
            A tuple of `(time, species)`.
        """

        # Begin the simulation.
        # Check if the maximum number of allowed reaction has been reached.
        while self.rxn_count < self.max_sim_rxn:

            # Add the current species counts to the output list.
            self.species_out[self.rxn_count, :] = self.species
            # print("current_species", self.species)

            # Add the current time to the output list.
            self.time_out[self.rxn_count] = self.time

            # Calculate the Av values.
            Av = np.array(self.calc_av(self.species))
            self.av_out[self.rxn_count, :] = Av
            # print("Av", Av)

            # Sum the Av values.
            Av_0 = sum(Av)
            # print("Av_0", Av_0)

            # Generate two random numbers.
            r_one, r_two = np.random.random(), np.random.random()

            # Calculate tau, the time until a reaction occurs. Then
            # increment the current time by that value.
            tau = self.calc_tau(Av_0, r_one)
            self.time += tau

            # Determine mu, the index of which reaction occurs.
            mu = self.calc_mu(Av, Av_0, r_two)
            self.mu_out[self.rxn_count] = mu

            # Increment the reaction counter by one.
            self.rxn_count += 1

            # Calculate the new population of species based on the mu index.
            self.species = self.species_changes[mu](self.species)

        # When the loop is over (the maximum number of reactions to be
        # simulated has been reached) return a tuple of the time and
        # species values.
        output_dict = {
            "time": self.time_out,
            "species": self.species_out,
            "av": self.av_out,
            "mu": self.mu_out
        }
        return output_dict


def pandas_output(out_dict):
    """
    Creates a pandas dataframe by iterating over selected dictionary
    entries that have more than one dimensions in their output arrays.

    :param out_dict:
        A dictionary provided by `CompleteGillespie.simulate()`

    :returns:
        A pandas dataframe, with enumerated columns generated from the
        multidimensional arrays `species` and `av`.
    """

    # Create the output data frame, and append the 1-dimensional
    # arrays and associated keys.
    df = pd.DataFrame()
    df['time'] = out_dict['time']
    df['mu'] = out_dict['mu']

    # Get the length of the species sub-lists.
    numb_species = len(out_dict['species'][0])

    # Iterate over this length and pull out individual species.
    for i in range(numb_species):

        # Construct the key for the output dictionary.
        new_species_key = "species_{}".format(i)

        # Use the key to assign the data to the output dictionary.
        df[new_species_key] = [s[i] for s in out_dict['species']]

    # Get the length of the av sub-lists.
    numb_av = len(out_dict['av'][0])

    # Iterate over this length and pull out individual species.
    for i in range(numb_av):

        # Construct the key for the output dictionary.
        new_av_key = "av_{}".format(i)

        # Use the key to assign the data to the output dictionary.
        df[new_av_key] = [s[i] for s in out_dict['av']]

    # Return the pandas dataframe constructed above.
    return df
