import numpy as np


class gillespie:

    def __init__(self, X_values, c_values, X_change, h_form, max_time, max_count):

        self.X_values = X_values
        self.c_values = c_values
        self.X_change = X_change
        self.h_form = h_form

        # some stopper variables so the simulation does not go on forever.
        self.max_time = max_time
        self.max_count = max_count

        # start at one with the count.
        self.count = 1
        self.time = 0

    # first thing to calculate is av = hv * cv
    # where v is an index
    def av_calc(self):
        return [self.h_form[v](self.X_values) * self.c_values[v - 1] for v in self.h_form]

    # Then calculate a0, which is the sum of the av values
    def a0_calc(self, avals):
        return sum(avals)

    # then tau, which is the probable length of time before any one
    # of the monitored reactions occurs.
    # r1 will be a random number we supply the function with
    # r1 will be between 0 and 1.0
    def tau_calc(self, asum, r1):
        return (1. / asum) * np.log(1. / r1)

    # next is mu, the random value
    # r2 is also between 0 and 1.0
    def mu_calc(self, avals, asum, r2):
        # start an index at 1
        u = 1
        # while r2 * ao is less than the successive sums of av
        # increase the u value by 1
        while sum(avals[:u]) < r2 * asum:
            u = u + 1
        return u

    # now define a simulation function
    def simulate(self, ):

        # an a list for output with the starting time, should be 0
        time_out = [0]
        # an empty list that contains the X values, include the starting here, as it will change after
        # the first pass of the algorithm.
        X_out = [self.X_values]

        # while we haven't hit the heat death of the universe...
        while self.time < self.max_time or self.max_time == None and self.count < self.max_count:

            # calculate av
            av = self.av_calc()
            # print "1av: ", av

            # calculate a0
            a0 = self.a0_calc(av)
            # print "2a0: ", a0

            # generate random numbers
            r1, r2 = np.random.random(), np.random.random()

            # calculate tau
            tau = self.tau_calc(a0, r1)
            # print "3tau", tau

            # calculate mu
            mu = self.mu_calc(av, a0, r2)
            # print "4mu: ", mu

            # call the change function, and set
            self.X_values = self.X_change[mu](self.X_values)

            self.time = self.time + tau
            # print 'time: ', self.time
            time_out.append(self.time)
            X_out.append(self.X_values)

            self.count = self.count + 1

        return time_out, X_out
