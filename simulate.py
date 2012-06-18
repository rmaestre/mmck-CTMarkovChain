#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
    This class provides several methods to simulate a m/m/c/k tail

    @author: Roberto Maestre {rmaestre@paradigmatecnologico.com} Paradigmalabs 2012
"""
from __future__ import division
from math import pow
from math import factorial as fac
import nose
import time

class mmck:
    
    def __init__(self, lambda_v, mu_v, c, k):
        """
        Set initial parameters
        """
        # Init parameters to start simulation
        self.l = lambda_v
        self.m = mu_v
        self.c = c
        self.k = k
        self.p0_optimizer = False
        self.p0_aux = 0
        
        # Calculate Rho parameter
        self.rho = self.l / (self.c * self.m)
        
        # Calculate r parameter
        self.r = self.l / self.m
    
    def stability_condition(self):
        """
        Checking system stability condition
        """
        return self.rho < 1
        
        
    def mu(self, mu, n):
        """
        Mu parameter is dinamic
        """
        if 0 < n and n < self.c:
            return n * self.mu
        elif self.c <= n and n <= self.k:
            return c * self.mu
        
    def p0(self):
        """
        Probability in "0" state
        """
        if not self.p0_optimizer:
            self.p0_optimizer = True
            acum = 0
            for n in range(0, self.c):
                acum += pow(self.r, n) / fac(n)
            aux = ((pow(self.r, self.c))/(fac(self.c))) * (((1 - pow(self.rho, self.k - self.c + 1)) / (1 - self.rho)))
            self.p0_aux = 1 / (acum + aux)
            return self.p0_aux
        else:
            return  self.p0_aux
        
        
    def pn(self,n):
        """
        Probability in "n" state, holding n > 0 and n <= k
        """
        assert(n > 0 and n <= self.k)
        
        if n < self.c:
            return (pow(self.l, n) / (fac(n) * pow(self.m, n))) * self.p0()
        elif self.c <= n:
            return pow(self.l, n) / (pow(self.c, (n - self.c)) * fac(self.c) * pow(self.m, n)) * self.p0()
        
        
    def lq(self):
        """
        Mean length og the tail
        """
        a = (self.p0() * pow(self.r, self.c) * self.rho) / (fac(self.c) * pow(1-self.rho, 2))
        b = 1 - pow(self.rho, self.k - self.c + 1) - ((1-self.rho) * (self.k - self.c + 1) * pow(self.rho, self.k - self.c))
        return a * b
        
    def lm(self):
        """
        Mean number of clients into the tail
        """
        return self.lq() + (self.r * 1 - self.pn(self.k))
        
        
    def w(self):
        """
        Mean time of a client waiting into the tail
        """
        return self.lm() / (self.l * (1 - self.pn(self.k)))
        
        
    def wq(self):
        """
        Mean time of a client waiting into the tail
        """
        return (self.lm() / (self.l * (1 - self.pn(self.k)))) - (1 / self.m)


t_start = time.time()
# Create mmck simulator object 
# mmck(l,m,servers,tail limit)
simulator = mmck(2.7,1.5,4,12)

# Some debug options
print "\nM/M/c/K model simulation"
print "------------------------"

print "\n+ MODEL PARAMETERS"
print "\tLambda: %0.4f" % simulator.l
print "\tMu: %0.4f" % simulator.m
print "\tc: %0.4f" % simulator.c
print "\tK: %0.4f" % simulator.k

print "\tStability: %s (rho = %0.4f)" % (simulator.stability_condition(), simulator.rho)

print "\n+ TAIL"
print "\tMean number of clients (l) = %0.4f" % simulator.lm()
print "\tMean length (lq) = %0.4f" % simulator.lq()
print "\tMean time of a client waiting into the tail (w) = %0.4f" % simulator.w()

print "\n+ SYSTEM"
print "\tMean time of a client into the system (wq) = %0.4f" % simulator.wq()

print "\n+ PROBABILITY DSTRIUTION"
acum = 0
acum += simulator.p0()
print "\tP_0 = %0.10f" % simulator.p0()
for n in range(1, simulator.k + 1):
    a = simulator.pn(n)
    if n < 10 or n == simulator.k:
        print "\tP_%s = %0.10f" % (n, a)
    elif n == 10 and simulator.k != 11:
        print "\tP_%s = ..... " % n
        print "\t..... "
    acum += a
print "\t[Total Probability: %s]\n" % acum
nose.tools.assert_almost_equal(acum, 1.0)

print "Elapsed time: %0.8f" % (time.time() - t_start)















