## main_oop_lightPropMC.py
#
# monte carlo model of photon transport in infinite 
# homogeneous medium and calculation of fluence rate 
# for an isotropic point source
#
# prepared for  
#
# Proseminar Biophotonic - IQ/HOT; WiSe 2015/2016
#
# author: O. Melchert
# date:   2015-11-29
import sys
from random import expovariate,random, seed
from math import sqrt, log, pi


class FluenceRate(object):
        """class defining data structure for measurement
        """
        def __init__(self,nBins,nP):
                """initialize an instance of fluence rate data structure
                """
                self.nBins      = nBins 
                self.nPhotons   = nP
                self.shellWidth = 50
                self.h          = [0. for i in range(nBins)] 
                
        def pos2bin(self,pos):
                """convert packet position to bin id
                """
                x,y,z   = pos
                r       = sqrt(x*x+y*y+z*z)
                return min( int(r/self.shellWidth/1e-4),self.nBins-1)

        def measure(self,p):
                """accumulate absorbed photon weight in corresponding bin
                """
                self.h[self.pos2bin(p.pos)]+=p.dw

        def dump(self):
                """list fluence rate
                """
                fac = 4*pi*pow(self.shellWidth,3)*self.nPhotons*1e-12
                for i in range(self.nBins-1):
                        print i*self.shellWidth*1e-4, self.h[i]/(fac*(i*i+i+1./3.))



class PhotonPacket(object):
        """class defining a photon packet
        """
        def __init__(self,wgt=1., pos=[0.,0.,0.], dirCos=[0.,0.,1.]):
                """initialize an instance of a photon packet
                """
                self.wgt    = wgt
                self.pos    = pos 
                self.dirCos = dirCos
                self.dw     = 0.

        def advance(self,mt):
                """variable stepsize method to propagate photon packet 
                
                implements variable step size method to model 
                light propagation in medium, wherein probability
                density of packet stepsize s is given by

                p(s) = (\mu_a+\mu_s) \exp{-(\mu_a+\mu_s)s}
                
                Sampling s this way is equivalent to say that
                a photon is forced to be absorbed or scattered
                after each step with no interaction along the
                propagated distance (NOTE: small s have high p(s)),
                see Prahl et al., Proc SPIE IS 5 (1989) 102--111

                \param[in] mt Total attenuation 
                """
                # interaction-free propagation distance
                s = expovariate(mt)
                # compute new positions using stepsize and directional cosines
                self.pos = map(lambda x,dirCos: x+s*dirCos,self.pos,self.dirCos)

        def scatter(self):
                """isotropic scattering
                
                implements isotropic scattering of non-absorbed part
                of the split wave packet
                """
                x1 = x2 = 1.
                while(x1*x1+x2*x2 >= 1):
                        x1 = 2.*random() - 1.
                        x2 = 2.*random() - 1.
                
                t2 = x1*x1 + x2*x2
                fac = sqrt(1.-t2)

                self.dirCos = [ 2.*x1*fac, 2.*x2*fac, 1.-2.*t2] 

        def absorb(self,a):
                """absorption of photon packet

                After a propagation step, a photon packet is split
                into two parts. One part is scattered (see method scatter())
                the other part is absorbed. Therefore, the  weight of
                the photon packet after the ith scattering event is 
                decreased according to the Beer-Lambert law for turbid 
                media:
                w_{i+1} = w_i (1-frac{\mu_s}{\mu_s+\mu_a})
                where frac{\mu_s}{\mu_s+\mu_a} is the single particle albedo
                
                \param[in] a Single particle albedo
                """
                dw        = (1.-a)*self.wgt
                self.dw   = dw
                self.wgt -= dw

        def roulette(self,pm=0.1):
                """photon termination

                the roulette technique gives a photon a chance 
                of surviving with probabilty pm and weight self.wgt/pm 
                or else its weight is set to zero. 
                This means that for, say, pm=0.1, a photon packet is 
                terminated 9 out of 10 times, but once it is kept with a 
                weight that is increase by a factor of 10. Hence, packets 
                are terminated most of the time but energy is conserved 
                by occationally keeping packets with increased weight.
                This removes the photon in an unbiased fashion without 
                violating energy conservation and without propagating 
                a packet with negligible weight.

                \param[in] pm Packet survival probability
                """
                if self.wgt < 1e-3:
                   if random() < pm:
                        self.wgt /= pm
                   else: 
                        self.wgt = -1.

        def __str__(self):
                return "wgt = %lf,  pos = [%lf,%lf,%lf],  dirCos = [%lf,%lf,%lf]\n"\
                                %(self.wgt,
                                  self.pos[0],self.pos[1],self.pos[2],
                                  self.dirCos[0],self.dirCos[1],self.dirCos[2])
                


def main():
        ## simulation parameters
        s     = 100     # rng seed for MC simulation
        nP    = 10000   # number of photon packets
        ma    = 2.      # absorption coeff
        ms    = 20.     # scattering coeff
        nBins = 100     # number of bins for measurement

        ## derived quantities
        mt = ma + ms    # total attenuation
        a  = ms/mt      # single particle albedo

        ## initialize data structures for measurement 
        ## of absorbed power density
        h = FluenceRate(nBins,nP)

        seed(s)
        while(nP):
            p = PhotonPacket()
            while(p.wgt > 0.):
                p.advance(mt)
                p.scatter()
                p.absorb(a)
                h.measure(p)
                p.roulette()
            nP -= 1
        
        h.dump()


main()
# EOF: main_lightPropMC.py
