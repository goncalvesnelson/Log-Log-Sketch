"""
This module implements Hash Sketch probabilistic data structure which is able to calculate the cardinality of large multisets in a single pass using little auxiliary memory
"""
try:
    from bitarray import bitarray
except ImportError:
    raise ImportError('requires bitarray >= 0.3.4')

from hashlib import sha1
from math import ceil
from math import log
from math import pow
from sys import getsizeof
from Interface import Interface


def getSHA1Bin(word):
    """
    Returns the SHA1 hash binary representation of the word passed as argument
    Arguments:
    - `word`: argument to the SHA1 function
    """
    hash_s = sha1()
    hash_s.update(word)
    return bin(long(hash_s.hexdigest(),16))[2:].zfill(160)
    
    # returns the rank of the first 1 digit from the left
def getIndex(binString,endIndex=160):
    """
    Returns the possition of the first 1 bit from the left in the word until endIndex
    Arguments:
    - `binString`: String to be compared
    - `endIndex`: Maximum index tu use in the lookup for the 1 bit
    """
        #i = 0
        # hexInt = int(hexString,16)
        ## Is it worth compare if hexInt is still bigger than zero? 
        #while (i<endIndex and (hexInt & 0x01 != 1)): 
        #    hexInt = hexInt >> 1
        #    i+=1
        #return i
    res = -1
    try:
        res = binString.index('1')+1
    except(ValueError):
        res = endIndex
    return res


class Sketch(Interface):
    """
    This Class implements a Hash Sketch probabilistic data structure 
    """
    
    def __init__(self,maxCardinality):
        """Implementes a Hash Sketch
        maxCardinality
        this Sketch is able to count cardinalities up to cardinality *maxCardinality*
        """    
        self.__maxCardinality=maxCardinality
        self.__maxIndex = 0
        self.__bitarray = bitarray(160)
        self.__bitarray.setall(False)
        self.__name = "Sketch"

    def getName(self):
        return self.__name
        
    
    def add(self,item):
        """
        Adds the item to the Sketch
        Arguments:
        - `item`: item to be inserted
        """
        binword = getSHA1Bin(item)
        index = getIndex(binword)
        self.__bitarray[index]= True
        self.__maxIndex = max(self.__maxIndex,index)

    def getRawEstimate(self):
        """
        Returns the raw estimate of set's cardinality
        """
        
        return self.__maxIndex

    def getNumberEstimate(self):
       """Returns Linear Count estimate of the set's cardinality
       """
       zeroBits = self.__bitarray.count(False)
       fraction = float(zeroBits)/len(self.__bitarray)
       res = -len(self.__bitarray)*log(fraction)
       return res

class LogLogSketch(Interface):
    """Implements a LogLog Sketch
    """
    
    def __init__(self, maxCardinality, error_rate):
        """Implementes a LogLog Sketch
        *maxCardinality
        this Sketch is able to count cardinalities up to cardinality *maxCardinality*
        error_rate
        the error_rate of the sketch when calculating the cardinality of the set
        """ 
        if not (0 < error_rate < 1):
            raise ValueError("Error_Rate must be between 0 and 1.")
        if not maxCardinality > 0:
            raise ValueError("maxCardinality must be > 0")

        self._maxCardinality = maxCardinality
        #k     
        self._k = int(round(log(pow(1.30/error_rate,2),2)))
        # m = 2**k
        self._bucketNumber = 1<<self._k
        # Bucket size = loglog(MaxCardinality) bits
        self._bucketSize = self._wordSizeCalculator(self._maxCardinality)

        #M(1)... M(m) = 0
        
        self._bucketList =[bitarray(self._bucketSize) for _ in xrange(self._bucketNumber)]
        for barray in self._bucketList:
            barray.setall(False)

        self.__name = "LogLogSketch"
        
    def getName(self):
        return self.__name
        
    def add(self,item):
        """
        Adds the item to the LogLog Sketch
        Arguments:
        - `item`: Item to be added the the sketch
        """
        binword = getSHA1Bin(item)
        pos = int(binword[:self._k],2)
        #Retrives the position of leftmost 1 
        aux = getIndex(binword[self._k:],160-self._k)
        ##The position cannot be bigger than the maximum number that can be fitted in wordsize bits
        index = min(aux,(1<<self._bucketSize)-1)
        newValue = max(int(self._bucketList[pos].to01(),2),index)
        self._bucketList[pos] = bitarray(bin(newValue)[2:])
        #Perhaps it would be faster if operations were done in binary only
                                      
    def getNumberEstimate(self): 
        """
        Returns the estimate of the cardinality
        Arguments:
        """
        # E = am * m * 2**(1/m*sum M(j))
        m = self._bucketNumber
        e = 0.39701 * m*2**((1.0/m)*sum([int(x.to01(),2) for x in self._bucketList]))
        return e
        

    def __sizeof__(self):
        # return getsizeof(self._maxCardinality)+ getsizeof(self._k)+getsizeof(self._bucketNumber)+ \
        #        getsizeof(self._bucketList)
        return self._bucketNumber* self._bucketSize
    
    def _wordSizeCalculator(self,Nmax):
        """
        Estimates the size of the memory Units, using the maximum cardinality as an argument
        Arguments:
        - `Nmax`: Maximum Cardinality 
        """
        return int(ceil(log(log(Nmax,2),2)))
        
class SuperLogLogSketch(Interface):
    """Implements the improved version of LogLog Sketches, SuperLogLog Sketches
    """
    
    def __init__(self, maxCardinality, error_rate):
        """Implementes a SuperLogLog Sketch
        *maxCardinality
        this Sketch is able to count cardinalities up to cardinality *maxCardinality*
       error_rate
            the error_rate of the sketch when calculating the cardinality of the set
        """ 
        if not (0 < error_rate < 1):
            raise ValueError("Error_Rate must be between 0 and 1.")
        if not maxCardinality > 0:
            raise ValueError("maxCardinality must be > 0")

        self._maxCardinality = maxCardinality
        #k     
        self._k = int(round(log(pow(1.05/error_rate,2),2)))
        # m = 2**k
        self._bucketNumber = 1<<self._k
        # Bucket size = loglog(MaxCardinality/m +3) bits
        self._bucketSize = int(ceil(log(log(float(self._maxCardinality)/self._bucketNumber+3,2),2)))
        self._B = 1 << self._bucketSize
        self.__name = "SuperLogLogSketch"
        #M(1)... M(m) = 0
        self._bucketList = [0 for _ in xrange(self._bucketNumber)]

        
    def getName(self):
        return self.__name
        
    def getNumberEstimate(self,beta = 0.7):
        """
        Returns the estimate of the cardinality
        Arguments:
            beta= Used to get the truncated list. Keep only the *beta* smallest values and discard the rest
        """
        newList = sorted(self._bucketList)
        lastIndex = ceil(len(newList)*beta)
        nbeta = lastIndex/len(newList)
        newList = newList[:int(lastIndex)]
        # E = am * m * 2**(1/m*sum M(j))
        m = self._bucketNumber*nbeta
        e = 0.39701 * m*2**((1.0/m)*sum(newList))
        return e

    def _restrition_rule(self,unrestricted_value):
        return min(unrestricted_value,self._B)

    
    def add(self,item):
        """
        Adds the item to the LogLog Sketch
        Arguments:
        - `item`: Item to be added the the sketch
        """
        binword = getSHA1Bin(item)
        index = int(binword[:self._k],2)
        self._bucketList[index] = self._restrition_rule(max(self._bucketList[index],getIndex(binword[self._k:],160-self._k)))

class HyperLogLogSketch(Interface):
    """Implements a HyperLogLog Sketch
    """
    
    def __init__(self, maxCardinality, error_rate):
        """Implementes a HyperLogLog Sketch
        *maxCardinality
        this Sketch is able to count cardinalities up to cardinality *maxCardinality*
        error_rate
        the error_rate of the sketch when calculating the cardinality of the set
        """
        self.__ALPHA16=0.673
        self.__ALPHA32=0.697
        self.__ALPHA64=0.709
            
        if not (0 < error_rate < 1):
            raise ValueError("Error_Rate must be between 0 and 1.")
        if not maxCardinality > 0:
            raise ValueError("maxCardinality must be > 0")

        self._maxCardinality = maxCardinality
        #k     
        self._k = int(round(log(pow(1.04/error_rate,2),2)))
        # m = 2**k
        self._bucketNumber = 1<<self._k
        # Bucket size = loglog(MaxCardinality) bits
        self._bucketSize = self._wordSizeCalculator(self._maxCardinality)

        #M(1)... M(m) = 0
        
        self._bucketList =[0 for _ in xrange(self._bucketNumber)]
        
        self.__name = "HyperLogLogSketch"
        self._alpha = self.__getALPHA(self._bucketNumber)


    def __getALPHA(self,m):
        if m <=16:
            return self.__ALPHA16
        elif m<=32:      
            return self.__ALPHA32
        elif m<=64:      
            return self.__ALPHA64
        else:            
            return 0.7213/(1+1.079/float(m))
           
    def getName(self):
        return self.__name
        
    def add(self,item):
        """
        Adds the item to the LogLog Sketch
        Arguments:
        - `item`: Item to be added the the sketch
        """
        binword = getSHA1Bin(item)
        pos = int(binword[:self._k],2)
        #Retrives the position of leftmost 1 
        aux = getIndex(binword[self._k:],160-self._k)
        # Sets its own register value to maximum value seen so far
        self._bucketList[pos] = max(aux,self._bucketList[pos])
        

                                      
    def getNumberEstimate(self): 
        """
        Returns the estimate of the cardinality
        Arguments:
        """
        # E = am * m * 2**(1/m*sum M(j))
        m = self._bucketNumber
        raw_e = self._alpha*pow(m,2)/sum([pow(2,-x) for x in self._bucketList])
        if raw_e <= 5/2.0*m: # Small range correction
            #count number or registers equal to 0
            v = self._bucketList.count(0)
            if v!=0:
                return m*log(m/float(v),2)
            else:
                return raw_e
        elif raw_e <= 1/30.0*2**160: #intermidiate range correction -> No correction
            return raw_e
        else:
            return -2**160*log(1-raw_e/2.0**160,2)
        

    def join(self,*HyperLogLogList):
        """Joins the HyperLogLog Sketches passed as argument, with this HyperLogLog Sketch"""
        if HyperLogLogList:
            for sketch in HyperLogLogList:
                if type(sketch)!=type(self):
                    raise TypeError("all arguments must be HyperLogLog Sketches")           

            bucketLists = zip(self._bucketList,*[sketch._bucketList for sketch in HyperLogLogList])
            self._bucketList = map(max,bucketLists)
        else:
            raise TypeError("join expected at least 1 argument, got 0")
    
    def __sizeof__(self):
        return self._bucketNumber* self._bucketSize
    
    def _wordSizeCalculator(self,Nmax):
        """
        Estimates the size of the memory Units, using the maximum cardinality as an argument
        Arguments:
        - `Nmax`: Maximum Cardinality 
        """
        return int(ceil(log(log(Nmax,2),2)))    
 
def bs(number):
    """
    Returns the binary representation of the number given as argument

    - `number`: Number to be convertet to binary string
    """

    return str(number) if number<=1 else bs(number>>1) + str(number&1)


if __name__ == '__main__':
    import sys
    a = HyperLogLogSketch(2000000,0.05)
    b = HyperLogLogSketch(2000000,0.05)
    c = HyperLogLogSketch(2000000,0.05)
    
    for i in xrange(100000):
        a.add(str(i))
    for i in xrange(1500):
        b.add(str(i))
    for i in xrange(100000,200000):
        c.add(str(i))
    #print sys.getsizeof(a)
    print a.getNumberEstimate()
    print b.getNumberEstimate()
    print c.getNumberEstimate()
    c.join(a,b)
    print c.getNumberEstimate()
