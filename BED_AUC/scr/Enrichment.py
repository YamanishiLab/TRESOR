import math
import os
import string
import sys

class Enrichment:
   
    def __init__(self):
        self.ranks = []
        self.nbr_actives = 0
        self.nbr_total = 0
    
    def getNbrInactives(self): 
        return self.nbr_total - self.nbr_actives
    
    def getNbrActives(self):
        return self.nbr_actives
    
    def getNbrTotal(self):
        return self.nbr_total
    
    def readFromFile(self, file):
        self.ranks = []
        self.nbr_actives = -1
        self.nbr_total = -1
        if not os.access( file, os.R_OK ):
            sys.stderr.writelines( "Cannot access file '" + file + "'.\n" )
            sys.exit(1)
        lLines = open( file ).readlines()
        iRank = -1
        iCount = -1
        iCountPrev = -1
        for sLine in lLines:
            # sLine = string.strip( sLine )
            sLine = sLine.strip( )
            if len( sLine ) == 0:
                continue
            # lWords = string.split( sLine )
            lWords = sLine.split()
            iRank = int( lWords[0] ) # 予測された順位
            iCount = int( lWords[1] ) # 正例内での順位
            if iCount == 0:
                sys.stderr.write("In Enrichment.readFromFile(), the cummulative count cannot be zero.")
                sys.exit(1)
            if len( self.ranks ) > 0 and iRank < self.ranks[-1]:
                sys.stderr.write("In Enrichment.readFromFile(), the rank must be increasing going down.")
                sys.exit(1)
            if iRank <= 0:
                sys.stderr.write("In Enrichment.readFromFile(), the rank must be positive.")
                sys.exit(1)
            if iCount < iCountPrev:
                sys.stderr.write("In Enrichment.readFromFile(), the cummulative count must be increasing.")
            iCountPrev = iCount
            self.ranks.append( iRank )
        self.nbr_actives = iCount
        self.nbr_total = self.ranks[-1]
        if self.nbr_actives < self.nbr_total:
            pass
        del self.ranks[-1]

    def readListData(self, y_rank, y_count, total): # This method was added by Satoko Namba.
        # If you discover some misstakes, please contact the following e-mail address: namba.satoko775@mail.kyutech.jp.
        self.ranks = []
        self.nbr_actives = -1
        self.nbr_total = total
        iRank = -1
        iCount = -1
        iCountPrev = -1
        for iRank, iCount in zip(y_rank, y_count):
            if iCount == 0:
                sys.stderr.write("In Enrichment.readListData(), the cummulative count cannot be zero.")
                sys.exit(1)
            if len( self.ranks ) > 0 and iRank < self.ranks[-1]:
                sys.stderr.write("In Enrichment.readListData(), the rank must be increasing going down.")
                sys.exit(1)
            if iRank <= 0:
                sys.stderr.write("In Enrichment.readListData(), the rank must be positive.")
                sys.exit(1)
            if iCount < iCountPrev:
                sys.stderr.write("In Enrichment.readListData(), the cummulative count must be increasing.")
            iCountPrev = iCount
            self.ranks.append( iRank )
        self.nbr_actives = iCount
        if self.nbr_actives < self.nbr_total:
            pass

    def toString(self):
        sReturn = ""
        iCount = 1
        for iRank in self.ranks:
            sReturn += "%d %d\n" % (iRank, iCount)
            iCount += 1
        sReturn += "%d %d\n" % ( self.getNbrTotal(), self.getNbrActives() )
        return sReturn
    
    def calculateBEDROC(self, alpha = 20.0 ):
        if alpha < 0.00001:
            os.stderr.write( "In method calculatBEDROC, the alpha parameter argument must be greater than zero." )
            sys.exit(1)
        N = float( self.getNbrTotal() )
        n = float( self.getNbrActives() )
        sum = 0.0
        for rank in self.ranks:
            sum += math.exp( -alpha * rank / N )
        ra = n/N   
        factor1 = ra * math.sinh( alpha/2.0 )/( math.cosh(alpha/2.0) - math.cosh(alpha/2.0 - ra*alpha ) )
        factor2 = 1.0 / ra * (math.exp(alpha/N) - 1.0)/( 1.0 - math.exp(-alpha))
        constant = 1.0 / ( 1.0 - math.exp( alpha * ( 1.0 - ra ) ) )
        bedroc = sum * factor1 * factor2 + constant
        return bedroc
    
