# Return a tuple indicating if Fixman potential or torque is present
def FPFTLabel(inpString):
    try:
        if inpString.find('NFPNFT') > 0:
            return (False, False)
        elif inpString.find('FPNFT') > 0:
            return (True, False)
        elif inpString.find('NFPFT') > 0:
            return (False, True)
        elif inpString.find('FPFT') > 0:
            return (True, True)
        else:
            raise Exception('FPFTLabel: string not found')
    except Exception as exception:
        print exception.args
#

#
def FPFTLabel2String(FPFTTuple):
    try:
        if FPFTTuple == (False, False):
            return 'NFPNFT'
        elif FPFTTuple == (True, False):
            return 'FPNFT'
        elif FPFTTuple == (False, True):
            return 'NFPFT'
        elif FPFTTuple == (True, True):
            return 'FPFT'
        else:
            raise TypeError
    except TypeError:
        print "FPFTLabel2String: need a tuple of two bools"
#
