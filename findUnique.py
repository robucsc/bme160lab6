#!/usr/bin/env python3
# Name: Rob Lodes (rlodes)
# Group Members: got help from various tutors, followed David's suggestions, and Jamie Moore (jmoore7)

"""
Read in FastA sequences via standard in (STDIN), find the unique subsequences that occur in each single
tRNA such that no members of this set occur among any of the other tRNA sets, then output in a text
graphical manner to illustrate the positioning of the essential elements.


usage: python findUniue.py < filename1 > output.txt
input: standard in (STDIN), in FastA format
output: A header for the sequence from the FastA, the sequence, and a text graphic of the essential elements

Example usage:
    python findUnique.py < dataFile.fa > output.txt
"""
# sys.stdout.reconfigure(encoding='utf-8')

class tRNA:
    '''
    Employs sets to find the essential segments in a sequence.
    Input is a file, or sequence via the FastA reader
    Output this class creates an outputList that can be used in the final output in main
    Methods: __init__, _buildPset, findUniques, findEssentials, findPosition
    '''
    tRNAlist = []
    def __init__(self, head, seq):
        '''
        Initialize the class, build the pSet, and set up the public sets and list
        input: head, seq
        '''
        self.head = head
        self.seq = seq
        self.pSet = self._buildPset()
        self.uniques = set()
        self.essentials = set()
        self.outputList = []
        tRNA.tRNAlist.append(self)

    def _buildPset(self):
        '''
        method called by inti to build the pSet
        output: myPset
        '''
        myPset = set()                                          # build the set
        for start in range(len(self.seq)):                      # outer loop find the start
            for end in range(start + 1, len(self.seq) + 1):     # find the end
                myPset.add(self.seq[start:end])                 # add to the set
        return myPset

    def findUniques(self):
        '''
        Add unique substrings to the pSet
        '''
        superPset = set()
        for current in tRNA.tRNAlist:           # travers the tRNAList
            if current is not self:             # check to see if self
                for substring in current.pSet:  # iterate over all substring of the object's powerSet
                    superPset.add(substring)    # add the substring to the superPset
        self.uniques = self.pSet - superPset    # remove the superPset from the pSet, and set to uniques


    def findEssentials(self):       # how to create a 'power set'
        '''
         Find and remove the non essential substrings
        '''
        nonEssentials = set()
        for element in self.uniques:                                        # traverse uniques
            if element[:-1] in self.uniques or element[1:] in self.uniques: # check the front and back to trim
                nonEssentials.add(element)                                  # add the element to non essentials
        self.essentials = self.uniques - nonEssentials                      # remove the non essentials

    def findPosition(self):
        '''
        Find the position of the substring, and add it to the outputList
        '''
        for element in self.essentials:                                     # traverse essentials
            position = self.seq.find(element)                               # hold the position of the element
            while(position != -1):                                          # while there are elements
                self.outputList.append(position * '.' + element)            # append outputList with element and dots
                position = self.seq.find(element, position + 1)             # update position
        self.outputList.sort(key = lambda x:len(x))                         # sort on length

import sys

class FastAreader:

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                if not line:  # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence

########################################################################
# Main
# Here is the main program
#
########################################################################

# def main(inCL=None):
def main(inFile=None):
    '''
    The Main function.
    Input: FastA sequences from STDIN
    Output: a list of essential substrings graphically positioned in relation to their position in the sequence
    Example use: python findUnique.py < dataFile.fa > output.txt
    '''

    myFasta = FastAreader(inFile)  # use this for debugging.
    # myFasta = FastAreader() # use this one for the command line (also to turn in)

    for head, seq, in myFasta.readFasta():
        seq = seq.replace('.', '').replace('_', '').replace('-', '')    # clean up seq
        tRNA(head, seq)                                                 # instanciate a tRNA
    tRNA.tRNAlist.sort(key = lambda element:element.head)               # sort base on the header
    for current in tRNA.tRNAlist:                                       # traverse the tRNAlist
        # send unique messages to all objects
        current.findUniques()                                           # call findUniques
        # send essential messages to all objects
        current.findEssentials()                                        # call findEssentials
        current.findPosition()                                          # call findPosition
        print(current.head.replace(' ', ''))                            # clean up header, and print
        print(current.seq)                                              # print the seq
        for element in current.outputList:                              # travers the outputList
            print(element)                                              # print the element

if __name__ == "__main__":
    # main()
    main(inFile='bos-tRNA-7.fa')   # use this for debugging