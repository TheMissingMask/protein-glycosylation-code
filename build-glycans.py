#!/usr/bin/python


# imports #

import os
import argparse
import numpy as np
import re
import random
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import random
import MDAnalysis


# functions #

def unit(v):
    
    '''get a unit vector of v'''
    
    x=v[0]
    y=v[1]
    z=v[2]

    if x!=0:
        x=float(x)/np.linalg.norm(x)
    if y!=0:
        y=float(y)/np.linalg.norm(y)
    if z!=0:
        z=float(z)/np.linalg.norm(z)

    return np.divide(v,np.asarray([x,y,z]))

def rigidTransform(A,B):

    '''map coordinates of A onto those of B'''

    assert len(A)==len(B)

    # total number of points
    N=A.shape[0];

    # find the centroids i.e. the geometric mean of each set of coordinates
    cenA=np.mean(A,axis=0)
    cenB=np.mean(B,axis=0)

    # centre the points to remove the translation component
    AA=A-np.tile(cenA,(N,1))
    BB=B-np.tile(cenB,(N,1))

    # find the covariance matrix H
    H=np.transpose(AA)*BB

    # SVD will decompose a matrix into three other matrices such that [U,S,V]=SVD(E) and E=USV**T
    U,S,Vt=np.linalg.svd(H)

    R=Vt.T*U.T

    # special reflection case
    if np.linalg.det(R)<0:
        print('Reflection detected')
        Vt[2,:]*=-1
        R=Vt.T*U.T

    # find the translation t
    t=-R*cenA.T+cenB.T

    # return what we have found
    return R,t

def linesFromFile(fName):

    '''get the lines from a file and return them as a list'''

    f=open(fName,'r')
    lines=f.readlines()
    f.close()

    return lines

def linesBetween(lines,startString,stopString):

    '''collect a list of the lines between startString and stopString'''

    outLines=[]
    copy=False

    if stopString=='':
        for l in lines:
            if l=='\n':
                next
            elif l.strip()==startString:
                copy=True
            elif copy:
                if len(l.split())<1 or l.split()[0]==';':
                    next
                else:
                    outLines.append(l)
    else:
        for l in lines:
            if l=='\n':
                next
            elif l.strip()==startString:
                copy=True
            elif l.strip()==stopString or stopString in l:
                copy=False
            elif copy:
                if len(l.split())<1 or l.split()[0]==';':
                    next
                else:
                    outLines.append(l)

    return outLines


def updateITP(itpList,bondDict,angleDict,mapDict,bond0,VS0,pair,newID,newResnum,itpLines):

    '''update the itp file with the latest glycan residue'''

    # first, let us work out what is going on
    rRes=pair.split('-')[0]
    nrRes=pair.split('-')[2]
    bond1=pair.split('-')[1]
    bondAngle=bond1

    bondPair='%s-%s'%(mapDict[rRes],mapDict[nrRes])

    newLines=itpLines

    # atoms
    newAtoms=linesBetween(newLines,'[ atoms ]','[')
    tmp=[]
    for l in newAtoms:
        cols=l.split()
        tmp.append('%s  %s  %s  %s  %s  %s  %s\n'%(newID[cols[0]],cols[1],newResnum,cols[3],cols[4],newID[cols[0]],cols[6]))
    VS1=int(newID[cols[0]])
    VS0=int(VS0)
    newAtoms=tmp

    # bonds
    newBonds=linesBetween(newLines,'[ bonds ]','[ constraints ]')
    tmp=[]
    for l in newBonds:
        cols=l.split()
        tmp.append('%s  %s  %s  %s  %s\n'%(newID[cols[0]],newID[cols[1]],cols[2],cols[3],cols[4]))
    newBonds=tmp
    if bondPair in bondDict.keys():
        newBonds.append('%s  %s  2  %s  %s\n'%(VS1,VS0,bondDict[bondPair][0],bondDict[bondPair][1]))
    else:
        newBonds.append('%s  %s  2  %s  %s\n'%(VS1,VS0,bondDict[None][0],bondDict[None][1]))
   
    # constraints
    newConstraints=linesBetween(newLines,'[ constraints ]','[ angles ]')
    tmp=[]
    for l in newConstraints:
        cols=l.split()
        tmp.append('%s  %s  1  %s\n'%(newID[cols[0]],newID[cols[1]],cols[3]))
    newConstraints=tmp

    # angles
    newAngles=linesBetween(newLines,'[ angles ]','[ dihedrals ]')
    tmp=[]
    for l in newAngles:
        cols=l.split()
        tmp.append('%s  %s  %s  2  %s  %s\n'%(newID[cols[0]],newID[cols[1]],newID[cols[2]],cols[4],cols[5]))
    newAngles=tmp
    countingAtoms0=linesBetween(
            linesFromFile('%s.itp'%(nrRes)),
            '[ atoms ]',
            '[ bonds ]')[-1].split()[0]
    ### as for the previous monosaccharide, we need to consider something else entirely
    countingAtoms1=linesBetween(
            linesFromFile('%s.itp'%(rRes)),
            '[ atoms ]',
            '[ bonds ]')[-1].split()[0]
    ### now we can make the angles
    if bondAngle in angleDict.keys():
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(int(VS0)-(int(countingAtoms0))+1,VS0,VS1,angleDict[bondAngle][0][0],angleDict[bondAngle][0][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0-(int(countingAtoms0))+2,VS0,VS1,angleDict[bondAngle][1][0],angleDict[bondAngle][1][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0-(int(countingAtoms0))+3,VS0,VS1,angleDict[bondAngle][2][0],angleDict[bondAngle][2][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0,VS1,VS1-(int(countingAtoms1))+1,angleDict[bondAngle][3][0],angleDict[bondAngle][3][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0,VS1,VS1-(int(countingAtoms1))+2,angleDict[bondAngle][4][0],angleDict[bondAngle][4][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0,VS1,VS1-(int(countingAtoms1))+3,angleDict[bondAngle][5][0],angleDict[bondAngle][5][1]))

    # dihedrals -- at the moment we have none of these
    newDihedrals=linesBetween(newLines,'[ dihedrals ]','[ exclusions ]')

    # exclusions -- at the moment we have none of these
    newExclusions=linesBetween(newLines,'[ exclusions ]','[ virtual_sitesn ]')

    newVirtuals=linesBetween(newLines,'[ virtual_sitesn ]','')
    tmp=[]
    for l in newVirtuals:
        cols=l.split()
        tmp.append('%s  1  %s  %s  %s\n'%(VS1,newID['1'],newID['2'],newID['3']))
    newVirtuals=tmp

    itpList[0]=itertools.chain(itpList[0],newAtoms)
    itpList[1]=itertools.chain(itpList[1],newBonds)
    itpList[2]=itertools.chain(itpList[2],newConstraints)
    itpList[3]=itertools.chain(itpList[3],newAngles)
    itpList[4]=itertools.chain(itpList[4],newDihedrals)
    itpList[5]=itertools.chain(itpList[5],newExclusions)
    itpList[6]=itertools.chain(itpList[6],newVirtuals)

    return itpList,VS1,VS1

def loadPDB(filename):

    '''get the PDB coordinates and complete lines'''

    lines=linesFromFile(filename)

    coords=[]
    for l in lines:
        if l.split()[0]=='ATOM':
            coords.append(l)

    lineDict={}
    coordDict={}

    for c in coords:
        cols=c.split()
        lineDict[int(cols[4])]=[]
        coordDict[int(cols[4])]=[]

    for c in coords:
        cols=c.split()
        lineDict[int(cols[4])].append(c)
        coordDict[int(cols[4])].append([float(x) for x in cols[5:8]])

    return lineDict,coordDict

def updatePDB(glycanCoords,rRes,resnum0,maxResnum,maxID,pdbLines):

    '''update the coordinates of the glycan'''

    pos1=np.asarray(glycanCoords[1][-1])
    pos2=np.asarray(glycanCoords[2][-1])
    #if len(glycanCoords.keys())>=3:
    #    pos3=np.asarray(glycanCoords[3][-1])
    #else:
    #    pos3=np.asarray(pos2+(pos2-pos1))
    pos3=np.asarray(pos2+(pos2-pos1))

    translateList=np.asarray([
        np.subtract(pos2,pos1),
        np.subtract(pos3,pos2)])
    translate1=np.asarray([np.mean(translateList[:,0]),np.mean(translateList[:,1]),np.mean(translateList[:,2])])
    position=np.asarray(glycanCoords[int(resnum0)][-1])

    resLines=linesFromFile('%s.pdb'%(rRes))
    VS=np.asarray([float(x) for x in resLines[-1].split()[5:8]])
    if rRes=='Fuc':
        translate1=np.asarray([translate1[1],translate1[2],translate1[0]])
    position=np.add(translate1,position)
    translate2=np.subtract(position,VS)

    pdbLines[maxResnum+1]=[]
    glycanCoords[maxResnum+1]=[]

    for l in resLines:
        maxID+=1
        cols=l.split()
        coords=np.asarray([float(x) for x in cols[5:8]])
        newCoords=np.add(translate2,coords)
        pdbLines[maxResnum+1].append('ATOM     %3d  %s %s    %2s     %.3f %.3f %.3f  1.00  1.00\n'%(maxID,cols[2],cols[3],maxResnum+1,newCoords[0],newCoords[1],newCoords[2]))
        glycanCoords[maxResnum+1].append([newCoords[0],newCoords[1],newCoords[2]])

    return pdbLines,glycanCoords

def theMaestro(itpList,pair,bond0,VS0,maxResnum,maxID,resnum0,pdbLines,glycanCoords):

    '''orchestrates much of the rest of the glycosylation script'''

    print('I lost my shoe')

    rRes=re.split('\d',pair)[0][:-1] # this should get the residue to be added, and remove the 'a' or 'b' from the end
    bond1='%s%s'%(re.split('\d',pair)[0][-1],re.split('[A-z]',pair)[0]) # this should get the 'a'/'b' and the bond number
    nrRes=re.split('\d',pair)[1] # this should get the residue already present in the chain, with no numbers or anomeric identifiers

    # itp stuff
    itpLines=linesFromFile('%s.itp'%(rRes)) # get the topology information for the residue being added, which is named the same as in the library
    newResnum=int(maxResnum)+1
    newID={}
    for a in linesBetween(itpLines,'[ atoms ]','[ bonds ]'):
            newID[a.split()[0]]=int(a.split()[0])+maxID
    itpList,maxID,VS0=updateITP(itpList,bondDict,angleDict,mapDict,bond0,VS0,pair,newID,newResnum,itpLines)

    # pdb stuff
    pdbLines,glycanCoords=updatePDB(glycanCoords,rRes,resnum0,maxResnum,maxID,pdbLines)

    maxResnum=newResnum

    return itpList,maxResnum,maxID,pdbLines,glycanCoords,VS0


# set up the data dictionaries #

# interresidue bonds between the virtual sites
bondDict={
        'DHX-a2-DHX':0.47,
        'DHX-a3-DHX':0.51,
        'HXN-a4-HXA':0.47,
        'HXS-a1-HXS':0.49,
        'HXS-a2-HXS':0.52,
        'HXS-a3-HXS':0.52,
        'HXS-a4-HXS':0.49,
        'HXS-a6-HXS':0.58,
        'HXS-b1-HXS':0.56,
        'HXS-b2-HXS':0.55,
        'HXS-b3-HXS':0.56,
        'HXS-b4-HXS':0.54,
        'HXS-b6-HXS':0.64,
        'NHX-a3-HXS':0.52,
        'NHX-b1-HXS':0.58,
        'NHX-b2-HXS':0.57,
        'NHX-b3-HXS':0.56,
        'NHX-b6-HXS':0.65,
        'SIA-a3-NHX':0.51,
        'SIA-a4-NHX':0.46,
        'SIA-a6-NHX':0.57,
        'SIA-b6-NHX':0.52,
        'NHX-b4-SIA':0.53,
        'HXA-a3-NHX':0.51,
        'HXA-a4-HXN':0.50,
        'HXA-b4-HXN':0.54,
        'HXS-a3-DHX':0.51,
        'HXS-b4-DHX':0.49,
        'HXS-a1-NHX':0.51,
        'HXS-a3-NHX':0.53,
        'HXS-b3-NHX':0.55,
        'HXS-b6-NHX':0.62,
        'HXA-a1-HXS':0.51,
        'NHX-a3-NHX':0.53,
        'NHX-a4-NHX':0.49,
        'NHX-b3-NHX':0.56,
        'NHX-b4-NHX':0.53,
        'XYL-b4-XYL':0.57,
        'HXS-b4-XYL':0.51,
        'DHX-a3-NHX':0.53,
        'DHX-b3-NHX':0.57,
        'DHX-a6-HXS':0.63,
        'SIA-b6-HXS':0.51,
        'HXN-a4-HXN':0.50,
        'HXN-b4-HXN':0.53,
        'HXN-b6-HXN':0.63,
        'NHX-b3-DHX':0.56,
        'SIA-a8-SIA':0.65,
        'XXX-a1-XXX':0.50,
        'XXX-a2-XXX':0.50,
        'XXX-a3-XXX':0.52,
        'XXX-a4-XXX':0.49,
        'XXX-a5-XXX':0.54, 
        'XXX-a6-XXX':0.59,
        'XXX-b1-XXX':0.57,
        'XXX-b2-XXX':0.56,
        'XXX-b3-XXX':0.56,
        'XXX-b4-XXX':0.53,
        'XXX-b5-XXX':0.55,
        'XXX-b6-XXX':0.60,
        'XXX-b8-XXX':0.65
        }

# interresidue bonds, between virtual sites and ring beads
angleDict={
        'a3':[(85,60),(31,60),(145,50),(61,60),(134,50),(78,50)],
        'a4':[(142,40),(46,60),(82,50),(61,60),(124,50),(85,50)],
        'b3':[(102,30),(24,30),(129,30),(52,60),(145,60),(77,50)],
        'a6':[(128,10),(119,1),(16,30),(64,10),(122,20),(85,30)],
        'b6':[(131,1),(113,1),(15,20),(53,20),(132,30),(87,20)],
        'a2':[(45,40),(93,50),(130,20),(71,40),(125,40),(78,30)],
        'b4':[(125,50),(52,30),(88,40),(53,60),(139,60),(82,50)],
        'b2':[(38,20),(85,30),(144,30),(43,50),(138,50),(82,50)],
        'b1':[(47,50),(150,30),(86,30),(47,50),(139,40),(83,40)],
        'a1':[(36,60),(130,50),(99,50),(49,70),(121,50),(97,50)],
        'a8':[(76,1),(124,1),(48,1),(72,1),(102,1),(107,1)]
        }

# mapping of atomistic monosaccharides to Martini monosaccharides
mapDict={
        'GlcA':'HXA',
        'IdoA':'HXA',
        'Glc':'HXS',
        'Man':'HXS',
        'Gal':'HXS',
        'All':'HXS',
        'Alt':'HXS',
        'GlcNAc':'NHX',
        'GalNAc':'NHX',
        'GlcN':'HXN',
        'Neu5Ac':'SIA',
        'Neu5Gc':'SIA',
        'Fuc':'DHX',
        'Rha':'DHX',
        'Qui':'DHX',
        'Xyl':'XYL',
        'Gal3S':'HXS3S',
        'GalNAc3S':'NHX3S',
        'GlcA1S':'HXA1S',
        'IdoA1S':'HXA1S',
        'GlcNAc3S':'NHX3S',
        'GlcNAc2S':'NHX2S',
        'GalNAc2S':'NHX2S',
        'GalNAc4S':'NHX4S',
        'GlcNAc4S':'NHX4S'
        }

# the conditional probability dictionary --  built using glytoucan data
f=open('conP.dat')
lines=f.readlines()
f.close()

pairDict={}
for l in lines:
    cols=l.split()
    pairDict[cols[0]]=float(cols[1])

    
# get user input #

parser=argparse.ArgumentParser()
parser.add_argument('-t','--glycanType',default='N',help='Type of glycosylation (N-linked [N], mucin-type [O], glycosaminoglycan [GAG], monosaccharide [mono], domain-specific [D])')
parser.add_argument('-nt','--nType',default='complex',help='Options are complex, hybrid, and oligomannose')
parser.add_argument('--glycosaminoglycan',default=None,help='Options are ...')
parser.add_argument('--mono',default=None)
parser.add_argument('--domain',default=None) # need to include options (e.g. EGF, TSR, collagen), and corresponding bits of code, with residue search function
parser.add_argument('-o','--outName',default='glycan')
args=parser.parse_args()


# N-glycosylation #

if args.glycanType=='N':
    
    def checkN(pair):
        rRes=re.split('[\?|a|b][\d|\?]',pair)[0]
        if rRes not in ['Man','Gal','GlcNAc','GalNAc','Fuc','Neu']:
            return False
        else:
            return True

    capList=['Neu','a2','a3','a6','a4','Fuc']
    addList=['Gal','GalNAc','GlcNAc','Fuc','Neu']

    f=open('conP.dat')
    lines=f.readlines()
    f.close()

    pairDict={}
    for l in lines:
        cols=l.split()
        pairDict[cols[0]]=float(cols[1])

    core='GlcNAcb2Mana6(GlcNAcb2Mana3)Manb4GlcNAcb4GlcNAc'
    base=core.split(')')[1]

    newGlycans=[]

    branches=re.findall('\(.*?\)',core)
    branches.append(re.sub('\(.*?\)','',core))

    for branch in branches:

        branch=re.sub('[\)|\(]','',branch)
        tmp=branch
        nrRes=re.split('[a|b]\d',branch)[1]
        rRes=re.split('[a|b]\d',branch)[0]
        bond0=re.findall('[a|b]\d',branch)[0]

        counter=0
        broken=False
        while counter<=3:
            nrRes=rRes
            for pair in pairDict.keys():
                if checkN(pair)==False:
                    next
                if re.split('[\?|a|b][\d|\?]',pair)[1]==nrRes:
                    p=pairDict[pair]
                    if p>=random.random():
                        rRes=re.split('[\?|a|b][\d|\?]',pair)[0]
                        bond0=re.findall('[\?|a|b][\d|\?]',pair)[0]
                        if '?' in bond0:
                            next
                        elif rRes not in addList:
                            next
                        else:
                            tmp=('%s%s'%(rRes,bond0))+tmp
                            counter+=1
                            if rRes in capList or bond0 in capList:
                                broken=True
                                break
                            elif counter>=3:
                                broken=True
                                break
            if broken:
                break
        newGlycans.append(tmp)
    cpx=True
    if cpx:
        glycan=re.sub(base,'',tmp)+'['+newGlycans[0]+']'+base
    elif hybrid:
        glycan='Mana6[Mana3]Mana6[%s]'%(newGlycans[0])+base
    elif man:
        glycan=random.choice([
            'Mana6[Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana2Mana6[Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana6[Mana2Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana6[Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana2Mana6[Mana2Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana2Mana6[Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana6[Mana2Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana2Mana6[Mana2Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana2Mana6[Mana2Mana3]Mana6[Mana2Mana2Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana6[Mana2Mana3]Mana6[Mana2Mana2Mana3]Manb4GlcNAcb4GlcNAc',
            'Mana2Mana6[Mana3]Mana6[Mana2Mana2Mana3]Manb4GlcNAcb4GlcNAc'
            ])

    glycan=re.sub('Neu','Neu5Ac',glycan)
    a=glycan
    numbers=re.compile(r'([a|b]\d)')
    a=numbers.sub(r'(\1)',a)
    f=open('%s.txt'%(args.outName),'w')
    f.write(a)
    f.close()
    glycanToAdd=re.sub('[\(|\)]','',a)
    glycanToAdd=re.sub('\[','(',glycanToAdd)
    glycanToAdd=re.sub('\]',')',glycanToAdd)
    ####test-N.py ends here####
    
    #how do we use this string to build the glycan?
    #test case example is:   Galb4GlcNAcb2Mana6(Galb4GlcNAcb2Mana3)Manb4GlcNAcb4GlcNAc
    preBranch=glycanToAdd.split(')')[-1]
    branches=[
        glycanToAdd.split('(')[0],
        glycanToAdd.split(')')[0].split('(')[1]
              ]
    preBranchResidues=re.split('[a|b]\d',preBranch)
    preBranchBonds=re.findall('[a|b]\d',preBranch)
    preBranchResidues.reverse()
    preBranchBonds.reverse()

    itpDict={
     'atoms':[], # list of strings to print directly out into the file
     'bonds':[],
     'constraints':[],
     'angles':[],
     'dihedrals':[],
     'virtual_sitesn':[]
    }
    
    atomN=0
    resN=0
    newatoms=0
    beadNumbers={}
    for i in range(len(preBranchResidues)):
        beadNumbers[i]=[]
    i=0
    for res in preBranchResidues:
        cgRes=mapDict[res]
        resLines=linesFromFile('%s.itp'%(cgRes))
        coreEndNumbers=[]
        # atoms
        for l in linesBetween(resLines,'[ atoms ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+resN,cols[3],cols[4],int(cols[5])+atomN,cols[6])
            itpDict['atoms'].append(newLine)
            beadNumbers[i].append(int(cols[0])+atomN) # so now we have a list of all the atom ids of that residue
            if int(cols[0])<=3:
                coreEndNumbers.append(int(cols[0])+atomN)
        newatoms=int(cols[0])
        # constraints
        for l in linesBetween(resLines,'[ constraints ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,cols[2],cols[3])
            itpDict['constraints'].append(newLine)
        # angles
        for l in linesBetween(resLines,'[ angles ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,int(cols[2])+atomN,cols[3],cols[4],cols[5])
            itpDict['angles'].append(newLine)
        # virtual sites
        for l in linesBetween(resLines,'[ virtual_sitesn ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+atomN,int(cols[3])+atomN,int(cols[4])+atomN)
            itpDict['virtual_sitesn'].append(newLine)
        atomN+=newatoms
        i+=1
        resN+=1
        # now we need to add the bonds between the residues
    for i in range(len(preBranchBonds)):
        bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
        bondType='%s-%s-%s'%(mapDict[preBranchResidues[i]],preBranchBonds[i],mapDict[preBranchResidues[i+1]])
        if bondType not in bondDict.keys():
            bondLength=bondDict['XXX-%s-XXX'%(preBranchBonds[i])]
        else:
            bondLength=bondDict[bondType]
        newLine=' %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],2,bondLength,1250)
        itpDict['bonds'].append(newLine)
    # now we need to add the angle
    for i in range(len(preBranchBonds)):
        bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
        bondBeads=[beadNumbers[i][0],beadNumbers[i][1],beadNumbers[i][2],beadNumbers[i+1][0],beadNumbers[i+1][1],beadNumbers[i+1][2]]
        bondType='%s'%(preBranchBonds[i])
        angles=angleDict[bondType]
        newLines=[]
        for j in range(3):
            newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondBeads[j],bondVS[0],bondVS[1],2,angles[j][0],angles[j][1]))
        for j in range(3,6):
            newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],bondBeads[j],2,angles[j][0],angles[j][1]))
        for l in newLines:
            itpDict['angles'].append(l)
    coreEndAtom=atomN # this will be linked to each branch
    coreEndRes=cgRes

    for branch in branches:

        branchResidues=re.split('[a|b]\d',branch)[:-1]
        bondToCore=re.findall('[a|b]\d',branch)[-1]
        branchBonds=re.findall('[a|b]\d',branch)[:-1]
        branchResidues.reverse()
        branchBonds.reverse()

        newatoms=0
        beadNumbers={}
        for i in range(len(branchResidues)):
            beadNumbers[i]=[]
        i=0
        for res in branchResidues:
            cgRes=mapDict[res]
            resLines=linesFromFile('%s.itp'%(cgRes))
            # atoms
            for l in linesBetween(resLines,'[ atoms ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+resN,cols[3],cols[4],int(cols[5])+atomN,cols[6])
                itpDict['atoms'].append(newLine)
                beadNumbers[i].append(int(cols[0])+atomN) # so now we have a list of all the atom ids of that residue
            newatoms=int(cols[0])
            # constraints
            for l in linesBetween(resLines,'[ constraints ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,cols[2],cols[3])
                itpDict['constraints'].append(newLine)
            # angles
            for l in linesBetween(resLines,'[ angles ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,int(cols[2])+atomN,cols[3],cols[4],cols[5])
                itpDict['angles'].append(newLine)
            # virtual sites
            for l in linesBetween(resLines,'[ virtual_sitesn ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+atomN,int(cols[3])+atomN,int(cols[4])+atomN)
                itpDict['virtual_sitesn'].append(newLine)
            atomN+=newatoms
            i+=1
            resN+=1
        # now we need to add the bonds between the residues
        for i in range(len(branchBonds)):
            bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
            bondType='%s-%s-%s'%(mapDict[branchResidues[i]],branchBonds[i],mapDict[branchResidues[i+1]])
            if bondType not in bondDict.keys():
                bondLength=bondDict['XXX-%s-XXX'%(branchBonds[i])]
            else:
                bondLength=bondDict[bondType]
            newLine=' %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],2,bondLength,1250)
            itpDict['bonds'].append(newLine)
        # now we need to add the angle
        for i in range(len(branchBonds)):
            bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
            bondBeads=[beadNumbers[i][0],beadNumbers[i][1],beadNumbers[i][2],beadNumbers[i+1][0],beadNumbers[i+1][1],beadNumbers[i+1][2]]
            bondType='%s'%(branchBonds[i])
            angles=angleDict[bondType]
            newLines=[]
            for j in range(3):
                newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondBeads[j],bondVS[0],bondVS[1],2,angles[j][0],angles[j][1]))
            for j in range(3,6):
                newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],bondBeads[j],2,angles[j][0],angles[j][1]))
            for l in newLines:
                itpDict['angles'].append(l)
        branchEndAtom=atomN # this will be linked to each branch
        branchEndRes=cgRes
        bondType='%s-%s-%s'%(branchEndRes,bondToCore,coreEndRes) # correct way round?
        if bondType not in bondDict.keys():
            bondLength=bondDict['XXX-%s-XXX'%(branchBonds[i])]
        else:
            bondLength=bondDict[bondType]
        newLine=' %s  %s  %s  %s  %s\n'%(branchEndAtom,coreEndAtom,2,bondLength,1250)
        itpDict['bonds'].append(newLine)
        if branch==branches[0]:
            branch0Res=branchResidues
        elif branch==branches[1]:
            branch1Res=branchResidues
            newLines=[]
            angles=angleDict[bondType.split('-')[1]]
            for j in range(3):
                newLines.append(' %s  %s  %s  %s  %s  %s\n'%(beadNumbers[0][j],branchEndAtom,coreEndAtom,2,angles[j][0],angles[j][1]))
            for j in range(3,6):
                newLines.append(' %s  %s  %s  %s  %s  %s\n'%(branchEndAtom,coreEndAtom,coreEndNumbers[j-3],2,angles[j][0],angles[j][1]))

    coords=[]
    for i in range(len(coreEndRes)):
        coords.append(np.array([5,i+50,0]))
    coordCount=0
    res=preBranchResidues[0]
    newLines=[]
    lines=linesFromFile('%s.pdb'%(mapDict[res]))
    for l in lines:
        newLines.append(l)
    u=MDAnalysis.Universe('%s.pdb'%(mapDict[res]))
    ag=u.select_atoms('all')
    resCoords=[]
    k=2
    for j in range(len(preBranchResidues[1:])):

        res=preBranchResidues[j+1]
        v=MDAnalysis.Universe('%s.pdb'%(mapDict[res]))
        newRes=v.select_atoms('all')
        coords1=newRes.atoms.positions[:3]
        coords2=newRes.atoms.positions

        coords0=ag.atoms.positions[:3]
        group=ag
        v1=group.atoms.positions[0]-group.atoms.positions[1]
        v1/=unit(v1)
        newCoords=coords0+v1
        vector=group.atoms.positions[1]-group.atoms.positions[0]
        normalisedVector=vector/np.linalg.norm(vector)
        newCoords=np.asarray([newCoords[0]+(9*normalisedVector),newCoords[1]+(9*normalisedVector)])

        A=np.matrix(coords1[:2])
        B=np.matrix(newCoords[:2])
        R,t=rigidTransform(A,B)
        A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
        A2=A2.T
        C=np.matrix(coords2)
        C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
        C2=C2.T
        C2=np.asarray(C2,dtype=np.float64)
        C2=np.around(C2,decimals=3)

        resCoords=C2
        res=preBranchResidues[j+1]
        for i in range(len(resCoords)):
            newLines.append('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',i+1,'XXX','',mapDict[res],'',k,'',resCoords[i][0],resCoords[i][1],resCoords[i][2],1.00,1.00,'',''))

        f=open('test.pdb','w')
        for l in newLines:
            f.write(l)
        f.close()

        u=MDAnalysis.Universe('test.pdb')
        ag=u.select_atoms('resid %s'%(k))
        k+=1

    # now, test.pdb gives the basis for each of the branches
    for branch in branches:
        if branch==branches[1]:
            minus=True
        else:
            minus=False
        
        u=MDAnalysis.Universe('test.pdb')
        ag=u.select_atoms('resid 3')
        resCoords=[]

        branchResidues=re.split('[a|b]\d',branch)[:-1]
        branchResidues.reverse()

        k=4
        for j in range(len(branchResidues)):

            res=branchResidues[j] # need to define this
            v=MDAnalysis.Universe('%s.pdb'%(mapDict[res]))
            newRes=v.select_atoms('all')
            coords1=newRes.atoms.positions[:3]
            coords2=newRes.atoms.positions

            coords0=ag.atoms.positions[:3]+np.array([5,5,5])
            newCoords=coords0+v1
            vector=group.atoms.positions[1]-group.atoms.positions[0]
            normalisedVector=vector/np.linalg.norm(vector)
            if minus:
                newCoords=np.asarray([newCoords[0]+(-7*normalisedVector),newCoords[1]+(-7*normalisedVector)])
            else:
                newCoords=np.asarray([newCoords[0]+(7*normalisedVector),newCoords[1]+(7*normalisedVector)])

            A=np.matrix(coords1[:2])
            B=np.matrix(newCoords[:2])
            R,t=rigidTransform(A,B)
            A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
            A2=A2.T
            C=np.matrix(coords2)
            C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
            C2=C2.T
            C2=np.asarray(C2,dtype=np.float64)
            C2=np.around(C2,decimals=3)

            resCoords=C2
            res=branchResidues[j]
            for i in range(len(resCoords)):
                newLines.append('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',i+1,'XXX','',mapDict[res],'',k,'',resCoords[i][0],resCoords[i][1],resCoords[i][2],1.00,1.00,'',''))


            f=open('test1.pdb','w')
            for l in newLines:
                f.write(l)
            f.close()

            u=MDAnalysis.Universe('test1.pdb')
            ag=u.select_atoms('resid %s'%(k))
            k+=1

    f=open('test.itp','w')
    f.write('[ moleculetype ]\n%s  3\n\n'%(args.outName))
    for k in itpDict.keys():
        f.write('[ %s ]\n'%(k))
        for l in itpDict[k]:
            f.write(l)
        f.write('\n')
    f.close()

    os.system('gmx editconf -f test1.pdb -o %s.pdb -resnr 1'%(args.outName))
    os.system('mv test.itp %s.itp'%(args.outName))

    
# mucin-type O-linked #

elif args.glycanType=='O':

    capList=['Neu','Fuc']
    addList=['Gala3','Galb3','Galb4','GlcNAcb6','GlcNAcb3','GlcNAca6']

    f=open('conP.dat')
    lines=f.readlines()
    f.close()

    pairDict={}
    for l in lines:
        cols=l.split()
        pairDict[cols[0]]=float(cols[1])

    core=np.random.choice(['Galb3GalNAc','Galb3(GlcNAcb6)GalNAc','GlcNAcb3GalNAc','GlcNAcb3(GlcNAcb6)GalNAc','GalNAca3GalNAc','GalNAcb6GalNAc','GalNAca6GalNAc','Gala3GalNAc'],p=[0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1])
    base='GalNAc'

    newGlycans=[]

    branches=re.findall('\(.*?\)',core)
    branches.append(re.sub('\(.*?\)','',core))
    if core=='Galb3GalNAc':
        branches=[core]

    for branch in branches:
        branch=re.sub('[\)|\(]','',branch)
        tmp=branch
        nrRes=re.split('[a|b]\d',branch)[1]
        rRes=re.split('[a|b]\d',branch)[0]
        bond0=re.findall('[a|b]\d',branch)[0]

        counter=0
        broken=False
        while counter<=3:
            nrRes=rRes
            for pair in pairDict.keys():
                if re.split('[\?|a|b][\d|\?]',pair)[1]==nrRes:
                    p=pairDict[pair]
                    if p>=random.random():
                        rRes=re.split('[\?|a|b][\d|\?]',pair)[0]
                        bond0=re.findall('[\?|a|b][\d|\?]',pair)[0]
                        if '%s%s'%(rRes,bond0) not in addList and rRes not in ['Fuc','Neu']:
                            next
                        else:
                            tmp=('%s%s'%(rRes,bond0))+tmp
                            counter+=1
                            if rRes in capList:
                                broken=True
                                break
            if broken:
                break
        newGlycans.append(tmp)

    if len(newGlycans)>1:
        glycan=re.sub(base,'',tmp)+'['+newGlycans[0]+']'+base
    else:
        glycan=newGlycans[0]
    glycan=re.sub('Neu','Neu5Ac',glycan)
    a=glycan
    numbers=re.compile(r'([a|b]\d)')
    a=numbers.sub(r'(\1)',a)
    f=open('%s.txt'%(args.outName),'w')
    f.write(a)
    f.close()
    glycanToAdd=re.sub('[\(|\)]','',a)
    glycanToAdd=re.sub('\[','(',glycanToAdd)
    glycanToAdd=re.sub('\]',')',glycanToAdd)
    preBranch=glycanToAdd.split(')')[-1]
    if len(glycanToAdd.split('('))>1:
        branches=[
            glycanToAdd.split('(')[0],
            glycanToAdd.split(')')[0].split('(')[1]
                  ]
    else:
        branches=[]
    preBranchResidues=re.split('[a|b]\d',preBranch)
    preBranchBonds=re.findall('[a|b]\d',preBranch)
    preBranchResidues.reverse()
    preBranchBonds.reverse()

    itpDict={
     'atoms':[],
     'bonds':[],
     'constraints':[],
     'angles':[],
     'dihedrals':[],
     'virtual_sitesn':[]
    }
    
    atomN=0
    resN=0
    newatoms=0
    beadNumbers={}
    for i in range(len(preBranchResidues)):
        beadNumbers[i]=[]
    i=0
    k=0
    for res in preBranchResidues:
        k+=1
        cgRes=mapDict[res]
        resLines=linesFromFile('%s.itp'%(cgRes))
        # atoms
        for l in linesBetween(resLines,'[ atoms ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+resN,cols[3],cols[4],int(cols[5])+atomN,cols[6])
            itpDict['atoms'].append(newLine)
            beadNumbers[i].append(int(cols[0])+atomN) # so now we have a list of all the atom ids of that residue
        newatoms=int(cols[0])
        # constraints
        for l in linesBetween(resLines,'[ constraints ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,cols[2],cols[3])
            itpDict['constraints'].append(newLine)
        # angles
        for l in linesBetween(resLines,'[ angles ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,int(cols[2])+atomN,cols[3],cols[4],cols[5])
            itpDict['angles'].append(newLine)
        # virtual sites
        for l in linesBetween(resLines,'[ virtual_sitesn ]','['):
            cols=l.split()
            newLine=' %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+atomN,int(cols[3])+atomN,int(cols[4])+atomN)
            itpDict['virtual_sitesn'].append(newLine)
        atomN+=newatoms
        i+=1
        resN+=1
        # now we need to add the bonds between the residues
    for i in range(len(preBranchBonds)):
        bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
        bondType='%s-%s-%s'%(mapDict[preBranchResidues[i]],preBranchBonds[i],mapDict[preBranchResidues[i+1]])
        if bondType not in bondDict.keys():
            bondLength=bondDict['XXX-%s-XXX'%(preBranchBonds[i])]
        else:
            bondLength=bondDict[bondType]
        newLine=' %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],2,bondLength,1250)
        itpDict['bonds'].append(newLine)
    # now we need to add the angle
    for i in range(len(preBranchBonds)):
        bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
        bondBeads=[beadNumbers[i][0],beadNumbers[i][1],beadNumbers[i][2],beadNumbers[i+1][0],beadNumbers[i+1][1],beadNumbers[i+1][2]]
        bondType='%s'%(preBranchBonds[i])
        angles=angleDict[bondType]
        newLines=[]
        for j in range(3):
            newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondBeads[j],bondVS[0],bondVS[1],2,angles[j][0],angles[j][1]))
        for j in range(3,6):
            newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],bondBeads[j],2,angles[j][0],angles[j][1]))
        for l in newLines:
            itpDict['angles'].append(l)
    coreEndAtom=atomN # this will be linked to each branch
    coreEndRes=cgRes

    for branch in branches:

        branchResidues=re.split('[a|b]\d',branch)[:-1]
        bondToCore=re.findall('[a|b]\d',branch)[-1]
        branchBonds=re.findall('[a|b]\d',branch)[:-1]
        branchResidues.reverse()
        branchBonds.reverse()

        newatoms=0
        beadNumbers={}
        for i in range(len(branchResidues)):
            beadNumbers[i]=[]
        i=0
        for res in branchResidues:
            cgRes=mapDict[res]
            resLines=linesFromFile('%s.itp'%(cgRes))
            # atoms
            for l in linesBetween(resLines,'[ atoms ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+resN,cols[3],cols[4],int(cols[5])+atomN,cols[6])
                itpDict['atoms'].append(newLine)
                beadNumbers[i].append(int(cols[0])+atomN) # so now we have a list of all the atom ids of that residue
            newatoms=int(cols[0])
            # constraints
            for l in linesBetween(resLines,'[ constraints ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,cols[2],cols[3])
                itpDict['constraints'].append(newLine)
            # angles
            for l in linesBetween(resLines,'[ angles ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,int(cols[1])+atomN,int(cols[2])+atomN,cols[3],cols[4],cols[5])
                itpDict['angles'].append(newLine)
            # virtual sites
            for l in linesBetween(resLines,'[ virtual_sitesn ]','['):
                cols=l.split()
                newLine=' %s  %s  %s  %s  %s\n'%(int(cols[0])+atomN,cols[1],int(cols[2])+atomN,int(cols[3])+atomN,int(cols[4])+atomN)
                itpDict['virtual_sitesn'].append(newLine)
            atomN+=newatoms
            i+=1
            resN+=1

        for i in range(len(branchBonds)):
            bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
            bondType='%s-%s-%s'%(mapDict[branchResidues[i]],branchBonds[i],mapDict[branchResidues[i+1]])
            if bondType not in bondDict.keys():
                bondLength=bondDict['XXX-%s-XXX'%(branchBonds[i])]
            else:
                bondLength=bondDict[bondType]
            newLine=' %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],2,bondLength,1250)
            itpDict['bonds'].append(newLine)

        for i in range(len(branchBonds)):
            bondVS=(beadNumbers[i][-1],beadNumbers[i+1][-1])
            bondBeads=[beadNumbers[i][0],beadNumbers[i][1],beadNumbers[i][2],beadNumbers[i+1][0],beadNumbers[i+1][1],beadNumbers[i+1][2]]
            bondType='%s'%(branchBonds[i])
            angles=angleDict[bondType]
            newLines=[]
            for j in range(3):
                newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondBeads[j],bondVS[0],bondVS[1],2,angles[j][0],angles[j][1]))
            for j in range(3,6):
                newLines.append(' %s  %s  %s  %s  %s  %s\n'%(bondVS[0],bondVS[1],bondBeads[j],2,angles[j][0],angles[j][1]))
            for l in newLines:
                itpDict['angles'].append(l)
        branchEndAtom=atomN # this will be linked to each branch
        branchEndRes=cgRes
        bondType='%s-%s-%s'%(branchEndRes,bondToCore,coreEndRes) # correct way round?
        if bondType not in bondDict.keys():
            bondLength=bondDict['XXX-%s-XXX'%(branchBonds[i])]
        else:
            bondLength=bondDict[bondType]
        newLine=' %s  %s  %s  %s  %s\n'%(branchEndAtom,coreEndAtom,2,bondLength,1250)
        itpDict['bonds'].append(newLine)
        if branch==branches[0]:
            branch0Res=branchResidues
        elif branch==branches[1]:
            branch1Res=branchResidues

    coords=[]
    for i in range(len(coreEndRes)):
        coords.append(np.array([5,i+50,0]))

    coordCount=0
    res=preBranchResidues[0]
    newLines=[]
    lines=linesFromFile('%s.pdb'%(mapDict[res]))
    for l in lines:
        newLines.append(l)
    u=MDAnalysis.Universe('%s.pdb'%(mapDict[res]))
    ag=u.select_atoms('all')
    resCoords=[]
    for j in range(len(preBranchResidues[1:])):

        res=preBranchResidues[j+1]
        v=MDAnalysis.Universe('%s.pdb'%(mapDict[res]))
        newRes=v.select_atoms('all')
        coords1=newRes.atoms.positions[:3]
        coords2=newRes.atoms.positions

        coords0=ag.atoms.positions[:3]
        group=ag
        v1=group.atoms.positions[0]-group.atoms.positions[1]
        v1/=unit(v1)
        newCoords=coords0+v1
        vector=group.atoms.positions[1]-group.atoms.positions[0]
        normalisedVector=vector/np.linalg.norm(vector)
        newCoords=np.asarray([newCoords[0]+(9*normalisedVector),newCoords[1]+(9*normalisedVector)])

        A=np.matrix(coords1[:2])
        B=np.matrix(newCoords[:2])
        R,t=rigidTransform(A,B)
        A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
        A2=A2.T
        C=np.matrix(coords2)
        C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
        C2=C2.T
        C2=np.asarray(C2,dtype=np.float64)
        C2=np.around(C2,decimals=3)

        resCoords=C2
        res=preBranchResidues[j+1]
        for i in range(len(resCoords)):
            newLines.append('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',i+1,'XXX','',mapDict[res],'',k,'',resCoords[i][0],resCoords[i][1],resCoords[i][2],1.00,1.00,'',''))

        f=open('test.pdb','w')
        for l in newLines:
            f.write(l)
        f.close()

        u=MDAnalysis.Universe('test.pdb')
        ag=u.select_atoms('resid %s'%(k))
        k+=1

    # now, test.pdb gives the basis for each of the branches
    for branch in branches:
        if branch==branches[1]:
            minus=True
        else:
            minus=False
        
        u=MDAnalysis.Universe('test.pdb')
        ag=u.select_atoms('resid %s'%(k))
        resCoords=[]

        ###
        branchResidues=re.split('[a|b]\d',branch)[:-1]
        branchResidues.reverse()

        for j in range(len(branchResidues)):

            res=branchResidues[j] # need to define this
            v=MDAnalysis.Universe('%s.pdb'%(mapDict[res]))
            newRes=v.select_atoms('all')
            coords1=newRes.atoms.positions[:3]
            coords2=newRes.atoms.positions

            coords0=ag.atoms.positions[:3]+np.array([5,5,5])
            newCoords=coords0+v1
            vector=group.atoms.positions[1]-group.atoms.positions[0]
            normalisedVector=vector/np.linalg.norm(vector)
            if minus:
                newCoords=np.asarray([newCoords[0]+(-7*normalisedVector),newCoords[1]+(-7*normalisedVector)])
            else:
                newCoords=np.asarray([newCoords[0]+(7*normalisedVector),newCoords[1]+(7*normalisedVector)])

            A=np.matrix(coords1[:2])
            B=np.matrix(newCoords[:2])
            R,t=rigidTransform(A,B)
            A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
            A2=A2.T
            C=np.matrix(coords2)
            C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
            C2=C2.T
            C2=np.asarray(C2,dtype=np.float64)
            C2=np.around(C2,decimals=3)

            resCoords=C2
            res=branchResidues[j]
            for i in range(len(resCoords)):
                newLines.append('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',i+1,'XXX','',mapDict[res],'',k,'',resCoords[i][0],resCoords[i][1],resCoords[i][2],1.00,1.00,'',''))

            f=open('test1.pdb','w')
            for l in newLines:
                f.write(l)
            f.close()

            u=MDAnalysis.Universe('test1.pdb')
            ag=u.select_atoms('resid %s'%(k))
            k+=1

    f=open('%s.itp'%(args.outName),'w')
    f.write('[ moleculetype ]\n%s  3\n\n'%(args.outName))
    for k in itpDict.keys():
        f.write('[ %s ]\n'%(k))
        for l in itpDict[k]:
            f.write(l)
        f.write('\n')
    f.close()

    os.system('gmx editconf -f test1.pdb -o %s.pdb -resnr 1'%(args.outName))

    
# glycosaminoglycans #

############### HERE ##################

# monosaccharide #

elif args.glycanType=='mono':

    monoToAdd=args.mono

    os.system('gmx editconf -f %s.pdb -o %s.pdb -resnr 1'%(mapDict[monoToAdd],args.outName))
    os.system('cp %s.itp %s.itp'%(mapDict[monoToAdd],args.outName))

print('An artist should create beautiful things, but should put nothing of his own life into them')
