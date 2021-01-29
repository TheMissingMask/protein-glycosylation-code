#!/usr/bin/python

# imports #

import subprocess
import argparse
import scipy
import numpy as np
import re
import os
import random
import glob
import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder
from scipy.spatial.distance import cdist


# define stuff #

def unit(v):
    
    '''get the unit vector of v'''
    
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

def findSites(prot,noncanonical=False,cull=True):

    '''locate the potential N-linked glycosylation sites (in NxS/T[/C] sequons, x!=P; excluding adjacent pairs of sequons if cull is True), returning their resnums'''

    # allow for noncanonical glycosylation sites
    if noncanonical==False:
        triplet3=['T','S']
    else:
        triplet3=['T','S','C']

    # find the PNGs
    siteList=[]
    for i in range(len(prot.residues.sequence())-2):
        if prot.residues.sequence()[i]=='N' and prot.residues.sequence()[i+1]!='P' and prot.residues.sequence()[i+2] in triplet3:
            siteList.append(i+prot.residues.resnums[0])
    # now we are going to remove glycosites that are immediately adjacent to each other
    if cull==True:
        for i in range(len(siteList)-1):
            if siteList[i+1]-siteList[i]==1:
                siteList.pop(i)

    return siteList

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
                if len(l.split())<1 or l.split()[0]==';' or l.split()[0][0]==';':
                    next
                else:
                    if 1==0:#';' in l:
                        next
                    else:
                        outLines.append(l)
    else:
        for l in lines:
            if l=='\n':
               next
            elif l.strip()==startString:
                copy=True
            elif l.strip()==stopString:
              copy=False
            elif copy:
                if len(l.split())<1 or l.split()[0]==';' or l.split()[0][0]==';':
                    next
                else:
                    if 1==0:# ';' in l:
                        next
                    else:
                        outLines.append(l)

    return outLines

def isSurface(pdb,c=5.0):
    
    '''identify the surface-exposed residues based on cutoff c'''

    os.system('gmx editconf -f %s.pdb -o boxed.pdb -bt cubic -c -d 2.0'%(pdb))
    os.system('gmx solvate -cp boxed.pdb -cs wat.gro -o solvated.pdb -radius 0.24')

    u=MDAnalysis.Universe('solvated.pdb')
    prot=u.select_atoms('not name W')
    sol=u.select_atoms('name W')

    saDict={}
    for res in prot.residues:
        minD=np.min(cdist(res.atoms.positions,sol.atoms.positions))
        if minD<=c:
            saDict[res.resnum]=1
        else:
            saDict[res.resnum]=0

    os.system('rm boxed.pdb solvated.pdb')

    return saDict



# set stuff up #

# there are better ways to do this -- need to edit
aaDict={
        'ARG':'R',
        'HIS':'H',
        'LYS':'K',
        'ASP':'D',
        'GLU':'E',
        'SER':'S',
        'THR':'T',
        'ASN':'N',
        'GLN':'Q',
        'CYS':'C',
        'GLY':'G',
        'PRO':'P',
        'ALA':'A',
        'ILE':'I',
        'LEU':'L',
        'MET':'M',
        'PHE':'F',
        'TRP':'W',
        'TYR':'Y',
        'VAL':'V'
        }

# for the topology prediction
mapDict={
        'O':0,
        'i':1,
        'M':2,
        'o':3
        }

# get user input #

parser=argparse.ArgumentParser()
parser.add_argument('-p','--prot',default='Protein',help='prefix of the itp and pdb files for the CG protein [default: Protein]')
parser.add_argument('-t','--glycanType',default='N',help='Type of glycosylation (N-linked [N], mucin-type [O], glycosaminoglycan [GAG], monosaccharide [mono])')
parser.add_argument('-nt','--nType',default='complex',help='Options are complex, hybrid, and oligomannose')
parser.add_argument('-s','--sites',default=None,nargs='+',help='sites to be glycosylated; if not provided, they will be predicted from the sequence (providing sites is necessary for all types of glycosylation other than N-linked)')
parser.add_argument('-o','--outName',default='glycoprotein',help='name for output files [default: glycoprotein]')
parser.add_argument('--bilayer',action='store_true',default=False,help='is there a bilayer present in the pdb file? [default: no]')
parser.add_argument('--cull',action='store_true',default=True,help='only glycosylate one site where two are adjacent to each other? [default: yes]')
parser.add_argument('--noncanonical',action='store_true',default=False,help='include noncanonical sequons? [default: no]')
parser.add_argument('--glycosaminoglycan',default=None,help='type of glycosaminoglycan to add if type is glycosaminoglycan')
parser.add_argument('--mono',default=None,help='identify of the monosaccharide if the glycan type is a monosaccharide')
args=parser.parse_args()


# get the protein coordinates, topology etc. #

print('preparing to glycosylate')

u=MDAnalysis.Universe('%s.pdb'%(args.prot))
prot=u.select_atoms('protein')

chainIDs=list(set(prot.chainIDs))

if len(chainIDs)>1:
    for chain in chainIDs:
        prot=u.select_atoms('protein and segid %s'%(chain))

    if args.sites:
        sites=args.sites
    else:
        sites=list(set(findSites(prot,noncanonical=args.noncanonical,cull=args.cull)))
    siteList=sorted([int(x) for x in sites])
    sites=siteList

    # for TM proteins, make sure no sites are in the TM or intracellular regions
    if args.bilayer==True:

        seq=[]
        fasta='>protein sequence\n'
        for r in prot.residues.resnames:
            fasta+=aaDict[r]
            seq.append(aaDict[r])

        f=open('prot.fasta','w')
        for l in fasta:
            f.write(l)
        f.close()

        os.system('tmhmm -f prot.fasta -m TMHMM2.0.model -p')

        topology=[]
        txt=glob.glob('protein.annotation')
        f=open(txt[0])
        lines=f.readlines()
        f.close()
        for l in lines[1:]:
            for k in l.strip():
                topology.append(k)

        topologyDict={}
        for i in range(len(prot.residues)):
            topologyDict[prot.residues[i].resid]=topology[i]

        for site in sites:
            if topologyDict[site]=='i' or topologyDict[site]=='M':
                sites.remove(site)

    # make sure all sites are solvent-exposed
    saDict=isSurface(args.prot)

    for site in sites:
        if saDict[site]!=1:
            sites.remove(site)

    # check that N-linked only on Asn and Mucin-type only on Ser or Thr
    for site in sites:
        if args.glycanType=='N':
            for res in prot.residues:
                if res.resnum==site and res.resname !='ASN':
                    sites.remove(site)
        elif args.glycanType=='O':
            for res in prot.residues:
                if res.resnum==site and res.resname not in ['SER','THR']:
                    sites.remove(site)

    os.system('rm boxed.pdb solvated.pdb')

    # prepare the topology file #
    outTop=open('topol-%s.top'%(chain),'w')
    outTop.write(
        '''#include "martini_v3.0.3.itp"
        #include "%s.itp"
        #ifdef POSRES
        #include "posre.itp"
        #endif
    '''%(args.prot))

    for site in sites:
        outTop.write('#include "glycan-%s.itp\n'%(site))

    outTop.write(
    '''#include "martini_v3.0_phospholipids.itp"
    #include "martini_v3.0_solvents.itp"
    #include "martini_v3.0_ions.itp"

    [ system ]
    ; name
    Glycoprotein

    [ molecules ]
    ; name  number
    Protein 1
    ''')

    for site in sites:
        outTop.write('Glycan-%s 1\n'%(site))

    outTop.write('\n\n[ intermolecular_interactions ]\n[ bonds ]\n')

    # Onward!

    print('the glycosites are:\n')
    print(sites)

    glycanbead=prot.residues[-1].atoms.ids[-1]+1
    protRes1=prot.residues.resnums[0]

    # now we loop through the sites
    
    for site in sites:

        print('doing site %s'%(site))

        for r in prot.residues:
            if r.resid==site:
                sidebead=r.atoms.ids[-1]

        outTop.write('%s %s   6     0.52   1000\n'%(sidebead,glycanbead))
        outTop.write('[ angles ]\n%s %s %s   2  152   500\n'%(glycanbead,sidebead,int(sidebead)-1))

        # before anything else, we should build the glycan
        if args.glycanType=='N':
            subprocess.Popen('python3 build-glycans.py -t %s -nt %s -o glycan-%s'%(args.glycanType,args.nType,site),shell=True).wait()
        elif args.glycanType=='O':
            subprocess.Popen('python3 build-glycans.py -t %s -o glycan-%s'%(args.glycanType,site),shell=True).wait()
        elif args.glycanType=='mono':
            subprocess.Popen('python3 build-glycans.py -t %s --mono %s -o glycan-%s'%(args.glycanType,args.mono,site),shell=True).wait()
        elif args.glycanType=='glycosaminoglycan':
            subprocess.Popen('python3 build-glycans.py -t %s --gag %s -o glycan-%s'%(args.glycanType,args.gag,site),shell=True).wait()

        print('glycan-%s'%(site))
        for r in prot.residues:
            if r.resid==site:
                ID=r.atoms.ids[-1]
        i=1
        addID=ID
        for r in prot.residues:
            if r.resid==site:
                addRes=r.resnum
        v=MDAnalysis.Universe('glycan-%s.pdb'%(site))
        agCore=v.select_atoms('all')
        newAtoms=[]

        siteID=ID

        pdbCore=MDAnalysis.Universe('glycan-%s.pdb'%(site))
        agCore=pdbCore.select_atoms('all')

        coords1=agCore.atoms.positions[:2]
        coords2=agCore.atoms.positions
        for r in prot.residues:
            if r.resid==site:
                coords0=r.atoms.positions

        rSele=''
        for r in prot.residues:
            if r.resid==site:
                resi=r
        print(resi.resname)
        for r in prot.residues:
            d=(np.linalg.norm(r.atoms.center_of_geometry()-resi.atoms.center_of_geometry()))
            if d<=10:
                rSele+='%s '%(r.resnum)
        group=prot.select_atoms('resnum %s'%(rSele))

        v1=group.atoms.positions[0]-group.atoms.positions[1]
        v1/=unit(v1)
        newCoords=coords0+v1#(v1*1.5) ################### AHOY

        vector=group.atoms.positions[1]-group.atoms.positions[0]
        normalisedVector=vector/np.linalg.norm(vector)
        newCoords=np.asarray([newCoords[0]+(0.55*normalisedVector),newCoords[1]+(0.55*normalisedVector)])

        A=np.matrix(coords1[:2])
        B=np.matrix(newCoords[:2]) # transforming the glycan onto the side chain, based on newCoords
        R,t=rigidTransform(A,B)
        A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
        A2=A2.T
        C=np.matrix(coords2)
        C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
        C2=C2.T
        C2=np.asarray(C2,dtype=np.float64)
        C2=np.around(C2,decimals=3)
        agCore.atoms.positions=C2
        if np.sum(cdist(C2,prot.atoms.positions)<4.0) >=5:
            next

        f=open('tmp.pdb','w')

        for a in prot.atoms:
            f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',a.id,a.name,' ',a.resname,' ',a.resnum,' ',a.position[0],a.position[1],a.position[2],1.00,0.00,' ',' ')) # write everything before the addition

        for a in agCore.atoms:
            f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',a.id+ID,a.name,' ',a.resname,' ',a.resnum+addRes,' ',a.position[0],a.position[1],a.position[2],1.00,0.00,' ',' '))

        f.close()

        subprocess.Popen('gmx editconf -f tmp.pdb -o renumbered.pdb -resnr %s'%(protRes1),shell=True).wait()
        subprocess.Popen('mv renumbered.pdb tmp.pdb',shell=True).wait()
        
        # before we move onto the next site, check for steric clashes
        siteidx=sites.index(site)
        siteidx1=sites.index(site)+1
        saDict=isSurface('tmp')
        if saDict[sites[siteidx1]]!=1:
            sites.remove(sites[siteidx1])

        u=MDAnalysis.Universe('tmp.pdb')
        prot=u.select_atoms('all')
        glycanbead=prot.residues[-1].atoms.ids[-1]+1

    subprocess.Popen('mv tmp.* %s-%s.pdb'%(args.outName,chain),shell=True).wait()
    outTop.close()

    print('He is all my art now\n')

else:

    if args.sites:
        sites=args.sites
    else:
        sites=list(set(findSites(prot,noncanonical=args.noncanonical,cull=args.cull)))
    siteList=sorted([int(x) for x in sites])
    sites=siteList

    # if a bilayer is present, make sure we don't add glycans to TM or intracellular regions
    if args.bilayer==True:

        seq=[]
        fasta='>protein sequence\n'
        for r in prot.residues.resnames:
            fasta+=aaDict[r]
            seq.append(aaDict[r])

        f=open('prot.fasta','w')
        for l in fasta:
            f.write(l)
        f.close()

        os.system('tmhmm -f prot.fasta -m TMHMM2.0.model -p')

        topology=[]
        txt=glob.glob('protein.annotation')
        f=open(txt[0])
        lines=f.readlines()
        f.close()
        for l in lines[1:]:
            for k in l.strip():
                topology.append(k)

        topologyDict={}
        for i in range(len(prot.residues)):
            topologyDict[prot.residues[i].resid]=topology[i]

        for site in sites:
            if topologyDict[site]=='i' or topologyDict[site]=='M':
                sites.remove(site)

    # ensure the sites are solvent-exposed
    saDict=isSurface(args.prot)

    for site in sites:
        if saDict[site]!=1:
            sites.remove(site)

    # check that N-linked only on Asn and Mucin-type only on Ser or Thr
    for site in sites:
        if args.glycanType=='N':
            for res in prot.residues:
                if res.resnum==site and res.resname !='ASN':
                    sites.remove(site)
        elif args.glycanType=='O':
            for res in prot.residues:
                if res.resnum==site and res.resname not in ['SER','THR']:
                    sites.remove(site)


    outTop=open('topol.top','w')
    outTop.write(
    '''#include "martini_v3.0.3.itp"
#include "%s.itp"
    '''%(args.prot))

    for site in sites:
        outTop.write('#include "glycan-%s.itp\n'%(site))

    outTop.write(
    '''#include "martini_v3.0_phospholipids.itp"
#include "martini_v3.0_solvents.itp"
#include "martini_v3.0_ions.itp"

[ system ]
; name
Glycoprotein

[ molecules ]
; name  number
Protein 1
''')

    for site in sites:
        outTop.write('Glycan-%s 1\n'%(site))

    outTop.write('\n\n[ intermolecular_interactions ]\n[ bonds ]\n')

# Onward!

    print('the glycosites are:\n')
    print(sites)

    glycanbead=prot.residues[-1].atoms.ids[-1]+1
    protRes1=prot.residues.resnums[0]
#addition=0

    for site in sites:

        print('doing site %s'%(site))

#    site+=addition

    # lets set up for the topology
        for r in prot.residues:
            if r.resid==site:
                sidebead=r.atoms.ids[-1]

        for tmpRes in prot.residues:
            if tmpRes.resnum in [site-1,site+1,site]:
                next
            else:
                for r in prot.residues:
                    if r.resid==site:
                        separation=np.linalg.norm(cdist(tmpRes.atoms.positions,r.atoms.positions))
                if separation<=10.0:
                    clash=True
                else:
                    clash=False
        if clash==True:
            next

        outTop.write('%s %s   6     0.52   1000\n'%(sidebead,glycanbead)) # 2 or 6?
        outTop.write('[ angles ]\n%s %s %s   2  152   500\n'%(glycanbead,sidebead,int(sidebead)-1))

        # before anything else, we should build the glycan
        if args.glycanType=='N':
            subprocess.Popen('python3 build-glycans.py -t %s -nt %s -o glycan-%s'%(args.glycanType,args.nType,site),shell=True).wait()
        elif args.glycanType=='O':
            subprocess.Popen('python3 build-glycans.py -t %s -o glycan-%s'%(args.glycanType,site),shell=True).wait()
        elif args.glycanType=='mono':
            subprocess.Popen('python3 build-glycans.py -t %s --mono %s -o glycan-%s'%(args.glycanType,args.mono,site),shell=True).wait()
        print('glycan-%s'%(site))
        for r in prot.residues:
            if r.resid==site:
                ID=r.atoms.ids[-1]
        i=1
        addID=ID
        for r in prot.residues:
            if r.resid==site:
                addRes=r.resnum
        v=MDAnalysis.Universe('glycan-%s.pdb'%(site))
        agCore=v.select_atoms('all')
        newAtoms=[]

        siteID=ID

        pdbCore=MDAnalysis.Universe('glycan-%s.pdb'%(site))
        agCore=pdbCore.select_atoms('all')

        coords1=agCore.atoms.positions[:2]
        coords2=agCore.atoms.positions
        for r in prot.residues:
            if r.resid==site:
                coords0=r.atoms.positions

        rSele=''
        for r in prot.residues:
            if r.resid==site:
                resi=r
        print(resi.resname)
        for r in prot.residues:
            d=(np.linalg.norm(r.atoms.center_of_geometry()-resi.atoms.center_of_geometry()))
            if d<=10:
                rSele+='%s '%(r.resnum)
        group=prot.select_atoms('resnum %s'%(rSele))

        v1=group.atoms.positions[0]-group.atoms.positions[1]
        v1/=unit(v1)
        newCoords=coords0+v1

        vector=group.atoms.positions[1]-group.atoms.positions[0]
        normalisedVector=vector/np.linalg.norm(vector)
        newCoords=np.asarray([newCoords[0]+(0.55*normalisedVector),newCoords[1]+(0.55*normalisedVector)])

        A=np.matrix(coords1[:2])
        B=np.matrix(newCoords[:2]) # transforming the glycan onto the side chain, based on newCoords
        R,t=rigidTransform(A,B)
        A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
        A2=A2.T
        C=np.matrix(coords2)
        C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
        C2=C2.T
        C2=np.asarray(C2,dtype=np.float64)
        C2=np.around(C2,decimals=3)
        agCore.atoms.positions=C2
        if np.sum(cdist(C2,prot.atoms.positions)<4.0) >=5:
            next

        f=open('tmp.pdb','w')

        for a in prot.atoms:
            f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',a.id,a.name,' ',a.resname,' ',a.resnum,' ',a.position[0],a.position[1],a.position[2],1.00,0.00,' ',' ')) # write everything before the addition

        for a in agCore.atoms:
            f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',a.id+ID,a.name,' ',a.resname,' ',a.resnum+addRes,' ',a.position[0],a.position[1],a.position[2],1.00,0.00,' ',' '))

        f.close()

        subprocess.Popen('gmx editconf -f tmp.pdb -o renumbered.pdb -resnr %s'%(protRes1),shell=True).wait()
        subprocess.Popen('mv renumbered.pdb tmp.pdb',shell=True).wait()
        
        siteidx=sites.index(site)
        siteidx1=sites.index(site)+1
        saDict=isSurface('tmp')
        if saDict[sites[siteidx1]]!=1:
            sites.remove(sites[siteidx1])
            
        u=MDAnalysis.Universe('tmp.pdb')
        prot=u.select_atoms('all')
        glycanbead=prot.residues[-1].atoms.ids[-1]+1

    subprocess.Popen('mv tmp.* %s.pdb'%(args.outName),shell=True).wait()
    outTop.close()

print('He is all my art now\n')
