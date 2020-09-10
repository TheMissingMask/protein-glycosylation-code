#!/usr/bin/python


def linesFromFile(fName):

    '''just get the lines from a file and return them as a list'''

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
    
itpDict={
'[ atoms ]':[],
'[ bonds ]':[],
'[ constraints ]':[],
'[ angles ]':[],
'[ dihedrals ]':[],
'[ virtual_sitesn ]':[]
}

f=open('Protein.itp')
lines=f.readlines()
f.close()

for k in itpDict.keys():
    for l in linesBetween(lines,k,'['):
        itpDict[k].append(l)
        
glycanList=[]

f=open('topol.top')
lines=f.readlines()
f.close()

for l in lines:
    if 'Glycan' in l:
        glycanList.append('g'+l.split()[0][1:])

sidebeads=[]
for l in linesBetween(lines,'[ bonds ]','['):
    sidebeads.append(l.split()[0])

maxID=int(itpDict['[ atoms ]'][-1].split()[0])
maxRes=int(itpDict['[ atoms ]'][-1].split()[2])

for i in range(len(glycanList)):

    glycan=glycanList[i]
    sidebead=sidebeads[i]

    itpDict['[ bonds ]'].append(' %s  %s  %s  %s  %s\n'%(sidebead,maxID+1,2,0.52,1000))
    itpDict['[ angles ]'].append(' %s  %s  %s  %s  %s  %s\n'%(maxID+1,sidebead,str(int(sidebead)-1),2,152,500))

    f=open('%s.itp'%(glycan))
    lines=f.readlines()
    f.close()
    print(lines)

    for l in linesBetween(lines,'[ atoms ]','['):
        cols=l.split()
        newline=' %s  %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+maxID,cols[1],int(cols[2])+maxRes,cols[3],cols[4],int(cols[0])+maxID,cols[6])
        itpDict['[ atoms ]'].append(newline)
    
    for l in linesBetween(lines,'[ bonds ]','['):
        cols=l.split()
        newline=' %s  %s  %s  %s  %s\n'%(int(cols[0])+maxID,int(cols[1])+maxID,cols[2],cols[3],cols[4])
        itpDict['[ bonds ]'].append(newline)
    
    for l in linesBetween(lines,'[ constraints ]','['):
        cols=l.split()
        newline=' %s  %s  %s  %s \n'%(int(cols[0])+maxID,int(cols[1])+maxID,cols[2],cols[3])
        itpDict['[ constraints ]'].append(newline)
    
    for l in linesBetween(lines,'[ angles ]','['): 
        cols=l.split()
        newline=' %s  %s  %s  %s  %s  %s\n'%(int(cols[0])+maxID,int(cols[1])+maxID,int(cols[2])+maxID,cols[3],cols[4],cols[5])
        itpDict['[ angles ]'].append(newline)
    
    for l in linesBetween(lines,'[ virtual_sitesn ]','['):
        cols=l.split()
        newline=' %s  %s  %s  %s  %s\n'%(int(cols[0])+maxID,cols[1],int(cols[2])+maxID,int(cols[3])+maxID,int(cols[4])+maxID)
        itpDict['[ virtual_sitesn ]'].append(newline)
    
    maxID=int(itpDict['[ atoms ]'][-1].split()[0])
    maxRes=int(itpDict['[ atoms ]'][-1].split()[2])

f=open('new.itp','w')
f.write('[ moleculetype ]\nGlycoprotein 3\n\n')
for k in itpDict.keys():
    f.write('%s\n'%(k))
    for l in itpDict[k]:
        f.write(l)
    f.write('\n')
        
f.close()
