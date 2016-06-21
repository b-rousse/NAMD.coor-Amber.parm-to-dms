#!/usr/bin/env python 
# -*- coding: utf-8 -*-
########################################################
#
# this script converts a snapshot from a namd trajectory
# to a Desmond restart file  for use on Anton. Please
# note that the script requires VMD version 1.9 or higher.
# It is based on the script convertNAMDtoDMS.py found on the Anton Wiki
# It was modified by Phil Blood during the 2015 workshop.
# Marcela Madrid added the boxa info assuming a cubic box, March 3, 2016.
#
# (C) 2010-2012 Markus Dittrich, NRBSC, PSC, CMU
#
# call with:
#
# vmd -dispdev text -python -e convertNAMDtoRST-final-system.py \
#      -args -p amberparmfile.top \
#      -c namd-system.restart.coor -v namd-system.restart.vel \
#      -x namd-system.xsc -o outfile -s -S "protein"
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
########################################################
#I converted psf to top file, so we can now use Amber-FF,NAMD-traj combinations.
#just view the variable "psfFile" as "topFile"
import sys
import optparse
from atomsel import *
from AtomSel import AtomSel
from Molecule import Molecule
from VMD import evaltcl


def parse_cmdline(cmdlineArgs):
    """
    This function initializes the command line parser.
    """

    parser = optparse.OptionParser("Usage: vmdt -python -e "
            "readVelocities.py -args [options]")

    parser.add_option("-p", "--psffile", action="store", dest="psfFile")
    parser.add_option("-c", "--coorfile", action="store", dest="coorFile")
    parser.add_option("-v", "--velfile", action="store", dest="velFile")
    parser.add_option("-x", "--xscfile", action="store", dest="xscFile")
    parser.add_option("-o", "--outputfile", action="store", dest="outFile")
    parser.add_option("-s", "--centerSystem", action="store_true", dest="doCenter")
    parser.add_option("-S", "--centerSelection", action="store", dest="centerSel")

    parser.set_defaults(doCenter = False, centerSel = "all")

    opts, args = parser.parse_args(cmdlineArgs)
    psfFile   = opts.psfFile
    coorFile  = opts.coorFile
    velFile   = opts.velFile
    xscFile   = opts.xscFile
    outFile   = opts.outFile
    doCenter  = opts.doCenter
    centerSel = opts.centerSel

    # all filenames are required
    if (psfFile == None) or (coorFile == None) or (velFile == None) \
       or (xscFile == None) or (outFile == None):

                parser.print_help()
                exit()


    return psfFile, coorFile, velFile, xscFile, outFile, doCenter, \
           centerSel



def load_velocities(psfFile, velFile):
    """
    Load the binary velocity file and extract velocities.
    """

    mol = Molecule()
    mol.load(psfFile, "parm7")
    mol.load(velFile, "namdbin")

    allVelocities = atomsel('all')
    xVel = allVelocities.get('x')
    yVel = allVelocities.get('y')
    zVel = allVelocities.get('z')

    # conversion from binvel units to A/ps
#    convFactor = 20.4582651391 
#    xVel = [v * convFactor for v in xVel]
#    yVel = [v * convFactor for v in yVel]
#    zVel = [v * convFactor for v in zVel]

    mol.delete()
    return xVel, yVel, zVel


def load_Coor(psfFile, coorFile):
    """
    Load the binary coor file and extract coors.
    """

    mol = Molecule()
    mol.load(psfFile, "parm7")
    mol.load(coorFile, "namdbin")

    allCoors = atomsel('all')
    xCor = allCoors.get('x')
    yCor = allCoors.get('y')
    zCor = allCoors.get('z')

    mol.delete()
    return xCor, yCor, zCor

def load_system(psfFile, coorFile):
    """
    Load the main system.
    """

    mol = Molecule()
    mol.load(psfFile, "parm7")						
    mol.load(coorFile, "namdbin")
    return mol



def set_velocities(mol, xVel, yVel, zVel):
    """
    Add the molecule velocities to the system.
    """

    allAtoms = atomsel("all")
    allAtoms.set("vx", xVel)
    allAtoms.set("vy", yVel)
    allAtoms.set("vz", zVel)

def set_Coor(mol, xCor, yCor, zCor):
    """
    Add the molecule Coors to the system.
    """

    allAtoms = atomsel("all")
    allAtoms.set("cx", xCor)
    allAtoms.set("cy", yCor)
    allAtoms.set("cz", zCor)

def save_mol_as_dms(mol, fileName):
    """
    Save the current molecule as dms file.
    """

    mol.save(fileName + ".rst7")



def set_pbc(xscFile):
    """
    Sets the systems periodic boundaries."
    """
    xscFile = open(xscFile,"r")
        
    for line in xscFile:
       continue

    items = line.split()
    xDim  = items[1]
    yDim  = items[5]
    zDim  = items[9]

    #set pbd
    pbcCommand = ("package require pbctools; pbc set { %s %s %s }" 
                % (xDim, yDim, zDim))
    evaltcl(pbcCommand)

    xscFile.close()
        

def load_pbc(xscFile):
    """
    Sets the systems periodic boundaries."
    """
    xscFile = open(xscFile,"r")

    for line in xscFile:
       continue

    items = line.split()
    xDim  = items[1]
    yDim  = items[5]
    zDim  = items[9]

    #set pbd
    #evaltcl(pbcCommand)

    xscFile.close()

    return xDim, yDim, zDim



def center_system(selection):
    """
    Center the system around the selection.
    """

    centerSel = atomsel(selection)
    center = centerSel.center()
    negCenter = [-1.0 * item for item in center]

    moveSel = atomsel("all")
    moveSel.moveby(negCenter)



def remove_tip3p_hh_bond():
    """
    This removes the bond between hydrogen atoms in
    TIP3P water if present since viparr will introduce
    the proper constraint.
    """

    # it looks like atomsel doesn't support set/getbonds
    # so we have to use the deprecated AtomSel for now
    oh2Sel = AtomSel("resname TIP3 and name OH2", 1)
    h1Sel  = AtomSel("resname TIP3 and name H1", 1)
  
    oh2Indices = oh2Sel.get("index")
    bondlist = []
    for i in oh2Indices:
        bondlist.append([i])
    h1Sel.setbonds(bondlist) 



#####################################################
# main routine
#####################################################
if __name__ == "__main__":

    # parse the command line
    psfFile, coorFile, velFile, xscFile, outfile, doCenter, \
        centerSel = parse_cmdline(sys.argv[1:])

    # transform NAMD to dms 
    x, y, z = load_Coor(psfFile, coorFile)
    vx, vy, vz = load_velocities(psfFile, velFile)
    box = load_pbc(xscFile)
    NATOM=#number of atoms in your system                                         #change for your system
    fo = open("output.rst",'w') 
    fo.write('\n')
    fo.write('%5i\n'% (NATOM))
    if NATOM %2==0:
        for i in range (0,NATOM-1,2):
            fo.write('%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n'% (x[i], y[i], z[i], x[i+1], y[i+1], z[i+1]))
        for i in range (0,NATOM-1,2):
            fo.write('%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n'% (vx[i],vy[i], vz[i], vx[i+1], vy[i+1], vz[i+1]))
    else:
        for i in range (0,NATOM-2,2):
            fo.write('%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n'% (x[i], y[i], z[i], x[i+1], y[i+1], z[i+1]))
        fo.write('%12.7f %12.7f %12.7f\n'% (x[NATOM-1], y[NATOM-1], z[NATOM-1]))
        for i in range (0,NATOM-2,2):
            fo.write('%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n'% (vx[i],vy[i], vz[i], vx[i+1], vy[i+1], vz[i+1]))
        fo.write('%12.7f %12.7f %12.7f\n'% (vx[NATOM-1], vy[NATOM-1], vz[NATOM-1]))
#The boxa part added by Marcela, assuming cubic box of water.
	boxa=90.0
    fo.write('%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f \n'% (float(box[0]),float(box[1]),float(box[2]),float(boxa),float(boxa),float(boxa)))
    fo.close()
 #   mol = load_system(psfFile, coorFile)

  #  if doCenter:
   #     center_system(centerSel)

   # set_velocities(mol, vx, vy, vz)
   # set_pbc(xscFile)
   # remove_tip3p_hh_bond()
   # save_mol_as_dms(mol, outfile)  
    
    exit()

