#!/usr/local/bin/python

# Name: AutoDoge (AUTODock Output Generator / Extractor)
# Author: xywu@PKUHSC wuxingyu@bjmu.edu.cn
# Version: 0.11 (3/23/2016)
# Written in Python 2.7.6
# Description: Generate AutoDock4 VS report and extract best-binding conformations
#              for best-binding or best-clustering cluster.
#              For more options and how to run this script, please have a look
#              at the readme file provided within the package.
# Operating system supported: Windows / Linux

import glob
import os
import sys
from sys import platform as _platform
import argparse
from argparse import RawTextHelpFormatter
import timeit
import subprocess
import fileinput

########################Open Debug Mode########################
DEBUG = 0
###############################################################

# Argparser:
def argParser():
    print " NOTE: -----------------------------welcome------------------------------"
    print " NOTE: This is AutoDoge (AutoDock Output Generator / Extractor) 0.11"
    print " NOTE: Written by xywu@PKUHSC in Python 2.7.6, 2016.3.23"
    print " NOTE: -----------------------------welcome------------------------------"
    print "  RUN: Running argparse ..."
    parser = argparse.ArgumentParser(prog='autodock_ops', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-m', type=int, default=3, choices = [1,2,3,4,5],
                        help='''MODE selection (int): 1, 2, 3, 4, 5
      1: Best-binding conformations only
      2: Best-clustering conformations only
      3: Best-binding + best-clustering conformations (default)
      4: Best binding conformation from each cluster
      5: All docking comformations''')
    parser.add_argument('-r', type=int, default=0, choices = [0,1,2],
                        help='''RANKBY selection (int): 0, 1, 2
      0: Rankby name (default)
      1: Rankby binding affinity (from low to high)
      2: Rankby rmsd from reference molecule
         (from low to high)''')
    parser.add_argument('-f', type=float,
                        help='''FILTER setting (float)
      RANKBY needs to be 1 or 2 for FILTER to perform''')
    parser.add_argument('-p', type=str, default="AutoDoge_Output",
                        help='''PREFIX of the output files (str)
      default = \"AutoDoge_Output\"''')
    parser.add_argument('--reponly', action='store_true', default=False,
                        help='''REPORT ONLY (bool): default = False''')
    parser.add_argument('--confonly', action='store_true', default=False,
                        help='''OUTPUT CONF ONLY (bool): default = False''')
    parser.add_argument('--allclus', action='store_true', default=False,
                        help='''ALLCLUS (bool): default = False
      Extract all binding conformations from best-clustering
      and/or best-binding cluster based on the argument -m.
      For minor conformation modifications for 3D-QSAR needs.
      ** -m must be 1, 2 or 3. mode 4 + allclus = mode 5''')
    parser.add_argument('--babel', type=str, default="pdbqt",
                        help='''BABEL (str): default = \"pdbqt\"
      Use OpenBabel to convert the output pdbqt to another
      format, enter valid OpenBabel file format here''')
    parser.add_argument('--doge', action='store_true', default=False,
                        help='''DOGE (bool): default = False
      Print out our generous sponsor at the end of console.''')
    args = parser.parse_args()
    MODE = args.m
    RANKBY = args.r
    FILTER = args.f
    PREFIX = args.p
    REPONLY = args.reponly
    CONFONLY = args.confonly
    BABEL = args.babel
    ALLCLUS = args.allclus
    DOGE = args.doge
    
    if (REPONLY) and (CONFONLY):
        print "ERROR: REPONLY and CONFONLY = True"
        print "ERROR: Enter anything to exit."
        raw_input()
        sys.exit()
    
    ARGCOMMANDOK = 0
    if (MODE == 3) and (RANKBY == 0) and (not FILTER) and (PREFIX == "AutoDoge_Output") and (not REPONLY) and (not CONFONLY) and (BABEL == "pdbqt") and (not ALLCLUS) and (not DOGE):
        print " NOTE: Argparse didn't receive any new parameters."
        print " NOTE: Do you want to use the default settings?"
        print " NOTE:   enter \"exit\" to exit,"
        print " NOTE:   enter anything else to use default settings."
        while (ARGCOMMANDOK == 0):
            ARGCOMMAND = str((raw_input(" USER: COMMAND = ")))
            if ARGCOMMAND == "exit":
                ARGCOMMANDOK = 1
                sys.exit()
            else:
                print "  RUN: --> go default"
                FILTER = 666
                ARGCOMMANDOK = 1
    else:
        if not FILTER:
            FILTER = 666
    return MODE, RANKBY, FILTER, PREFIX, REPONLY, CONFONLY, BABEL, ALLCLUS, DOGE

# Check operating system:
def checkOS():
    IsLinux = 0
    if _platform == "linux" or _platform == "linux2":
        IsLinux = 1; SEP = "//"
        print " USER: OS = Linux"
    else:
        SEP = "\\"
        print " USER: OS = Windows"
    return IsLinux, SEP

# Get dir_name using glob:
def getDlg(IsLinux):
    dlgname = []
    dlgcount = 0
    dlgprint = "none"
    print "  RUN: Detecting dlg_name using wildcard search ..."
    for dlgnm in glob.glob('*.dlg'):
        dlgcount += 1
        if dlgcount == 1:
            dlgprint = dlgnm
        dlgname.append(dlgnm)
    print "  RUN: " + str(dlgcount) + " dlg file(s) found within this folder:"
    if dlgcount != 0:
        print "  RUN: for example, " + dlgprint
    else:
        print "ERROR: No dlg file found."
        if IsLinux == 1:
            sys.exit()
        else:
            print "ERROR: Enter anything to exit ..."
            sys.exit()
    return dlgname

def AutoDoge(DEBUG, dlgname, SEP, MODE, RANKBY, FILTER, PREFIX, REPONLY, CONFONLY, BABEL, ALLCLUS, DOGE):
    # Create REPORT file & new folder(s):
    print "  RUN: -----------------------------running------------------------------"
    # init confdir for reponly debug, dir3 for MODE 4, dir4 for MODE 5:
    confdir1 = "REPONLY"; confdir2 = "REPONLY"; confdir3 = "REPONLY"; confdir4 = "REPONLY"
    def CreateAndWriteFirstLine(reportfilename):
        reportfilewrite = open(reportfilename, "w")
        writefirstline = "MOL_NAME\tCLUSTER_INDEX\tNUM_CONF_IN_CLUSTER\tBIND_AFF\tRUN_RMSREF\tAVG_CLUS_BIND_AFF"
        reportfilewrite.write(writefirstline + '\n')
        return reportfilewrite
    def CreateConfDir(confdir):
        try: os.makedirs(confdir)
        except: print "  RUN: Conformation directory already exists."
        return confdir
    if MODE == 1 or MODE == 3:
        if not CONFONLY: reportfilewrite1 = CreateAndWriteFirstLine(PREFIX + "_BestBindingConf_Report.log")
        if not REPONLY: confdir1 = CreateConfDir(PREFIX + "_BestBindingConf")
    if MODE == 2 or MODE == 3:
        if not CONFONLY: reportfilewrite2 = CreateAndWriteFirstLine(PREFIX + "_BestClusteringConf_Report.log")
        if not REPONLY: confdir2 = CreateConfDir(PREFIX + "_BestClusteringConf")
    if (MODE == 4) or (MODE == 5):
        if not CONFONLY: reportfilewrite3 = CreateAndWriteFirstLine(PREFIX + "_BestAllConf_Report.log")
        if (MODE == 4) and (not REPONLY): confdir3 = CreateConfDir(PREFIX + "_BestAllConf")
    if MODE == 5:
        if not REPONLY: confdir4 = CreateConfDir(PREFIX + "_AllConf")
    
    # Functions for reading file and extracting conformations:
    def ClusteringAnalysis(ExtractingMode):
        if ExtractingMode == 1:
            BBCR = 1
            BBnum = int(ClusterHistogram[9][18:22].replace(" ",""))
            BBLBE = float(ClusterHistogram[9][6:16].replace(" ",""))
            BBMBE = float(ClusterHistogram[9][25:35].replace(" ",""))
            BBNIC = int(ClusterHistogram[9][36:41].replace(" ",""))
            if DEBUG: print "DEBUG: [BB]  clusterrank_index = 1"
            if DEBUG: print "DEBUG: [BB]  model_index = " + str(BBnum)
            return BBCR, BBnum, BBLBE, BBMBE, BBNIC
        if ExtractingMode == 2:
            iterclus = 0
            numberofmolincluster = []
            for eachelement in ClusterHistogram:
                if iterclus >= 9:
                    numberofmolincluster.append(float(eachelement[36:41].replace(" ", "")))
                iterclus += 1
            BCCR = numberofmolincluster.index(max(numberofmolincluster)) + 1
            BCnum = int(ClusterHistogram[BCCR + 8][18:22].replace(" ",""))
            BCLBE = float(ClusterHistogram[BCCR + 8][6:16].replace(" ",""))
            BCMBE = float(ClusterHistogram[BCCR + 8][25:35].replace(" ",""))
            BCNIC = int(ClusterHistogram[BCCR + 8][36:41].replace(" ",""))
            if DEBUG: print "DEBUG: [BC]  clusterrank_index = " + str(BCCR)
            if DEBUG: print "DEBUG: [BC]  model_index = " + str(BCnum)
            return BCCR, BCnum, BCLBE, BCMBE, BCNIC
        if (ExtractingMode == 4) or (ExtractingMode == 5):
            iterclus = 0
            ABCR = []; ABnum = []; ABLBE = []; ABMBE = []; ABNIC = []
            for eachelement in ClusterHistogram:
                if iterclus >= 9:
                    ABCR.append(int(eachelement[0:4].replace(" ", "")))
                    ABnum.append(int(eachelement[18:22].replace(" ", "")))
                    ABLBE.append(float(eachelement[6:16].replace(" ", "")))
                    ABMBE.append(float(eachelement[25:35].replace(" ", "")))
                    ABNIC.append(int(eachelement[36:41].replace(" ", "")))
                iterclus += 1
            if DEBUG: print "DEBUG: [AB]  ClusteringAnalysis in Mode " + str(ExtractingMode) + "."
            if ExtractingMode == 4: return ABCR, ABnum, ABLBE, ABMBE, ABNIC
            else: return ABCR
    def RMSDAnalysis(ModelNumber):
        for eachelement in RMSDTable:
            if int(eachelement[13:18].replace(" ", "")) == ModelNumber:
                rms_ref = float(eachelement[43:50].rstrip().replace(" ",""))
                return rms_ref
                break
    def AllConfLogic():
        ACC = []; ACR = []; ACBE = []; ACRMSD = []
        iterrmsd = 0
        for eachelement in RMSDTable:
            ACC.append(int(eachelement[0:4].replace(" ", "")))
            ACR.append(int(eachelement[13:18].replace(" ", "")))
            ACBE.append(float(eachelement[21:30].replace(" ", "")))
            ACRMSD.append(float(eachelement[42:50].replace(" ", "")))
        return ACC, ACR, ACBE, ACRMSD
    def WriteReport(RankbyMode, ExtractingMode, dlgname_list, clusindex_list, NuminClus_List, LowBE_List, RMSD_List, MedBE_List):
        if RankbyMode == 0:
            for i in range(len(dlgname_list)):
                write_me = str(dlgname_list[i][:-4]) + "\t" + str(clusindex_list[i]) + "\t" + str(NuminClus_List[i]) + "\t" + str(LowBE_List[i]) + "\t" + str(RMSD_List[i]) + "\t" + str(MedBE_List[i])
                if (len(dlgname_list) <= 64) and (ExtractingMode <= 3): print "  RUN: " + write_me
                if ExtractingMode == 1: reportfilewrite1.write(write_me + '\n')
                elif ExtractingMode == 2: reportfilewrite2.write(write_me + '\n')
                elif ExtractingMode == 4: reportfilewrite3.write(write_me + '\n')
                else: print "ERROR: Function WriteReport: ExtractingMode should be 1, 2 or 4."
        else:
            if RankbyMode == 1: rerank = [seqn[0] for seqn in sorted(enumerate(LowBE_List), key=lambda x:x[1])]
            elif RankbyMode == 2: rerank = [seqn[0] for seqn in sorted(enumerate(RMSD_List), key=lambda x:x[1])]
            for eachelement in rerank:
                write_me = str(dlgname_list[eachelement][:-4]) + "\t" + str(clusindex_list[eachelement]) + "\t" + str(NuminClus_List[eachelement]) + "\t" + str(LowBE_List[eachelement]) + "\t" + str(RMSD_List[eachelement]) + "\t" + str(MedBE_List[eachelement])
                if (len(dlgname_list) <= 64) and (ExtractingMode <= 3): print "  RUN: " + write_me
                if ExtractingMode == 1: reportfilewrite1.write(write_me + '\n')
                elif ExtractingMode == 2: reportfilewrite2.write(write_me + '\n')
                elif ExtractingMode == 4: reportfilewrite3.write(write_me + '\n')
                else: print "ERROR: Function WriteReport: ExtractingMode should be 1, 2 or 4."
        return None
    def ExtractModel(ExtractingMode, inputdlgname, ModelNumber, ClusterNumber = 0, AllConfDir = "None"):
        if ExtractingMode == 1:
            conffilewrite = open(confdir1 + SEP + inputdlgname[:-4] + "_BestBindingConf.pdbqt", "w")
            for eachelement in ConfMatrix[ModelNumber]: conffilewrite.write(eachelement[8:])
        if ExtractingMode == 2:
            conffilewrite = open(confdir2 + SEP + inputdlgname[:-4] + "_BestClusteringConf.pdbqt", "w")
            for eachelement in ConfMatrix[ModelNumber]: conffilewrite.write(eachelement[8:])
        if ExtractingMode == 4:
            conffilewrite = open(confdir3 + SEP + inputdlgname[:-4] + SEP + inputdlgname[:-4] + "_CLUSTER" + str(ClusterNumber) + "_BestBindingConf.pdbqt", "w")
            for eachelement in ConfMatrix[ModelNumber]: conffilewrite.write(eachelement[8:])
        if ExtractingMode == 5:
            conffilewrite = open(AllConfDir + SEP + inputdlgname[:-4] + "_CLUSTER" + str(ClusterNumber) + "_RUN" + str(ModelNumber) + ".pdbqt", "w")
            for eachelement in ConfMatrix[ModelNumber]: conffilewrite.write(eachelement[8:])
        if ExtractingMode == 12:
            conffilewrite = open(AllConfDir + SEP + inputdlgname[:-4] + "_RUN" + str(ModelNumber) + ".pdbqt", "w")
            for eachelement in ConfMatrix[ModelNumber]: conffilewrite.write(eachelement[8:])
        return None
    debugprint = [["", "[BB]", "[BC]"], ["", "", "", "", "[AB]", "[AC]", "", "", "", "", "", "", "[BBCA]"]]
    def RankbyAndFilter(ExtractingMode, IsAllConf, RankbyMode, FilterValue, BAValue, RMSDValue, ConfIndex, ClusterNumber = 0, AllConfDir = "None"):
        if (ExtractingMode <= 2) and (IsAllConf == 1): ExtractingMode = 12
        if (not REPONLY) or (IsAllConf == 1): 
            if FilterValue == 666:
                if ExtractingMode == 4: ExtractModel(ExtractingMode, inputdlg, ConfIndex, ClusterNumber)
                elif (ExtractingMode == 5) or (ExtractingMode == 12): ExtractModel(ExtractingMode, inputdlg, ConfIndex, ClusterNumber, AllConfDir)
                else: ExtractModel(ExtractingMode, inputdlg, ConfIndex)
            else:
                if RankbyMode == 1:
                    if BAValue <= FilterValue:
                        if ExtractingMode == 4: ExtractModel(ExtractingMode, inputdlg, ConfIndex, ClusterNumber)
                        elif (ExtractingMode == 5) or (ExtractingMode == 12): ExtractModel(ExtractingMode, inputdlg, ConfIndex, ClusterNumber, AllConfDir)
                        else: ExtractModel(ExtractingMode, inputdlg, ConfIndex)
                        if DEBUG: print "DEBUG: " + debugprint[IsAllConf][ExtractingMode] + " RANKBY = BA -> Keep conformation."
                    else:
                        if DEBUG: print "DEBUG: " + debugprint[IsAllConf][ExtractingMode] + " RANKBY = BA -> Abandon conformation."
                elif RankbyMode == 2:
                    if RMSDValue <= FilterValue:
                        if ExtractingMode == 4: ExtractModel(ExtractingMode, inputdlg, ConfIndex, ClusterNumber)
                        elif (ExtractingMode == 5) or (ExtractingMode == 12): ExtractModel(ExtractingMode, inputdlg, ConfIndex, ClusterNumber, AllConfDir)
                        else: ExtractModel(ExtractingMode, inputdlg, ConfIndex)                    
                        if DEBUG: print "DEBUG: " + debugprint[IsAllConf][ExtractingMode] + " RANKBY = RMS -> Keep molecule conformation."
                    else:
                        if DEBUG: print "DEBUG: " + debugprint[IsAllConf][ExtractingMode] + " RANKBY = RMS -> Abandon molecule conformation."                    
                else: print "ERROR: " + debugprint[IsAllConf][ExtractingMode] + " RANKBY = " + str(RankbyMode) + ", should be 1 or 2 for FILTER to work!"
        if REPONLY:
            if (DEBUG) and (IsAllConf == 0): print "DEBUG: " + debugprint[IsAllConf][ExtractingMode] + " REPONLY enabled."    
        return None
    
    # Main conditions:
    i = 0
    # for rankby and filter:
    BBclusindex_List = []; BBconfindex_List = []; BBLowBE_List = []; BBMedBE_List = []; BBNuminClus_List = []; BBRMSD_List = []
    BCclusindex_List = []; BCconfindex_List = []; BCLowBE_List = []; BCMedBE_List = []; BCNuminClus_List = []; BCRMSD_List = []
    CorruptedFiles = []
    # Main condition:
    for i in range(len(dlgname)):
        inputdlg = dlgname[i]
        if (i % 100 == 1) and (i != 1): print "  RUN: " + str(i) + " dlg files scanned."
        if DEBUG: print "DEBUG: Processing name: " + dlgname[i][:-4]
        confid = 1
        confwrite = 0
        cluswrite = 0
        rmsdwrite = 0
        CorruptedFlag = 1
        # Get Number of GA Run to build the matrix:
        for eachline in fileinput.input(inputdlg):
            if "DPF> ga_run " in eachline:
                num_ga = int(eachline[12:16].replace(" ",""))
                if DEBUG: print "DEBUG: Number of GA Run = " + str(num_ga)
                break
        fileinput.close()
        # Build the matrix according to num_ga:
        ConfMatrix = [[] for i in range(num_ga + 1)] # Wrong: [[]] * (num_ga + 1) -> error on append elements afterwards.
        ClusterHistogram = []
        RMSDTable = []
        for eachline in fileinput.input(inputdlg):
            # Fill the the matrix using real confid:
            if ("DOCKED: MODEL" in eachline) and (confwrite == 0): confwrite = 1
            if confwrite == 1: ConfMatrix[confid].append(eachline) # Not using rstrip to make writing easier
            if ("DOCKED: ENDMDL" in eachline) and (confwrite == 1): confwrite = 0; confid += 1
            # Get the Clustering Histogram:
            if "CLUSTERING HISTOGRAM" in eachline: cluswrite = 1
            if "_____|___________|_____|___________|_____|__________" in eachline: cluswrite = 0
            if cluswrite == 1: ClusterHistogram.append(eachline.rstrip())
            # Get the RMSD Table:
            if "RMSD TABLE" in eachline: rmsdwrite = 1
            if "INFORMATION ENTROPY ANALYSIS FOR THIS CLUSTERING" in eachline: rmsdwrite = 0
            if rmsdwrite == 1: RMSDTable.append(eachline.rstrip())
            # Check file integrity:
            if "autodock4: Successful Completion" in eachline:
                CorruptedFlag = 0
                if DEBUG: print "DEBUG: File integrity = true"
        fileinput.close()
        RMSDTable = RMSDTable[8:-4]
        # Calculate the model number, append to list:
        BBconfindex = 0; BCconfindex = 0
        if MODE == 1 or MODE == 3:
            BBclusindex, BBconfindex, BBLowBE, BBMedBE, BBNuminClus = ClusteringAnalysis(1)
            BBRMSD = RMSDAnalysis(BBconfindex)
            BBclusindex_List.append(BBclusindex); BBconfindex_List.append(BBconfindex); BBLowBE_List.append(BBLowBE); BBMedBE_List.append(BBMedBE); BBNuminClus_List.append(BBNuminClus)
            BBRMSD_List.append(BBRMSD)
            RankbyAndFilter(1, 0, RANKBY, FILTER, BBLowBE, BBRMSD, BBconfindex)
        if MODE == 2 or MODE == 3:
            BCclusindex, BCconfindex, BCLowBE, BCMedBE, BCNuminClus = ClusteringAnalysis(2)
            BCRMSD = RMSDAnalysis(BCconfindex)
            BCclusindex_List.append(BCclusindex); BCconfindex_List.append(BCconfindex); BCLowBE_List.append(BCLowBE); BCMedBE_List.append(BCMedBE); BCNuminClus_List.append(BCNuminClus)
            BCRMSD_List.append(BCRMSD)
            RankbyAndFilter(2, 0, RANKBY, FILTER, BCLowBE, BCRMSD, BCconfindex)      
        if (MODE == 4) or (MODE == 5): # Codes for MODE 4 (Report + Extract Model)
            ABclusindex = []; ABconfindex = []; ABLowBE = []; ABMedBE = []; ABNuminClus = []; ABRMSD = []
            ABclusindex, ABconfindex, ABLowBE, ABMedBE, ABNuminClus = ClusteringAnalysis(4)
            ABName = [inputdlg] * len(ABclusindex)
            for eachelement in ABconfindex: ABRMSD.append(RMSDAnalysis(eachelement))
            if not REPONLY:
                if MODE == 4:
                    secondarydirs = confdir3 + SEP + inputdlg[:-4]
                    try: os.makedirs(secondarydirs)
                    except: pass
                    for i in range(len(ABclusindex)): RankbyAndFilter(4, 1, RANKBY, FILTER, ABLowBE[i], ABRMSD[i], ABconfindex[i], ABclusindex[i])
            # Write report for both MODE 4 and MODE 5:
            if not CONFONLY: WriteReport(RANKBY, 4, ABName, ABclusindex, ABNuminClus, ABLowBE, ABRMSD, ABMedBE)
        if (MODE == 5) and (not REPONLY): # Codes for MODE 5 (Extract Model Only, Report written in MODE 4)
            secondarydirs = confdir4 + SEP + inputdlg[:-4]
            try: os.makedirs(secondarydirs)
            except: pass
            ACclusindex = []; ACclusindex = ClusteringAnalysis(5)
            for eachelement in ACclusindex:
                threetierdirs = secondarydirs + SEP + "CLUSTER" + str(eachelement)
                try: os.makedirs(threetierdirs)
                except: pass
            ACclusscaledindex = []; ACrunindex = []; ACBindAff = []; ACRMSD = []
            ACclusscaledindex, ACrunindex, ACBindAff, ACRMSD = AllConfLogic()
            for i in range(len(ACclusscaledindex)):
                WriteInHere = secondarydirs + SEP + "CLUSTER" + str(ACclusscaledindex[i])
                RankbyAndFilter(5, 1, RANKBY, FILTER, ACBindAff[i], ACRMSD[i], ACrunindex[i], ACclusscaledindex[i], WriteInHere)
        AllConfDir1 = "REPONLY"; AllConfDir2 = "REPONLY"
        if ALLCLUS: # Codes for ALLCLUS in MODE 1 to 3:
            if (MODE == 4) or (MODE == 5):
                print "ERROR: ALLCLUS not available in MODE 4 or 5."
                print "ERROR: ALLCLUS: MODE 4 + ALLCLUS = MODE 5."
            else:
                ACclusscaledindex = []; ACrunindex = []; ACBindAff = []; ACRMSD = []
                ACclusscaledindex, ACrunindex, ACBindAff, ACRMSD = AllConfLogic()
                if (MODE == 1) or (MODE == 3):
                    AllConfDir1 = PREFIX + "_BestBindingCluster"
                    secondarydirs = AllConfDir1 + SEP + inputdlg[:-4]
                    try: os.makedirs(secondarydirs)
                    except: pass
                    for j in range(len(ACclusscaledindex)):
                        if ACclusscaledindex[j] == 1: RankbyAndFilter(1, 1, RANKBY, FILTER, ACBindAff[j], ACRMSD[j], ACrunindex[j], ACclusscaledindex[j], secondarydirs)
                if (MODE == 2) or (MODE == 3):
                    AllConfDir2 = PREFIX + "_BestClustingCluster"
                    secondarydirs = AllConfDir2 + SEP + inputdlg[:-4]
                    try: os.makedirs(secondarydirs)
                    except: pass
                    for j in range(len(ACclusscaledindex)):
                        if ACclusscaledindex[j] == BCclusindex: RankbyAndFilter(2, 1, RANKBY, FILTER, ACBindAff[j], ACRMSD[j], ACrunindex[j], ACclusscaledindex[j], secondarydirs)
    
    # Write report (MODE <= 3):
    def ReportLegends():
        print "  RUN: ------------------------------report------------------------------"
        print "  RUN: MN = Molecule Name; CI = Culster Index;"
        print "  RUN: NCIC = Number of Conformations in this Cluster;"
        print "  RUN: BA = (Lowest) Binding Affinity (of this Cluster);"
        print "  RUN: RR = RMSD to Reference molecule;"
        print "  RUN: ABAC = Average Binding Affinity of this Cluster;"  
        return None
    if (not CONFONLY) and (MODE <= 3):
        if (MODE == 1) or (MODE == 3):
            ReportLegends() 
            print "  RUN: Resorting best-binding conf report on RANKBY = " + str(RANKBY)
            if len(dlgname) <= 64:
                ReportLegends()
                print "  RUN: MN\tCI\tNCIC\tBA\tRR\tABAC"
            else: print "  RUN: Number of molecules processed > 64"
            WriteReport(RANKBY, 1, dlgname, BBclusindex_List, BBNuminClus_List, BBLowBE_List, BBRMSD_List, BBMedBE_List)
            reportfilewrite1.close() # close report file for ALLCONF to run properly:
        if (MODE == 2) or (MODE == 3):
            print "  RUN: Resorting best-clustering conf report on RANKBY = " + str(RANKBY)
            if len(dlgname) <= 64:
                ReportLegends()
                print "  RUN: MN\tCI\tNCIC\tBA\tRR\tABAC"
            else: print "  RUN: Number of molecules processed > 64"
            WriteReport(RANKBY, 2, dlgname, BCclusindex_List, BCNuminClus_List, BCLowBE_List, BCRMSD_List, BCMedBE_List)
            reportfilewrite2.close() # close report file for ALLCONF to run properly:
        if MODE == 4: reportfilewrite3.close() # close report file for ALLCONF to run properly:
    if CONFONLY:
        if DEBUG: print "DEBUG: CONFONLY enabled."
    
    # Babel:
    if (BABEL != "pdbqt") and (not REPONLY):
        # Babel init parameters:
        inputff = "-ipdbqt"
        outputff = "-o" + BABEL
        foronce1 = 0; foronce2 = 0; foronce3 = 0; foronce4 = 0
        # Functions for Babel:
        def printfolderproblem():
            print " NOTE: BABEL: output folder already exists."
        def printbabelproblem():
            print "ERROR: BABEL: while converting file format."
            print "ERROR: BABEL: Babel not found."
        def BabelMode123(confdir):
            confdirall = confdir + SEP + "*.pdbqt"
            confdirbabel = confdir + "_" + BABEL
            confdirbabelall = confdirbabel + SEP + "*." + BABEL
            try: os.makedirs(confdirbabel)
            except: printfolderproblem()
            if DEBUG: print "DEBUG: babel " + inputff + " " + confdirall + " " + outputff + " " + confdirbabelall
            try: subprocess.call(["babel", inputff, confdirall, outputff, confdirbabelall])
            except: printbabelproblem()
            return None
        def BabelMode4AllClus(AllConfDir):
            global foronce1; global foronce2 # Global the variant in function to avoid "local variable 'xxx' referenced before assignment"
            babeloutputroot = AllConfDir + "_" + BABEL
            try: os.makedirs(babeloutputroot)
            except: printfolderproblem()
            for nm in dlgname:
                babelname = nm[:-4]
                babelinput = AllConfDir + SEP + babelname + SEP + "*.pdbqt"
                babeloutputbranch = babeloutputroot + SEP + babelname
                try: os.makedirs(babeloutputbranch)
                except:
                    if foronce1 == 0: 
                        printfolderproblem()
                        foronce1 = 1
                babeloutput = babeloutputbranch + SEP + "*." + BABEL
                try: subprocess.call(["babel", inputff, babelinput, outputff, babeloutput])
                except: 
                    if foronce2 == 0: 
                        printbabelproblem()
                        foronce2 = 1
            return None
        def BabelMode5(AllDir):
            global foronce3; global foronce4 # Global the variant in function to avoid "local variable 'xxx' referenced before assignment"
            babeloutputroot = AllDir + "_" + BABEL
            try: os.makedirs(babeloutputroot) # level 1 folder
            except: printfolderproblem()
            for nm in dlgname:
                babelname = nm[:-4]
                babelinputbranch = AllDir + SEP + babelname
                babeloutputbranch = babeloutputroot + SEP + babelname
                try: os.makedirs(babeloutputbranch) # level 2 folder
                except:
                    if foronce3 == 0:
                        printfolderproblem()
                        foronce3 = 1
                dirroot = os.getcwd() + SEP + babelinputbranch
                tridirs = [x[1] for x in os.walk(dirroot)]
                quadfiles = [x[2] for x in os.walk(dirroot)]
                k = 0 # k loop: to ameliorate a babel bug.
                for eachdir in tridirs[0]:
                    k += 1
                    babeloutputtree = babeloutputbranch + SEP + eachdir
                    try: os.makedirs(babeloutputtree) # level 3 folder
                    except: pass # one error line is enough.
                    if len(quadfiles[k]) > 1:
                        babelinput = babelinputbranch + SEP + eachdir + SEP + "*.pdbqt"
                        babeloutput = babeloutputbranch + SEP + eachdir + SEP + "*." + BABEL
                    if len(quadfiles[k]) == 1:
                        babelinput = babelinputbranch + SEP + eachdir + SEP + quadfiles[k][0]
                        babeloutput = babeloutputbranch + SEP + eachdir + SEP + quadfiles[k][0][:-6] + "." + BABEL
                    try: subprocess.call(["babel", inputff, babelinput, outputff, babeloutput])
                    except: 
                        if foronce4 == 0:
                            printbabelproblem()
                            foronce4 = 1
            return None
        # Main conditions:
        print "  RUN: BABEL: from pdbqt to " + BABEL
        if confdir1 != "REPONLY": BabelMode123(confdir1)
        if confdir2 != "REPONLY": BabelMode123(confdir2)
        if AllConfDir1 != "REPONLY": BabelMode4AllClus(AllConfDir1)
        if AllConfDir2 != "REPONLY": BabelMode4AllClus(AllConfDir2)
        if confdir3 != "REPONLY": BabelMode4AllClus(confdir3)
        if confdir4 != "REPONLY": BabelMode5(confdir4)
    if (BABEL != "pdbqt") and (REPONLY): print "ERROR: BABEL: REPONLY should be False."
    
    # DOGE:
    if DOGE:
        print " DOGE:                      Y.                      _   "
        print " DOGE:                      YiL                   .```.  "
        print " DOGE:                      Yii;                .; .;;`.    "
        print " DOGE:                      YY;ii._           .;`.;;;; :    "
        print " DOGE:                      iiYYYYYYiiiii;;;;i` ;;::;;;;    "
        print " DOGE:                  _.;YYYYYYiiiiiiYYYii  .;;.   ;;; "
        print " DOGE:               .YYYYYYYYYYiiYYYYYYYYYYYYii;`  ;;;;    "
        print " DOGE:             .YYYYYYY$$YYiiYY$$$$iiiYYYYYY;.ii;`..   "
        print " DOGE:            :YYY$!.  TYiiYY$$$$$YYYYYYYiiYYYYiYYii.    "
        print " DOGE:            Y$MM$:   :YYYYYY$!\"``\"4YYYYYiiiYYYYiiYY.    "
        print " DOGE:         `. :MM$$b.,dYY$$Yii\" :'   :YYYYllYiiYYYiYY    "
        print " DOGE:      _.._ :`4MM$!YYYYYYYYYii,.__.diii$$YYYYYYYYYYY"
        print " DOGE:      .,._ $b`P`     \"4$$$$$iiiiiiii$$$$YY$$$$$$YiY;"
        print " DOGE:         `,.`$:       :$$$$$$$$$YYYYY$$$$$$$$$YYiiYYL"
        print " DOGE:          \"`;$$.    .;PPb$`.,.``T$$YY$$$$THANKSiiiYYU:  "
        print " DOGE:          ;$P$;;: ;;;;i$y$\"!Y$$$b;$$$Y$YY$$FORiiiYYiYY "
        print " DOGE:          $Fi$$ .. ``:iii.`-\":YYYYY$$YYUSINGYYYiiYiYYY    "
        print " DOGE:          :Y$$rb ````  `_..;;i;YYY$YY$$$$$$$AUTODOGEYY:    "
        print " DOGE:           :$$$$$i;;iiiiidYYYYYYYYYY$$$$$$YYYYYYYiiYYYY. "
        print " DOGE:            `$$$$$$$YYYYYYYYYYYYY$$$$$$YYYYYYYYiiiYYYYYY    "
        print " DOGE:            .i!$$$$$$YYYYYYYYY$$$$$$YYY$$YYiiiiiiYYYYYYY    "
        print " DOGE:           :YYiii$$$$$$$YYYYYYY$$$$YY$$$$YYiiiiiYYYYYYi'    "
        print " DOGE:           Such completion.       Enter anything to wow."
    if not DOGE: print "  RUN: Task complete. Enter anything to exit."

def main():
    MODE, RANKBY, FILTER, PREFIX, REPONLY, CONFONLY, BABEL, ALLCLUS, DOGE = argParser()
    IsLinux, SEP = checkOS()
    dlgname = getDlg(IsLinux)
    start = timeit.default_timer()
    AutoDoge(DEBUG, dlgname, SEP, MODE, RANKBY, FILTER, PREFIX, REPONLY, CONFONLY, BABEL, ALLCLUS, DOGE)
    stop = timeit.default_timer()
    print " NOTE: -----------------------------complete-----------------------------"
    print " NOTE: TASK RUNNING TIME = " + str(round(stop - start, 2)) + " sec."

if __name__ == '__main__':
    main()