#!/usr/local/bin/python

# Name: AutoDoge (AUTODock Output Generator / Extractor)
# Author: xywu@PKUHSC (wuxingyu@bjmu.edu.cn)
# Version: 0.12 (7/12/2016)
# Written in Python 2.7.12 (Anaconda)
# Description: Generate AutoDock4 VS report and extract best-binding conformations
#              for best-binding or best-clustering cluster.
# Operating systems supported: Windows / Linux

import Tkinter
from Tkinter import *
import tkFileDialog
import glob
import os
import sys
from sys import platform as _platform
import subprocess
import fileinput

####################### SET DEBUG #######################
DEBUG = 0
####################### SET DEBUG #######################

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
        if CorruptedFlag:
            CorruptedFiles.append(inputdlg)
            if DEBUG: print "DEBUG: File corrupted."
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
        if ALLCLUS:
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
    if not DOGE: print "  RUN: Task complete."
    
    return CorruptedFiles

###############################################################
############################# GUI #############################
###############################################################
if __name__ == '__main__':
    form = Tkinter.Tk()

    getFld = Tkinter.IntVar()

    form.wm_title('AutoDoge 0.12 GUI')
    
    # labelFrame
    stepOne = Tkinter.LabelFrame(form, text=\
                " 1. Working Directory, Mode Selection, and Filter Criteria: ")
    stepOne.grid(row=0, columnspan=10, sticky='WE', \
                 padx=5, pady=5, ipadx=5, ipady=5)

    stepTwo = Tkinter.LabelFrame(form, text=" 2. Output Configuration: ")
    stepTwo.grid(row=2, columnspan=10, sticky='WE', \
                 padx=5, pady=5, ipadx=5, ipady=5)

    stepThree = Tkinter.LabelFrame(form, text=" 3. Run AutoDoge: ")
    stepThree.grid(row=3, columnspan=10, sticky='WE', \
                   padx=5, pady=5, ipadx=5, ipady=5)
    
    # stepOne: Set Working Directory and Running Mode:
    dname = StringVar(stepOne, value=os.getcwd())
    def browseDir():
        dname = tkFileDialog.askdirectory(initialdir = os.getcwd())
        inDirEnt.delete(0, END)
        inDirEnt.insert(0, dname)
        return
    inDirLbl = Tkinter.Label(stepOne, text="Set Working Directory:")
    inDirLbl.grid(row=0, column=0, sticky='E', padx=5, pady=2)
    inDirEnt = Tkinter.Entry(stepOne, textvariable=dname)
    inDirEnt.grid(row=0, column=1, columnspan=7, sticky="WE", pady=3)
    inDirBtn = Tkinter.Button(stepOne, text="Browse ...", command=browseDir)
    inDirBtn.grid(row=0, column=8, sticky='W', padx=5, pady=2)

    def getMode(event):
        if modeSelection.get() == "4 - Best-Binding Conformation From Each Cluster"\
        or modeSelection.get() == "5 - All Conformations, Sorted by Cluster":
            allconfChk.deselect()
            allconfChk.configure(state='disabled')
        else:
            allconfChk.configure(state='normal')
        return
    modeSelectionLbl = Tkinter.Label(stepOne, text="Set Running Mode:")
    modeSelectionLbl.grid(row=1, column=0, sticky='E', padx=5, pady=2)
    modeSelection = StringVar()
    modeSelection.set("3 - Best-Binding + Best-Clustering Conformations")
    modeSelectionOM = OptionMenu(stepOne, modeSelection,
                        "1 - Best-Binding Conformations",
                        "2 - Best-Clustering conformations",
                        "3 - Best-Binding + Best-Clustering Conformations",
                        "4 - Best-Binding Conformation From Each Cluster",
                        "5 - All Conformations, Sorted by Cluster",
                        command=getMode)
    modeSelectionOM.grid(row=1, column=1, columnspan=6, sticky="WE", pady=2)

    def getRankby(event):
        if rankbySelection.get() == "0 - By File Name":
            filterLblVar.set("Filter Criteria: Not Applicapable")
            filterEnt.grid_forget()
        else:
            filterEnt.grid(row=2, column=6, columnspan=4, sticky="WE", pady=3)
            filterLblVar.set("Filter Criteria: " + rankbySelection.get()[7:] + " Less Than")
        return
    rankbySelectionLbl = Tkinter.Label(stepOne, text="Set Sorting Mode:")
    rankbySelectionLbl.grid(row=2, column=0, sticky='E', padx=5, pady=2)
    rankbySelection = StringVar()
    rankbySelection.set("0 - By File Name")
    rankbySelectionOM = OptionMenu(stepOne, rankbySelection,
                          "0 - By File Name",
                          "1 - By Binding Affinity",
                          "2 - By RMSref",
                          command=getRankby)
    rankbySelectionOM.grid(row=2, column=1, columnspan=4, sticky="W", pady=2)
    
    filterLblVar = StringVar()
    filterLblVar.set("Filter Criteria: Not Applicapable")
    filterLbl = Tkinter.Label(stepOne, textvariable=filterLblVar)
    filterLbl.grid(row=2, column=5, sticky='E', padx=5, pady=2)
    filterEnt = Tkinter.Entry(stepOne)
    filterEnt.grid(row=2, column=6, columnspan=3, sticky="WE", pady=3)
    filterEnt.grid_forget()
    
    allconfSelection = IntVar()
    allconfSelection.set(0)
    dogeSelection = IntVar()
    dogeSelection.set(0)
    allconfChk = Tkinter.Checkbutton(stepOne, \
                           text="Output All Conformations from the Cluster",\
                           onvalue=1, offvalue=0, variable=allconfSelection)
    allconfChk.grid(row=3, column=0, columnspan=3, pady=2, sticky='W')
    dogeChk = Tkinter.Checkbutton(stepOne, \
                           text="Need Doge?",\
                           onvalue=1, offvalue=0, variable=dogeSelection)
    dogeChk.grid(row=3, column=4, columnspan=3, pady=2, sticky='W')
    
    def ronly():
        if confonlySelection: confonlyChk.deselect()
        if reponlySelection.get() == 1:
            boffLblVar.set("Babel Output File Format: Not Applicapable")
            boffSelection.set("pdbqt")
            boffEnt.grid_forget()
        if reponlySelection.get() == 0:
            boffLblVar.set("Babel Output File Format:")
            boffEnt.grid(row=5, column=6, columnspan=4, sticky="WE", pady=3)
        return
    def conly():
        if reponlySelection:
            reponlyChk.deselect()
            boffLblVar.set("Babel Output File Format:")
            boffEnt.grid(row=5, column=6, columnspan=4, sticky="WE", pady=3)    
        return
    reponlySelection = IntVar()
    reponlySelection.set(0)
    confonlySelection = IntVar()
    confonlySelection.set(0)
    reponlyChk = Tkinter.Checkbutton(stepOne, \
                           text="Output Report Only",\
                           onvalue=1, offvalue=0, variable=reponlySelection,\
                           command=ronly)
    reponlyChk.grid(row=4, column=0, columnspan=3, pady=2, sticky='W')
    confonlyChk = Tkinter.Checkbutton(stepOne, \
                           text="Output Conformations Only",\
                           onvalue=1, offvalue=0, variable=confonlySelection,\
                           command=conly)
    confonlyChk.grid(row=4, column=4, columnspan=3, pady=2, sticky='W')
    
    # stepTwo: Output Configuration:
    prefixLbl = Tkinter.Label(stepTwo, text="Output Prefix:")
    prefixLbl.grid(row=5, column=0, sticky='E', padx=5, pady=2)
    prefixSelection = StringVar()
    prefixSelection.set("AutoDoge_Output")
    prefixEnt = Tkinter.Entry(stepTwo, textvariable=prefixSelection)
    prefixEnt.grid(row=5, column=1, columnspan=4, sticky="WE", pady=3)
    
    boffLblVar = StringVar()
    boffLblVar.set("Babel Output File Format:")
    boffLbl = Tkinter.Label(stepTwo, textvariable=boffLblVar)
    boffLbl.grid(row=5, column=5, sticky='WE', padx=5, pady=2)
    boffSelection = StringVar()
    boffSelection.set("pdbqt")
    boffEnt = Tkinter.Entry(stepTwo, textvariable=boffSelection)
    boffEnt.grid(row=5, column=6, columnspan=4, sticky="WE", pady=3)

    #stepThree: Run AutoDoge
    def runAutoDoge():
        if DEBUG:
            print "  RUN: ------------------------------input-------------------------------"
            print "  RUN: WORKINGDIR = " + inDirEnt.get()
            print "  RUN: MODE = " + modeSelection.get()[0]
            print "  RUN: RANKBY = " + rankbySelection.get()[0]
            print "  RUN: FILTER = " + str(filterEnt.get())
            print "  RUN: PREFIX = " + prefixSelection.get()
            print "  RUN: BABEL = " + boffSelection.get()
            print "  RUN: ALLCONF = " + str(allconfSelection.get())
            print "  RUN: REPONLY = " + str(reponlySelection.get())
            print "  RUN: CONFONLY = " + str(confonlySelection.get())
            print "  RUN: DOGE = " + str(dogeSelection.get())
            print "  RUN: ------------------------------input-------------------------------"
        IsLinux, SEP = checkOS()
        os.chdir(inDirEnt.get())
        dlgname = getDlg(IsLinux)
        if dlgname == []:
            autoDogeLblVar.set("ERROR: No dlg found in folder.")
            return
        else:
            FILTER = 666
            if filterEnt.get() != "":
                FILTER = filterEnt.get()
                FILTER = float(FILTER)
            corruptedFiles = AutoDoge(DEBUG, dlgname, SEP,\
                             int(modeSelection.get()[0]),\
                             int(rankbySelection.get()[0]),\
                             FILTER, prefixSelection.get(),\
                             reponlySelection.get(), confonlySelection.get(),\
                             boffSelection.get(), allconfSelection.get(),\
                             dogeSelection.get())
            if len(corruptedFiles) != 0:
                autoDogeLblVar.set("WARNING: Task completed. Have corrupted file(s).")
                if DEBUG:
                    print "DEBUG: Following file(s) is/are corrupted:"
                    print corruptedFiles
            else:
                if dogeSelection.get():
                    autoDogeLblVar.set("Wow such complete. Sanic speed.")
                else:
                    autoDogeLblVar.set("Task completed.")
    autoDogeBtn = Tkinter.Button(stepThree, text="RUN AUTODOGE", command=runAutoDoge)
    autoDogeBtn.grid(row=6, column=0, sticky='E', padx=5, pady=2)
    autoDogeLblVar = StringVar()
    autoDogeLblVar.set("Ready to run.")
    autoDogeLbl = Tkinter.Label(stepThree, textvariable=autoDogeLblVar)
    autoDogeLbl.grid(row=6, column=1, sticky='E', padx=5, pady=2)

    form.mainloop()