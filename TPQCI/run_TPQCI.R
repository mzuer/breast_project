##################################################################
#Script to calculate TPQCI and detect CNV influenced modules     #
#Code By Yusong Liu                                              # 
#Mar 24, 2021                                                    #
#                                                                #
#If any part of this code is used in publishable works,          #
#please citing:                                                  #
#                                                                #
#   Y. Liu et al., "TPQCI: A Topology Potential-based Method to  #
#   Quantify Functional Influence of Copy Number Variations",    #
#   Methods (2021)                                               #
#                                                                #
##################################################################   

source("getlv1net.r")
source("calcModule.R")

#######
#Input#
#######

#Cautions: All this part are virtual declarations, need to be replaced to actual data as need

#PPI is a iGraph class object
ppi<-make_graph()
#CNV is a data frame that each row is a gene, each column is a sample
#the first column is gene name with a name "Gene"
cnv<-data.frame()
#RNA is also a data frame that each row is a gene, each column is a sample
#the first column is gene name with a name "Gene"
rna<-data.frame()

############
#Parameters#
############

#if RNA-seq is logarithmic, set this flag to T
is.RNA.log2ed<-F

#parameters while calculating TPQCI
cnv.mass.threshold<-0.3
pcc.threshold<-0.3

#parameters while module detecting

#Obsoleted parameters
#high.degree.ratio.threshold<-0.5
#lm.min.neighbors<-0

tau.g<-0.15
tau.l<-0.3
beta<-0.5
mu<-20


path.output<-"modules.csv"

#######################
#Run:Calculating TPQCI#
#######################

if(!is.RNA.log2ed)
  rna<-rnaLog2Trans(rna)

cr.cor<-getLV1cor(cnv=cnv, rna=rna, ppn=ppi)
mass<-getCNVMass(cnv,cnv.mass.threshold)

tpqci<-calTopoPot(mass,cr,cor,pcc.threshold)


######################
#Run:Module detection#
######################
tp.weight<-tpqci$topot
names(tp.weight)<-tpqci$node_c

tw.ppi<-getWeightedNetwork(ppi,tp.weight)

#Remove high degree nodes if uncommented (Obsoleted)
#tw.ppi<-removeHighDegNodes(tw.ppi,high.degree.ratio.threshold)

lm.genes<-findLocalMax(tw.ppi)
#switch to this if need limit n.neighbors of local maximal nodes(Obsoleted)
#lm.genes<-findLocalMax(tw.ppi, lm.min.neighbors)

modules<-findModules(tw.ppi,lm.genes,tau.g,tau.l,mu)
modules<-moduleMerge(modules,beta)

#if need to output the modules, uncomment next line
#writeModToCSV(modules, path.output)

