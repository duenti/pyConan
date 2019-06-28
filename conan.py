import xml.etree.ElementTree as ET
from Levenshtein import ratio
import numpy as np
import os
import math
from decimal import Decimal
from scipy.stats import hypergeom
from scipy.misc import logsumexp
import scipy.spatial.distance
import random
import time
import sys
import pandas as pd
import networkx as nx
from itertools import combinations
import itertools
from operator import itemgetter
import getopt

MIN_JC = 0.5
FREQ_WINDOW = 0.10

properties = {}
properties[0] = "A"
properties[1] = "C"
properties[2] = "D"
properties[3] = "E"
properties[4] = "F"
properties[5] = "G"
properties[6] = "H"
properties[7] = "I"
properties[8] = "K"
properties[9] = "L"
properties[10] = "M"
properties[11] = "N"
properties[12] = "P"
properties[13] = "Q"
properties[14] = "R"
properties[15] = "S"
properties[16] = "T"
properties[17] = "V"
properties[18] = "W"
properties[19] = "Y"
properties[20] = "Amide"
properties[21] = "Aliphatic"
properties[22] = "Basic"
properties[23] = "Hydroxyl"
properties[24] = "Sulfur"
properties[25] = "Non-Polar"
properties[26] = "Polar"
properties[27] = "Hydrophobic"
properties[28] = "Hydrophilic"
properties[29] = "Pos.Charged"
properties[30] = "Neg.Charged"
properties[31] = "VerySmall"
properties[32] = "Small"
properties[33] = "MediumA"
properties[34] = "MediumB"
properties[35] = "Aromatic"
properties[36] = "ND"
properties[37] = "QE"

setsIdList = {}
setsIdList["A"] = [0,21,25,31]
setsIdList["C"] = [1,24,26,32]
setsIdList["D"] = [2,28,30,32,36]
setsIdList["E"] = [3,30,33,37]
setsIdList["F"] = [4,25,27,35]
setsIdList["G"] = [5,21,25,31]
setsIdList["H"] = [6,22,33]
setsIdList["I"] = [7,25,27,34]
setsIdList["K"] = [8,22,28,29,34]
setsIdList["L"] = [9,21,25,27,34]
setsIdList["M"] = [10,24,25,27,34]
setsIdList["N"] = [11,20,26,28,32,36]
setsIdList["P"] = [12,25,28,32]
setsIdList["Q"] = [13,20,26,28,33,37]
setsIdList["R"] = [14,22,28,29,34]
setsIdList["S"] = [15,23,26,31]
setsIdList["T"] = [16,23,26,32]
setsIdList["V"] = [17,21,25,27,33]
setsIdList["W"] = [18,25,27,35]
setsIdList["Y"] = [19,21,23,26,35]

setsAA = {}
setsAA[0] = "A"
setsAA[1] = "C"
setsAA[2] = "D"
setsAA[3] = "E"
setsAA[4] = "F"
setsAA[5] = "G"
setsAA[6] = "H"
setsAA[7] = "I"
setsAA[8] = "K"
setsAA[9] = "L"
setsAA[10] = "M"
setsAA[11] = "N"
setsAA[12] = "P"
setsAA[13] = "Q"
setsAA[14] = "R"
setsAA[15] = "S"
setsAA[16] = "T"
setsAA[17] = "V"
setsAA[18] = "W"
setsAA[19] = "Y"
setsAA[20] = "NQ"
setsAA[21] = "GAVLY"
setsAA[22] = "HKR"
setsAA[23] = "STY"
setsAA[24] = "CM"
setsAA[25] = "FGVLAIPMW"
setsAA[26] = "YSNTQC"
setsAA[27] = "LIFWVM"
setsAA[28] = "RKNQPD"
setsAA[29] = "KR"
setsAA[30] = "DE"
setsAA[31] = "GAS"
setsAA[32] = "CDNPT"
setsAA[33] = "EVQH"
setsAA[34] = "MILKR"
setsAA[35] = "FYW"
setsAA[36] = "ND"
setsAA[37] = "QE"

alphabet = {}
alphabet["A"] = 1
alphabet["C"] = 1
alphabet["D"] = 1
alphabet["E"] = 1
alphabet["F"] = 1
alphabet["G"] = 1
alphabet["H"] = 1
alphabet["I"] = 1
alphabet["K"] = 1
alphabet["L"] = 1
alphabet["M"] = 1
alphabet["N"] = 1
alphabet["P"] = 1
alphabet["Q"] = 1
alphabet["R"] = 1
alphabet["S"] = 1
alphabet["T"] = 1
alphabet["V"] = 1
alphabet["W"] = 1
alphabet["Y"] = 1
alphabet["Amide"] = 2
alphabet["Aliphatic"] = 5
alphabet["Basic"] = 3
alphabet["Hydroxyl"] = 3
alphabet["Sulfur"] = 2
alphabet["Non-Polar"] = 9
alphabet["Polar"] = 6.1
alphabet["Hydrophobic"] = 6
alphabet["Hydrophilic"] = 6
alphabet["Pos.Charged"] = 2
alphabet["Neg.Charged"] = 2
alphabet["VerySmall"] = 3.2
alphabet["Small"] = 5.1
alphabet["MediumA"] = 4
alphabet["MediumB"] = 5.1
alphabet["Aromatic"] = 3.1
alphabet["ND"] = 2.2
alphabet["QE"] = 2.2


class Node:
	def __init__(self,aaid,pos):
		self.id = [aaid]
		self.position = [pos]
	def __hash__(self):
		return hash((str(self.id),str(self.position)))
	def __eq__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.id == other.id and self.position == other.position
	def add2Cluster(self,residue):
		self.id.append(residue)
	def getClusters(self):
		cluster = []
		for i in range(0,len(self.id)):
			cluster.append(properties[self.id[i]] + str(self.position[i]))
		return cluster
	def toString(self):
		if len(self.id) > 1:
			return str(self.id)
		return properties[self.id[0]] + str(self.position[0])
	def toStackedString(self):
		if len(self.id) > 1:
			return str(self.id)
		return properties[self.id[0]] + "_" + str(self.position[0])
	def getTuple(self):
		return (self.id[0],self.position[0])
	def getIndex(self):
		return self.position[0]-1

def isHMM():
	for name,sequence in msa.items():
		if '.' in sequence:
			return True
	return False	

def isAA(aa,case):
	if case:
		if aa != "." and aa != "-" and aa.isupper():
			return True
		else:
			return False
	else:
		if aa != "." and aa != "-":
			return True
		else:
			return False

def hmmFiltering():
	new_msa = {}
	hmmValidPositions = []
	seq = msa.items()[0][1]

	for i,aa in enumerate(seq):
		if aa == '-' or aa.isupper():
			hmmValidPositions.append(i)
	for name,seq in msa.items():
		nValidPos = 0
		for i in hmmValidPositions:
			if isAA(seq[i],False):
				nValidPos+=1
		if float(nValidPos)/float(len(hmmValidPositions)) >= minocc:
			new_msa[name] = seq
	return new_msa

def sortedSizeSequences():
	new_seqs = []

	for seqname,sequence in msa.items():
		new_seq = sequence.replace('.','').replace('-','')
		new_seqs.append((seqname,new_seq))

	new_seqs = sorted(new_seqs,key=lambda x: len(x[1]),reverse=True)
	return new_seqs

def maxIdFiltering():
	sortedSeqs = sortedSizeSequences()
	new_msa = {}
	cluster_N = 0
	total = len(sortedSeqs)
	count = 0

	for seqname1,sequence1 in sortedSeqs:
		print("Max identity filtering... (" + str(count) + "/" + str(total) + ")")
		count += 1
		if(cluster_N == 0):
			cluster_N += 1
			new_msa[seqname1] = sequence1 #AInda sem gaps
		else:
			found = False
			for seqname2,sequence2 in new_msa.items():
				similarity = ratio(sequence1,sequence2)
				if similarity >= maxid:
					found = True
					break
			if(not found):
				new_msa[seqname1] = sequence1

	for seqname in new_msa.keys():
		new_msa[seqname] = msa[seqname]
	return new_msa

def writeUnalignedFasta():
	unal_file = outputdir + "/unaligned.fa"
	fw = open(unal_file,'w')
	for seqname,sequence in msa.items():
		fw.write('>' + seqname + "\n")
		fw.write(sequence.replace('.','').replace('-','') + "\n")
	fw.close()
	return unal_file


def maxIdCDhit():
	unal_file = writeUnalignedFasta()
	out_file = outputdir + "/cluster"
	n = 5
	if maxid > 0.7:
		n = 5
	elif maxid > 0.6:
		n = 4
	elif maxid > 0.5:
		n = 3
	else:
		n = 2

	os.system('./cd-hit -i ' + unal_file + ' -o ' + out_file + ' -c ' + str(maxid) + ' -n ' + str(n) + ' -M 3000 -T 2')

	new_msa = {}
	fr = open(out_file)

	for line in fr:
		line = line.strip()
		if len(line) > 1:
			if line[0] == '>':
				seqname = line[1:]
				sequence = msa[seqname]
				new_msa[seqname] = sequence
	return new_msa

	fr.close()


def aa2id(aa):
	aa = aa.upper()
	if aa == 'A':
		return 0
	elif aa == 'C':
		return 1
	elif aa == 'D':
		return 2
	elif aa == 'E':
		return 3
	elif aa == 'F':
		return 4
	elif aa == 'G':
		return 5
	elif aa == 'H':
		return 6
	elif aa == 'I':
		return 7
	elif aa == 'K':
		return 8
	elif aa == 'L':
		return 9
	elif aa == 'M':
		return 10
	elif aa == 'N':
		return 11
	elif aa == 'P':
		return 12
	elif aa == 'Q':
		return 13
	elif aa == 'R':
		return 14
	elif aa == 'S':
		return 15
	elif aa == 'T':
		return 16
	elif aa == 'V':
		return 17
	elif aa == 'W':
		return 18
	elif aa == 'Y':
		return 19
	elif aa == '-':
		return 20
	elif aa == '.':
		return 20
	else:
		return 21

def residueFiltering():
	freqList = []
	N = len(msa)
	seq = msa.items()[0][1]
	removed = [0,0]

	for i,aa in enumerate(seq):
		frequencies = np.zeros(22)
		for sequence in msa.values():
			aa = sequence[i]
			aaid = aa2id(aa)
			frequencies[aaid] += 1.0
		frequencies /= N
		freqList.append(frequencies)

	for seqname,sequence in msa.items():
		new_seq = ""
		for i,aa in enumerate(sequence):
			aaid = aa2id(aa)
			freq = freqList[i][aaid]
			if aa == '.':
				new_seq += '-'
			#elif freq > maxfreq or freq < minfreq:
			elif freq > maxfreq:
				new_seq += '-'
				removed[0] += 1
			elif freq < minfreq:
				new_seq += '-'
				removed[1] += 1
			else:
				new_seq += aa.upper()
		msa[seqname] = new_seq
	return (removed,freqList)

def nodesFiltering(nodes):
	freqList = []
	N = len(msa)
	seq = msa.items()[0][1]
	newNodes = set()
	removed = [0,0]

	for i,aa in enumerate(seq):
		frequencies = np.zeros(38)
		for sequence in msa.values():
			aa = sequence[i].upper()
			if aa in "ACDEFGHIKLMNPQRSTVWY":
				for propid in setsIdList[aa]:
					frequencies[propid] += 1
		frequencies /= N
		freqList.append(frequencies)

	for node in nodes:
		propid = node.id[0]
		position = node.position[0]
		freq = freqList[position-1][propid]
		if freq > maxfreq:
			removed[0] += 1
			continue
		elif freq < minfreq:
			removed[1] += 1
			continue
		else:
			newNodes.add(node)

	return (newNodes,removed,freqList)

def getNodes():
	nodeSet = set()
	if marginal == 1:
		for sequence in msa.values():
			for i,aa in enumerate(sequence):
				if aa in "ACDEFGHIKLMNPQRSTVWY":
					pos = i+1
					for v in setsIdList[aa]:
						node = Node(v,pos)
						nodeSet.add(node)
	else:
		for sequence in msa.values():
			for i,aa in enumerate(sequence):
				if aa in "ACDEFGHIKLMNPQRSTVWY":
					node = Node(setsIdList[aa][0],i+1)
					nodeSet.add(node)

	return nodeSet

def lnfact(x):
	result = np.int(10)
	if x==0 or x==1:
		return 0.0
	else:
		for i in range(1,x+1):
			result = result*i
	return math.log(result)

def stirling(x):
	if x<=8:
		return lnfact(x)
	else:
		return ((x+0.5)*math.log(x))-x+0.918938533205

def lnbdf(N,nx,px):
	if N==nx:
		return N*math.log(px)
	if nx==0:
		return N*math.log(1.0-px)
	else:
		return stirling(N)-stirling(nx)-stirling(N-nx)+nx*math.log(px)+(N-nx)*math.log(1-px)

def eto10(x):
	return x/2.30258509299

def cbd_tietjen(N,n,freq,right):
	a = math.pow(1-freq,N)
	if not right and n== 0:
		return a
	if right and n == 0:
		return 1
	suma = Decimal(0.0)
	if not right:
		suma = Decimal(a)
	u = freq/(1-freq)
	for i in range(2,int(N+2)):
		a = float(a*u*(N+2-i))/float(i-1)
		if (not right and i <= n+1) or (right and i-1 >= n):
			suma += Decimal(a)
	return suma

def cbd(N,n,freq,right):
	if math.pow(1-freq,N) != 0:
		val = cbd_tietjen(N,n,freq,right)
		if val == 0.0:
			return 5e-324
		else:
			return math.ceil(math.log10(cbd_tietjen(N,n,freq,right)))
	suma = 0
	result = 0
	ps = []
	minP = eto10(lnbdf(N,n,freq))
	if right:
		for i in range(n,N+1):
			val = eto10(lnbdf(N,i,freq))
			ps.append(val)
			if val > minP:
				minP = val
		for val in ps:
			suma += math.pow(10,val)-math.floor(minP)
	else:
		for i in range(0,n+1):
			val = eto10(lnbdf(N,i,freq))
			ps.append(val)
			if val > minP:
				minP = val
		for val in ps:
			suma += math.pow(10,val)-math.floor(minP)
	#print(suma)
	return math.floor(math.log(suma))+math.floor(minP)

def DRCN(res1,res2):
	aa1_list = setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	Na = 0.0
	Nb = 0.0
	Nba = 0.0
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0
	n = float(len(msa))

	for sequence in msa.values():
		if sequence[i1].upper() in aa1_list:
			Na +=1.0
			if sequence[i2].upper() in aa2_list:
				Nb += 1.0
				Nba += 1.0
				a += 1.0
			else:
				b += 1.0
		elif sequence[i2].upper() in aa2_list:
			Nb += 1.0
			c += 1.0
		else:
			d += 1.0

	jc = a/(a+b+c)

	if Nba > (Na*(Nb/n)):
		freq1 = float(Nb)/float(n)
		freq2 = float(Na)/float(n)
		pv1 = cbd(Na,Nba,freq1,True)*-1
		pv2 = cbd(Nb,Nba,freq2,True)*-1

		if pv1 > pv2:
			return (pv1,jc)
		else:
			return (pv2,jc)
	else:
		#Anti-Correlation
		return -1.0

def pearsonEdge(res1,res2):
	aa1_list = setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	sx = 0.0
	sy = 0.0
	sxx = 0.0
	syy = 0.0
	sxy = 0.0
	n = float(len(msa))

	for sequence in msa.values():
		x = 0.0
		y = 0.0

		if sequence[i1].upper() in aa1_list:
			x=1.0
		if sequence[i2].upper() in aa2_list:
			y=1.0
		sx += x
		sy += y
		sxx += x * x
		syy += y * y
		sxy += x * y

	#print(str(aa1_list) + " " + str(aa2_list) + " " + str(sx) + " " + str(sy) + " " + str(sxx) + " " + str(syy) + " " + str(sxy) + " " + str(n))
	cov = sxy / n - sx * sy / n / n
	sigmax = math.sqrt(sxx / n -  sx * sx / n / n)
	sigmay = math.sqrt(syy / n -  sy * sy / n / n)
	return cov / sigmax / sigmay

def borgatti(res1,res2,type):
	aa1_list = setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0

	for sequence in msa.values():
		bool1 = sequence[i1].upper() in aa1_list
		bool2 = sequence[i2].upper() in aa2_list

		if bool1 and bool2:
			a+=1.0
		elif bool1:
			b += 1.0
		elif bool2:
			c += 1.0
		else:
			d += 1.0

	if type == 0:#BN
		if a+b > a+c:
			return a/(a+c)
		else:
			return a/(a+b)
	elif type == 1:#JC
		return a/(a+b+c)
	elif type == 2:#Bonacich
		if a*d == b*c:
			return 0.5
		else:
			return (a*d - math.sqrt(a*b*c*d))/((a*d)-(b*c))

def survivorFixed(k,M,n,N):
	quant, tot, good, draw = k, M, n, N
	k2 = np.arange(quant + 1, draw + 1)
	logpmf = hypergeom._logpmf(k2, tot, good, draw)   # Evaluate the log-pmf instead
	logsf = logsumexp(logpmf)
	return logsf

def getTumminelloTuple(res1,res2):
	aa1_list = setsAA[res1.id[0]]
	i1 = res1.position[0]-1
	aa2_list = setsAA[res2.id[0]]
	i2 = res2.position[0]-1
	wProj = 0
	Di = 0
	Dj = 0
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0

	for sequence in msa.values():
		if sequence[i1].upper() in aa1_list:
			Di +=1.0
			if sequence[i2].upper() in aa2_list:
				Dj += 1
				wProj += 1
				a += 1.0
			else:
				b += 1.0
		elif sequence[i2].upper() in aa2_list:
			Dj += 1
			c += 1.0
		else:
			d += 1.0

	jc = a/(a+b+c)

	if Di > Dj:
		return ((wProj,Di,Dj),jc)
	else:
		return ((wProj,Dj,Di),jc)


def genCorrelationNetwork(method,netthr,freqList):
	N = len(nodes)
	nodeList = list(nodes)
	network = []
	tummineloDic = {}
	#temp = 0

	for i in range(0,N-1):
		print("Calculating Correlations... (" + str(i) + "/" + str(N) + ")")
		for j in range(i+1,N):
			n1 = nodeList[i]#Object
			n2 = nodeList[j]#Object
			jc = 0.0

			fr1 = freqList[n1.getIndex()][n1.id][0]
			fr2 = freqList[n2.getIndex()][n2.id][0]
			if n1.position != n2.position and abs(fr1-fr2) <= FREQ_WINDOW:
				#temp+=1
				if method == 1:
					w,jc = DRCN(n1,n2)
				else:		
					w = 0.0
					tummData,jc = getTumminelloTuple(n1,n2)
					if tummData in tummineloDic:
						w = tummineloDic[tummData]
					else:
						w = survivorFixed(tummData[0]-1,len(msa),tummData[1],tummData[2])*-1
						tummineloDic[tummData] = w

				if w >= netthr:
					#print(str(temp) + " " + n1.toString() + " " + n2.toString() + " " + str(w))
					network.append((n1,n2,w,jc))

	network = sorted(network,key=lambda x: x[2],reverse=True)
	if len(network) > 100000:
		network = network[0:100000]
	return network

def getNodeFrequency(node):
	aa_list = setsAA[node.id[0]]
	i = node.position[0]-1
	freq = 0.0
	for sequence in msa.values():
		if sequence[i] in aa_list:
			freq += 1.0
	return float(freq)/float(len(msa))

def pickColor(commIndex):
	colors = ["#9180ff","#7fffa1","#fffbbf","#ff80b3","#bffff2","#ffee00","#c8bfff","#ff2200","#4100f2","#81f200","#ffaa00","#ffbfbf","#3d9df2","#917399","#992654","#822699","#94994d","#258c00","#8c0000","#8c5e00"]
	if commIndex < len(colors):
		return colors[commIndex]
	else:
		r = lambda: random.randint(0,255)
		return '#%02X%02X%02X' % (r(),r(),r())

def filterNodesPriori(Gnx):
	for Ni in Gnx.nodes():
		if Ni in Gnx.nodes():
			neighbors = Gnx.neighbors(Ni)
			neighDic = {}
			neighDic[Ni.position[0]] = [Ni.id[0]]

			for node in neighbors:
				aaid = node.id[0]
				pos = node.position[0]
				if pos in neighDic:
					neighDic[pos].append(aaid)
				else:
					neighDic[pos] = [aaid]

			for pos, aaidList in neighDic.items():
				if len(aaidList) > 1:
					minsize = 999
					bestaaid = []
					for aaid in aaidList:
						value = len(setsAA[aaid])
						if value < minsize:
							bestaaid = [aaid]
							minsize = value
						elif value == minsize:
							bestaaid.append(aaid)
					for aaid in aaidList:
						if aaid not in bestaaid:
							node = Node(aaid,pos)
							Gnx.remove_node(node)
	return Gnx

#Usando Avg
def consineDistance(Ni,Nj):
	comm1 = Ni.getClusters()
	comm2 = Nj.getClusters()
	col1 = df[comm1[0]]
	col2 = df[comm2[0]]
	for res in comm1[1:]:
		col1+=df[res]
	for res in comm2[1:]:
		col2+=df[res]
	col1 /= len(comm1)
	col2 /= len(comm2)

	return scipy.spatial.distance.cosine(col1,col2)
	#return cc

def merge_nodes(G,nodes, new_node, attr_dict=None, **attr):
	G.add_node(new_node, attr_dict, **attr) # Add the 'merged' node

	for n1,n2,data in G.edges(data=True):
		if n1 in nodes:
			G.add_edge(new_node,n2,data)
		elif n2 in nodes:
			G.add_edge(n1,new_node,data)

	for n in nodes: # remove the merged nodes
		G.remove_node(n)

def communityDetection(Gnx,cutoff):
	communities = []
	count = 0
	distances = {}
	while True:
		count += 1
		mindist = 99999999
		#maxdist = 0
		mergecomm1 = ""
		mergecomm2 = ""
		
		for Ni,Nj in Gnx.edges():
			#dist = euclidianDistance(Ni,Nj)
			tupla = (Ni.toString(),Nj.toString())
			if tupla in distances:
				dist = distances[tupla]
			else:
				dist = consineDistance(Ni,Nj)
				distances[tupla] = dist
			if dist < mindist:
			#if dist > maxdist:
				mindist = dist
				#maxdist = dist
				mergecomm1 = Ni
				mergecomm2 = Nj
		if mergecomm1 == "" or mergecomm2 == "":
			break
		if mindist > cutoff:
			break
		#print(str(len(Gnx)) + " nodes. Distance = " + str(mindist))
		#new_node = mergecomm1 + "_" + mergecomm2
		new_node = Node(-1,-1)
		new_node.id = mergecomm1.id + mergecomm2.id
		new_node.position = mergecomm1.position + mergecomm2.position
		merge_nodes(Gnx,[mergecomm1,mergecomm2],new_node)

	#fw = open(outputcomm,'w')
	count = 0
	for comms in Gnx.nodes():
		idlist = comms.id
		poslist = comms.position
		if len(idlist) > 1:
			comm = []
			for i in range(len(idlist)):
				res = Node(idlist[i],poslist[i])
				#res = properties[idlist[i]] + str(poslist[i])
				#print(res + " "),
				#fw.write(node + " ")
				comm.append(res)
				count +=1
			communities.append(comm)
			#print("\n")
			#fw.write("\n")
	return count,communities

def communityDetection2(Gnx):
	communities = []
	count = 0
	components = nx.connected_component_subgraphs(Gnx)
	jc_all = 0.0
	N = 0.0

	for comp in components:
		Nedges = len(comp.edges())
		if Nedges > 0:
			jc_comm = 0.0
			for ni,nj in comp.edges():
				jc_comm += Gnx[ni][nj]['weight']
			jc_comm /= Nedges
			if jc_comm >= MIN_JC:
				communities.append(comp.nodes())
				N += 1.0
				jc_all += jc_comm
	if N == 0.0:
		return 0.0,[],0.0
	jc_all /= N

	return N,communities,jc_all

def filterRedundancy(comm):
	new_comm = set()
	positions = {}

	for node in comm:
		aa2_id,pos = node.getTuple()
		if pos in positions:
			aa_id = positions[pos]
			aa = properties[aa_id]
			aa2 = properties[aa2_id]
			size1 = alphabet[aa]
			size2 = alphabet[aa2]
			if size2 < size1:
				positions[pos] = aa2_id
		else:
			positions[pos] = aa2_id

	for pos,aa_id in positions.items():
		residue = Node(aa_id,pos)
		new_comm.add(residue)
	return new_comm

def filterNetwork(Gnx,comms):
	comm_nodes = set()

	for comm in comms:
		comm_nodes = comm_nodes.union(comm)

	for node in Gnx.nodes():
		if node not in comm_nodes:
			Gnx.remove_node(node)

	return Gnx

def writeCommunities(path,communities):
	communities = sorted(communities,key=lambda x: len(x),reverse=True)
	fw1 = open(path,'w')
	for comm in communities:
		if len(comm) > 1:
			for node in comm:
				fw1.write(node.toStackedString() + " ")
			fw1.write("\n")
	fw1.close()

def writeBackbone(path,backbone):
	avgJC = 0.0
	Nedges = 0.0
	fw2 = open(path,'w')
	for ni,nj in backbone.edges():
		w = backbone[ni][nj]['weight']
		pv = backbone[ni][nj]['pvalue']
		fw2.write(ni.toStackedString() + " " + nj.toStackedString() + " " + str(w)  + " " + str(pv) + "\n")
		avgJC+=w
		Nedges+=1.0
	fw2.close()
	return avgJC/Nedges

def readFasta(inputfile):
	msa = {}
	fr = open(inputfile,'r')
	header = None
	sequence = ""

	for line in fr:
		line = line.strip()
		if line.startswith('>'):
			if header is not None:
				msa[header] = sequence
				header = line[1:]
				sequence = ""
			else:
				header = line[1:]
		else:
			sequence += line
	if header is not None:
		msa[header] = sequence

	fr.close()
	return msa

def readStockholm(inputfile):
	msa = {}
	fr = open(inputfile,'r')

	for line in fr:
		if line[0] != "#":
			line = line.strip()
			temp = line.split()
			if len(temp) > 1:
				seqname = temp[0]
				sequence = temp[1]
				msa[seqname] = sequence

	fr.close()
	return msa

def printHelp():
	print("###########################################################################\n"
		"CONAN - Co-Variation Network Analyzer\n"
		"###########################################################################\n\n"
		"This software analyzes a residue co-variation network with the aim to emphasize local evolutionary constraints and "
		"consequently detect functional and specificity determinant sites. Its only mandatory input is a multiple sequence "
		"alignment in one of the following formats: Pfam, Stockholm or FASTA; and the output consists of several networks and "
		"communities files from several cuts of the network.\n\n"
		"The mandatory inputs consits of:\n"
		"-i <filename> - A multiple sequence alignment file\n"
		"-o <directory> - An output directory path\n\n"
		"Optional parameters:\n"
		"-p <value> - Minimum correlation p-value (in -log(x))\n"
		"-O <value> - Minimum Occupancy value. It is used to remove fragments of sequences. (0-1)\n"
		"-I <value> - Maximum Identity value. It is used to remove high identity sequences. (0-1)\n"
		"-f <value> - Minimum node frequency. It removes nodes (residues) that rarely occurs in the alignment. (0-1)\n"
		"-F <value> - Maximum node frequency. It removes higly conserved nodes (residues). (0-1)\n"
		"-m <0 or 1> - Method to Statistically validate the nodes.\n\t0 - Tumminello (Network based validation)\n\t1 - DRCN (Frequency based validation)\n"
		"-e <0 or 1> - Include marginally conservation properties.\n"
		"\t0 - Consider only co-variation between amino acids.\n"
		"\t1 - Also include stereochemical and structural amino acids properties.\n"
		"\nMore information at http://www.biocomp.icb.ufmg.br/conan")

######GENERAL VARIABLES##########
hmm = False

######PARSE PARAMETERS##########
if len(sys.argv) < 2:
   printHelp()
   sys.exit()

argv = sys.argv[1:]

inputfile = ''
outputdir = ''
minocc = 0.8
maxid = 0.8
minfreq = 0.05
maxfreq = 0.8
method = 5
netthr = 5
marginal = 1
min_pv = 15

try:
	opts, args = getopt.getopt(argv,"i:o:p:O:I:f:F:m:e:h")
	#print(opts)
	#print(args)
except getopt.GetoptError:
	print("###########################################################################\n"
		"CONAN - Co-Variation Network Analyzer\n"
		"###########################################################################\n\n"
		"This software analyzes a residue co-variation network with the aim to emphasize local evolutionary constraints and "
		"consequently detect functional and specificity determinant sites. Its only mandatory input is a multiple sequence "
		"alignment in one of the following formats: Pfam, Stockholm or FASTA; and the output consists of several networks and "
		"communities files from several cuts of the network.\n\n"
		"This software was still not pubblished, but further information about this methodology can be accessed at:\n\n"
		"da Fonseca, N. J., Afonso, M. Q. L., de Oliveira, L. C., & Bleicher, L. (2018). A new method bridging graph "
		"theory and residue co-evolutionary networks for specificity determinant positions detection. Bioinformatics.\n\n"
		"The mandatory inputs consits of:\n"
		"-i <filename> - A multiple sequence alignment file\n"
		"-o <directory> - An output directory path\n"
		"Use -h for access the optional parameters.\n"
		"You can also find more information at http://www.biocomp.icb.ufmg.br/conan")
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		printHelp()
		sys.exit()
	elif opt in ("-i"):
		inputfile = arg
	elif opt in ("-o"):
		outputdir = arg
	elif opt in ("-p"):
		min_pv = float(arg)
	elif opt in ("-O"):
		minocc = float(arg)
	elif opt in ("-I"):
		maxid = float(arg)
	elif opt in ("-f"):
		minfreq = float(arg)
	elif opt in ("-F"):
		maxfreq = float(arg)
	elif opt in ("-m"):
		method = int(arg)
	elif opt in ("-e"):
		marginal = int(arg)
if inputfile == "":
	print("The input file is a mandatory parameter.")
	sys.exit()
if outputdir == "":
	print("The output directory is a mandatory parameter.")
	sys.exit()

###########READ ALIGNMENT#############
##ATUALIZAR PARA LER FASTA
msa = {}
nodes = set()
nodes_frequencies = []

###Check MSA type
fr = open(inputfile,'r')

line = fr.readline()
fr.close()

if line[0] == '>':
	msa = readFasta(inputfile)
else:
	msa = readStockholm(inputfile)

N_msa = len(msa)
#########CREATE FOLDER############
if outputdir[len(outputdir)-1] == '/':#Remove /
	outputdir = outputdir[:-1]
if not os.path.exists(outputdir):
	os.makedirs(outputdir)

######FILTER MSA BY OCC############
print("Filtering the Alignment...")
hmm = isHMM()
if hmm and (minocc > 0.0):
	msa = hmmFiltering()
N_msaOCC = len(msa)
#print(N_msaOCC)

######FILTER BY MAXID#############
if maxid > 0.0:
	msa = maxIdCDhit()
N_msaOCC_ID = len(msa)

######FILTER RESIDUES###############
if marginal == 0:
	res_rem,nodes_frequencies = residueFiltering()

#######WRITE FILTERED MSA#########
fw = open(outputdir + "/filtered.txt",'w')
for name,sequence in msa.items():
	fw.write(name + "\t" + sequence + "\n")
fw.close()

######GENERATE NETWORK############
print("Calculating the correlations...")
nodes = getNodes()

if marginal == 1:
	nodes,res_rem,nodes_frequencies = nodesFiltering(nodes)

network = genCorrelationNetwork(method,netthr,nodes_frequencies)

#####COMMUNITY DETECTION##########
print("Detecting communities...")
#outputcomm = outputdir + "/comms.txt"

comPath = outputdir + "/communities/"
if not os.path.exists(comPath):
	os.makedirs(comPath)
backPath = outputdir + "/backbones/"
if not os.path.exists(backPath):
	os.makedirs(backPath)

fw = open(outputdir + "/cutoff.txt",'w')

#Filter by p-value
network = [(ni,nj,w,jc) for (ni,nj,w,jc) in network if w >= min_pv]
#Filter by Min Jacard Coeficient Value (Reduce complexity)
network = [(ni,nj,w,jc) for (ni,nj,w,jc) in network if jc >= MIN_JC]
#Sort by Jacard Coeficient
network = sorted(network, key=lambda tup: tup[3])

communities = []
detected_nodes = {}
Gnx = nx.Graph()

###Last State Variables###
lsv_communities = []
lsv_N_residues = 0

#First Iteration
if len(network) > 0:
	ni,nj,pv,jc = network.pop()
	comm = set([ni,nj])
	detected_nodes[ni] = 0
	detected_nodes[nj] = 0
	communities.append(comm)
	lsv_communities = communities
	Gnx.add_edge(ni,nj,weight=jc,pvalue=pv)
	lsv_N_residues = 2

#Post iterations
while(len(network) > 0):
	ni,nj,pv,jc = network.pop()
	Gnx.add_edge(ni,nj,weight=jc,pvalue=pv)
	if ni in detected_nodes and nj in detected_nodes:
		commI = detected_nodes[ni]
		commJ = detected_nodes[nj]
		if commI != commJ:#Merge communities
			if commI > commJ:
				temp = commI
				commI = commJ
				commJ = temp
			communities[commI] = communities[commI].union(communities[commJ])
			del communities[commJ]
			for res,comm in detected_nodes.items():
				if comm == commJ:
					detected_nodes[res] = commI
				elif comm > commJ:
					detected_nodes[res] = comm-1
		else:
			continue #Both nodes were already detected
	if ni in detected_nodes:#Add Nj to Ni community
		commI = detected_nodes[ni]
		communities[commI].add(nj)
		detected_nodes[nj] = commI
	elif nj in detected_nodes:#Add Ni to Nj community
		commJ = detected_nodes[nj]
		communities[commJ].add(ni)
		detected_nodes[ni] = commJ
	else:#New Cluster
		comm = set([ni,nj])
		i = len(communities)
		detected_nodes[ni] = i
		detected_nodes[nj] = i
		communities.append(comm)

	#Validate and Write State
	N_residues = 0
	comms = []
	for i,comm in enumerate(communities):
		fixed_comm = filterRedundancy(comm)
		N_residues+=len(fixed_comm)
		comms.append(fixed_comm)
	if N_residues > lsv_N_residues:#Write lsv
		Gnx = filterNetwork(Gnx,lsv_communities)

		####WRITE EVERYTHING
		writeCommunities(comPath + str(lsv_N_residues),lsv_communities)
		jc_avg = writeBackbone(backPath + str(lsv_N_residues),Gnx)
		Ncoms = len(lsv_communities)
		fw.write(str(lsv_N_residues) + " " + str(Ncoms) + " " + str(jc_avg) + "\n")
		print("N. Residues: " + str(lsv_N_residues) + "\tN. communities: " + str(Ncoms) + " JC avg: " + str(jc_avg))
		lsv_N_residues = N_residues
	else:#Store in lsv
		lsv_communities = comms
	