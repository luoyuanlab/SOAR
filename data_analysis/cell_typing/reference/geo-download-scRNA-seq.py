"""
This script queries the GEO using esearch and esummary for scRNA-seq datasets related to a certain organ/species
	* Note that this script does NOT check the validity of the arguments
	* Two intermediate files will be created during runtime: esearch.xml and esummary.xml (cleaned up at the end)
	* Please change any space in the keyword to a plus sign, e.g. spinal+cord instead of spinal cord

Usage: python geo-download-scRNA-seq.py <organ> <species>
	e.g. python geo-download-scRNA-seq.py lymph+node mouse

Author:
	Yiming Li
"""

import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time
import sys
import os

def loadRSS_esearch(organ, species, keyword):
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=' + keyword + '+' + organ + '+AND+' + species + '[organism]&retmax=100000&usehistory=y'
	resp = requests.get(url)
	with open('esearch.xml', 'wb') as f:
		f.write(resp.content)

def parseXML_esearch(xmlfile):
	tree = ET.parse(xmlfile)
	root = tree.getroot()
	items = []
	for item in root.findall('./IdList/Id'):
		items.append(item.text)
	
	return(items)

def loadRSS_esummary(gds_id):
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=' + str(gds_id)
	resp = requests.get(url)
	with open('esummary.xml', 'wb') as f:
		f.write(resp.content)

def parseXML_esummary(xmlfile):
	tree = ET.parse(xmlfile)
	root = tree.getroot()
	items = []
	for item in root.findall('./DocSum/Item'):
		if (item.attrib['Name'] in ["Accession", "title", "summary", "GPL", "taxon", "gdsType", "FTPLink"]):
			items.append(item.text)
	
	return(items)

organ = sys.argv[1] # organ = 'spinal+cord'
species = sys.argv[2] # species = "human"



############



keywords = ['scRNA-seq', 'single+cell+RNA-seq', 'single+cell+RNA+sequencing', 'single+cell+transcriptomics', 'single+cell+transcriptome'] # Change this to make the search broader/narrower

### Step 1. Get GDS IDs
gds_ids = []
for keyword in keywords:
	loadRSS_esearch(organ, species, keyword)
	time.sleep(1)
	gds_ids = gds_ids + parseXML_esearch('esearch.xml')

gds_ids = list(set(gds_ids)) # Remove duplicates

print(">>>>> [Organ: " + organ + "; Species: " + species + "] " + str(len(gds_ids)) + " GDS IDs found <<<<<\n")

### Step 2. Get meta-information
results = []
for gds_id in gds_ids:
	loadRSS_esummary(gds_id)
	time.sleep(1)
	results.append(parseXML_esummary('esummary.xml'))
	print("[" + gds_id + "] completed")

results = pd.DataFrame(results, columns=["Accession", "title", "summary", "GPL", "taxon", "gdsType", "FTPLink"])
results["GDS_ID"] = gds_ids
results["Species"] = species
results["Organ"] = organ
results.to_csv(species + "-" + organ + ".csv", sep='\t', index = False)

os.remove("esearch.xml")
os.remove("esummary.xml")
