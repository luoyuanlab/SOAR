### Usage: python3 geo-query.py
### Author: Yiming Li

import numpy as np
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time, os, shutil, sys

def fetch_species_GDS(species):
	# species = "mouse"
	urls = ["http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=spatial+transcriptomics+AND+" + species + "[organism]&retmax=100000&usehistory=y", "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=spatial+transcriptome+AND+" + species + "[organism]&retmax=100000&usehistory=y", "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=spatial+RNA-seq+AND+" + species + "[organism]&retmax=100000&usehistory=y", "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=spatial+RNA+sequencing+AND+" + species + "[organism]&retmax=100000&usehistory=y"]
	items = []
	query_keys = []
	WebEnvs = []
	for url in urls:
		resp = requests.get(url)
		with open('tmp2.xml', 'wb') as f:
			f.write(resp.content)
		
		with open('tmp2.xml', 'r') as file:
			xml_text = file.read()
		
		if "API rate limit exceeded" in xml_text:
			sys.exit('[ERROR] E-utils API rate limit exceeded')
		
		tree = ET.parse('tmp2.xml')
		root = tree.getroot()
		query_keys.append(root.findall('./QueryKey')[0].text)
		WebEnvs.append(root.findall('./WebEnv')[0].text)
		for item in root.findall('./IdList/Id'):
			items.append(item.text)
	
	return(items, query_keys, WebEnvs)

def fetch_PMIDs(query_keys, WebEnvs):
	pmids = []
	items = []
	for query_key, WebEnv in zip(query_keys, WebEnvs):
		url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gds&db=pubmed&query_key=' + str(query_key) + '&WebEnv=' + str(WebEnv)
		resp = requests.get(url)
		with open('tmp3.xml', 'wb') as f:
			f.write(resp.content)
		
		with open('tmp3.xml', 'r') as file:
			xml_text = file.read()
		
		if "API rate limit exceeded" in xml_text:
			sys.exit('[ERROR] E-utils API rate limit exceeded')
		
		tree = ET.parse('tmp3.xml')
		root = tree.getroot()
		
		for item in root.findall('./LinkSet/IdList/Id'):
			items.append(item.text)
		
		for pmid in root.findall('./LinkSet/LinkSetDb/Link/Id'):
			pmids.append(pmid.text)
	
	return(items, pmids)

def loadRSS(gds_id):
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=' + str(gds_id)
	resp = requests.get(url)
	with open('tmp.xml', 'wb') as f:
		f.write(resp.content)

def parseXML(xmlfile):
	with open(xmlfile, 'r') as file:
		xml_text = file.read()
	
	if "API rate limit exceeded" in xml_text:
		sys.exit('[ERROR] E-utils API rate limit exceeded')
	
	tree = ET.parse(xmlfile)
	root = tree.getroot()
	items = []
	for item in root.findall('./DocSum/Item'):
		if (item.attrib['Name'] in ["Accession", "title", "summary", "GPL", "GSE", "taxon", "gdsType", "FTPLink"]):
			items.append(item.text)
	
	return(items)

def fetch_species_meta(species, save_dir):
	print("### Starting initial query for [" + species + "]")
	species_ids, query_keys, WebEnvs = fetch_species_GDS(species)
	query_keys_df = pd.DataFrame({"query_key": query_keys, "WebEnv": WebEnvs})
	query_keys_df.to_csv(save_dir + "/" + species + "_query_keys.tsv", sep='\t', index = False)
	
	current_species_ids = [line.strip() for line in open(species + "_GDS_current.txt", 'r')]
	species_ids = list(set(species_ids) - set(current_species_ids))
	species_ids_concat = list(set(species_ids) | set(current_species_ids))
	print("### Initial query for [" + species + "] completed\n### Number of species GDS IDs: " + str(len(current_species_ids)) + " (previous), " + str(len(species_ids)) + " (this query), " + str(len(species_ids_concat)) + " (concatenated)")
	
	### Get meta-information
	print("### Starting query for [" + species + "] meta-information")
	results = []
	for gds_id in species_ids:
		loadRSS(gds_id)
		time.sleep(1)
		results.append(parseXML('tmp.xml'))
	
	results = pd.DataFrame(results, columns=["Accession", "title", "summary", "GPL", "GSE", "taxon", "gdsType", "FTPLink"])
	results["GDS_ID"] = species_ids
	results["Organism"] = species
	results.to_csv(save_dir + "/" + species + ".tsv", sep='\t', index = False)
	print("### Meta-information of [" + species + "] GDS IDs saved to " + save_dir + "/" + species + ".tsv")
	
	### Save GDS lists
	with open(save_dir + "/" + species + "_GDS_" + save_dir + ".txt", 'w') as f:
		for line in species_ids:
			f.write(f"{line}\n")
	
	shutil.copyfile(species + "_GDS_current.txt", species + "_GDS_current.txt.bk")
	
	with open(species + "_GDS_current.txt", 'w') as f:
		for line in species_ids_concat:
			f.write(f"{line}\n")
	
	print("### " + species + "_GDS_current.txt overwritten")
	
	print("### Starting query for [" + species + "] PMIDs")
	gds_ids, pmids = fetch_PMIDs(query_keys, WebEnvs)
	pmids_df = pd.DataFrame({"PMID": pmids})
	pmids_df.to_csv(save_dir + "/" + species + "_PMIDs.tsv", sep='\t', index = False)
	print("### PMIDs associated with [" + species + "] saved to " + save_dir + "/" + species + "_PMIDs.tsv")
	
	print("### [" + species + "] completed\n")
	return(results)

##########################################################################################

owd = os.getcwd()
os.chdir("/share/fsmresfiles/SpatialT/GEO_query")
save_dir = time.strftime("%Y%m%d")
os.makedirs(save_dir, exist_ok = True)

### Get mouse GDS IDs
results_mouse = fetch_species_meta("mouse", save_dir)

### Get human GDS IDs
results_human = fetch_species_meta("human", save_dir)

### Combine results
results_mouse = pd.read_csv(save_dir + "/mouse.tsv", sep = '\t')
results_human = pd.read_csv(save_dir + "/human.tsv", sep = '\t')
results_all = pd.concat([results_mouse, results_human])
results_all = results_all.sort_values(by=['GDS_ID', 'Organism'])

results_all['Accession_is_GSM'] = results_all['Accession'].str.startswith('GSM', na=False)
results_all['GSM'] = np.where(results_all['Accession_is_GSM'], results_all['Accession'], "")
results_all['Accession'] = np.where(results_all['Accession_is_GSM'], "", results_all['Accession'])
results_all['Is ST data'] = ""
results_all['Technology'] = ""
results_all['Platform'] = ""
results_all['Add to SOAR'] = ""
results_all['PMID'] = ""
results_all['GSE'] = "GSE" + results_all['GSE'].apply(str)

reordered_columns = ['GDS_ID', 'Organism', 'Accession', "GSE", "GSM", 'Is ST data', 'Technology', 'Platform', 'PMID', 'Add to SOAR', 'title', 'summary', 'GPL', 'taxon', 'gdsType', 'FTPLink']
results_all = results_all[reordered_columns]
results_all = results_all.sort_values(by=['GSE'])
results_all.to_csv(save_dir + "/all.tsv", sep='\t', index = False)

### Clean up
os.remove("tmp.xml")
os.remove("tmp2.xml")
os.remove("tmp3.xml")
os.chdir(owd)
