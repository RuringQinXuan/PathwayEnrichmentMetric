import os
from xml.dom import minidom

# class CTDPathway:
	# def __init__(self,pathway_name,pathway_ID,pathway_P,pathway_CorrectedP):
		# self.pathway_name=pathway_name
		# self.pathway_ID=pathway_ID
		# self.pathway_P=pathway_P
		# self.pathway_CorrectedP=pathway_CorrectedP

# class Get_CTD_pathway():
	# def __init__(self,train_fold,drug):
		# self.train_fold=train_fold
		# self.drug=drug
	# def GetElementsContent(self,name,mydoc):
		# ElementsContents=mydoc.getElementsByTagName(name)
		# contentList=[]
		# for item in ElementsContents:
			# contentList.append(item.firstChild.data)
		# return contentList
	# def Get_trainFileDirList(self):
		# testFileDirList=[]
		# drugList=[]
		# for root, dirs, files in os.walk(self.train_fold):  
			# for file in files:  
				# if os.path.splitext(file)[1] == '.xml' and self.drug in file:  
					# drug_file_dir=os.path.join(root, file)
					# drug_dic_name=os.path.splitext(file)[0]
		# return drug_file_dir,drug_dic_name
	# def Get_pathway_p(self,drug_file_dir,drug_dic_name):
		# filedir=drug_file_dir
		# drugname=drug_dic_name
		# mydoc = minidom.parse(filedir)
		# pathway_names=self.GetElementsContent('Pathway',mydoc)
		# pathway_IDs=self.GetElementsContent('PathwayId',mydoc)
		# pathway_Ps=self.GetElementsContent('P-value',mydoc)
		# pathway_CorrectedPs=self.GetElementsContent('CorrectedP-value',mydoc)
		# dic_pathway={}
		# for i in range(len(pathway_IDs)):
			# if 'KEGG' in pathway_IDs[i]:
				# pathway_ID=pathway_IDs[i].replace("KEGG","path")
				# dic_pathway[pathway_ID]=pathway_CorrectedPs[i]
		# return drugname,dic_pathway
	# def Get_all_pathway(self,testFileDirList,drugList):
		# dic_drug_pathway={}
		# for j in range(len(testFileDirList)):
			# if self.drug != None and self.drug in testFileDirList[j]:
				# drugname,dic_pathway=self.Get_pathway_p(j,testFileDirList,drugList)
				# return dic_pathway
			# else:
				# drugname,dic_pathway=self.Get_pathway_p(j,testFileDirList,drugList)
				# dic_drug_pathway[drugname]=dic_pathway
				# return dic_drug_pathway
	
def GetElementsContent(name,mydoc):
	ElementsContents=mydoc.getElementsByTagName(name)
	contentList=[]
	for item in ElementsContents:
		contentList.append(item.firstChild.data)
	return contentList	
	

def get_drugs(drugs_file_dir):
	drugs=[]
	with open(drugs_file_dir,'r') as drugs_file:
		for line in drugs_file:
			line= line.strip()
			drugs.append(line)
	return drugs
	
def get_CTD_pathway_dic_from_xml(drugs):
	CTDdic={}
	for drug in drugs:
		dic_pathway={}
		CTD_pathwayfile_dir=CTD_fold+'/CTD_'+drug+'_pathways.xml'
		mydoc = minidom.parse(CTD_pathwayfile_dir)
		pathway_names=GetElementsContent('Pathway',mydoc)
		pathway_IDs=GetElementsContent('PathwayId',mydoc)
		pathway_CorrectedPs=GetElementsContent('CorrectedP-value',mydoc)
		for i in range(len(pathway_IDs)):
			if 'KEGG' in pathway_IDs[i]:
				pathway_ID=pathway_IDs[i].replace("KEGG","path")
				dic_pathway[pathway_ID]=pathway_CorrectedPs[i]
		CTDdic[drug]=dic_pathway
	with open('./data/pathway/CTD/CTDdic','w') as file:
		file.write(str(CTDdic))
																										
CTD_fold="./data/pathway/CTD"
drugs_file_dir='./data/drug/drug_name'
drugs=get_drugs(drugs_file_dir)
for drug in drugs:
	get_CTD_pathway_dic_from_xml(drugs)

