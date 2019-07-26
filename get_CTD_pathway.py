import os
from xml.dom import minidom

def write_dic_to_file(file_dir,dic):
	with open(file_dir,'w') as file:
		file.write(str(dic))

class CTDPathway:
	def __init__(self,pathway_name,pathway_ID,pathway_P,pathway_CorrectedP):
		self.pathway_name=pathway_name
		self.pathway_ID=pathway_ID
		self.pathway_P=pathway_P
		self.pathway_CorrectedP=pathway_CorrectedP

class Get_CTD_pathway():
	def __init__(self,train_fold,drug):
		self.train_fold=train_fold
		self.drug=drug
	def GetElementsContent(self,name,mydoc):
		ElementsContents=mydoc.getElementsByTagName(name)
		contentList=[]
		for item in ElementsContents:
			contentList.append(item.firstChild.data)
		return contentList
	def Get_trainFileDirList(self):
		testFileDirList=[]
		drugList=[]
		for root, dirs, files in os.walk(self.train_fold):  
			for file in files:  
				if os.path.splitext(file)[1] == '.xml' and self.drug in file:  
					drug_file_dir=os.path.join(root, file)
					drug_dic_name=os.path.splitext(file)[0]
		return drug_file_dir,drug_dic_name
	def Get_pathway_p(self,drug_file_dir,drug_dic_name):
		filedir=drug_file_dir
		drugname=drug_dic_name
		mydoc = minidom.parse(filedir)
		pathway_names=self.GetElementsContent('Pathway',mydoc)
		pathway_IDs=self.GetElementsContent('PathwayId',mydoc)
		pathway_Ps=self.GetElementsContent('P-value',mydoc)
		pathway_CorrectedPs=self.GetElementsContent('CorrectedP-value',mydoc)
		dic_pathway={}
		for i in range(len(pathway_IDs)):
			if 'KEGG' in pathway_IDs[i]:
				pathway_ID=pathway_IDs[i].replace("KEGG","path")
				dic_pathway[pathway_ID]=pathway_CorrectedPs[i]
		return drugname,dic_pathway
	
train_fold="./data/pathway/CTD"	
# from Get_CTD_pathway import Get_CTD_pathway
c=Get_CTD_pathway(train_fold,"sirolimus")
drug_file_dir,drug_dic_name=c.Get_trainFileDirList()
drugname,dic_pathway=c.Get_pathway_p(drug_file_dir,drug_dic_name)
CTDdic_dir=train_fold+'/CTDdic'
write_dic_to_file(CTDdic_dir,dic_pathway)

