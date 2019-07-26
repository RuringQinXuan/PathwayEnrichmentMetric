'''
画AUC有两个思路
1，根据AUC文件画图
2，根据evaluation文件画图
这个文件先写根据evaluation文件画图
问题在于AUC文件得到的各个indicator的列表长度不一样
'''
import numpy as np
import matplotlib as mpl
mpl.use('Agg')	
import matplotlib.pyplot as plt	
plt.figure(num=1, figsize=(12, 8))
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_drugs(drugs_file_dir):
	drugs=[]
	with open(drugs_file_dir,'r') as drugs_file:
		for line in drugs_file:
			line= line.strip()
			drugs.append(line)
	return drugs
	
def get_dic_from_file(file_dir):
	with open(file_dir,'r') as file:
		dic=eval(file.read())
	return dic
	

def get_standard_pathwaylist(drugs,methods,indicators,KEGG_pathwaylist_dir):
	baseline_pathwaylist=get_dic_from_file(KEGG_pathwaylist_dir)
	standardPathwayList=set(baseline_pathwaylist)
	for drug in drugs:
		for method in methods:
			for indicator in indicators:
				filedir='./evaluation/rank_data/'+drug+'_'+method+'_'+indicator
				dic=get_dic_from_file(filedir)
				standardPathwayList=set(dic.keys())&standardPathwayList
	with open('./data/pathway/standardKEGG_pathwaylist','w') as file:
		file.write(str(standardPathwayList))
	return standardPathwayList

def get_evaluationResult_standard(drugs,methods,indicators,standardPathwayList):
	for drug in drugs:
		for method in methods:
			for indicator in indicators:
				filedir='./evaluation/rank_data/'+drug+'_'+method+'_'+indicator
				dic=get_dic_from_file(filedir)
				standarddic={}
				for pathway in dic.keys():
					if pathway in standardPathwayList:
						standarddic[pathway]=dic[pathway]
				with open('./evaluation/standard_rank_result/'+drug+'_'+method+'_'+indicator,'w') as file:
					file.write(str(standarddic))
					
# with open('./evaluation/picture/colorname_dic','w') as file:
	# file.write(str(cnames))
					
					
def get_picture_data(drugs,methods,indicators,standardPathwayList,CTD_pathway_list_dic):
	for drug in drugs:
		CTD_pathway_list=CTD_pathway_list_dic[drug].keys()
		for method in methods:
			for indicator in indicators:
				accumulate=0
				accumulateList=[]
				filedir='./evaluation/standard_rank_result/'+drug+'_'+method+'_'+indicator
				dic=get_dic_from_file(filedir)
				resultdir='./evaluation/picture/data/'+drug+'_'+method+'_'+indicator
				for pathway in dic.keys():
					if pathway in CTD_pathway_list:
						accumulate+=1
					accumulateList.append(accumulate)
				temp = np.array(accumulateList) 
				np.savetxt(resultdir,temp,fmt='%0.0f')
				
def draw_colorbar_picture(drug,CTD_pathway_list_dic):
	CTD_pathway_list=CTD_pathway_list_dic[drug]
	for method in methods:
		color_CTD_data=[]
		for data in CTD_pathway_list.values():
			color_CTD_data.append(data)
		for i in range(len(CTD_pathway_list.values()),302):
			color_CTD_data.append(1)
		color_CTD_data=np.array(color_CTD_data)
		a=color_CTD_data
		for indicator in indicators:
			filedir='./evaluation/standard_rank_result/'+drug+'_'+method+'_'+indicator
			dic=get_dic_from_file(filedir)
			b=[]
			for pathway in dic.keys():
				if pathway in CTD_pathway_list.keys():
					b.append(CTD_pathway_list[pathway])
				else:
					b.append(1)
			b=np.array(b)
			a=np.concatenate(a,b,axis=1)
		fig, ax = plt.subplots()
		im = ax.imshow(harvest)
		x_labels=indicators+['CTD']
		y_labels=range(1,302)
		ax.set_xticks(np.arange(len(x_labels)))
		ax.set_yticks(np.arange(len(y_labels)))
		ax.set_xticklabels(x_labels)
		ax.set_yticklabels(y_labels)
		plt.setp(ax.get_xticklabels(), rotation=45, ha="right",rotation_mode="anchor")
		ax.set_title("comparison the evaluation result pathway list with CTD known pathway list about "+drug)
		fig.tight_layout()
		plt.savefig('./evaluation/picture/picture/barPic/'+drug+'_'+method+'.png', dpi=600)


				
def draw_picture(drug):
	plt.title(drug+'method and indicator comparison')
	#画图属性
	color_names=['red','darkorange','black','darkgreen','navy','purple','deepskyblue','peachpuff']
	subfig=0
	layer=1
	font_legend = {
	'weight' : 'normal', 
	'size'  :8 , 
	} 
	font_axis = {
	'weight' : 'normal', 
	'size'  : 15, 
	} 
	font_title={
	'weight' : 'normal', 
	'size'  : 10, 
	}
	fig=plt.figure()
	fig.suptitle("Areas comparison of indicators in "+drug, fontsize=16)
	figsize = 8, 9
	plt.subplots(figsize=figsize)                # 设定整张图片大小
	#图像数据
	baseline=np.array(range(301))
	a=baseline
	# for drug in drugs:
	for method in methods:
		areas=[]
		subfig+=1
		ax = plt.subplot(len(methods),2,subfig)
		ax.yaxis.set_major_locator(MultipleLocator(20))
		color_i=-1
		for indicator in indicators:
			#获得数据
			pic_data_dir='./evaluation/picture/data/'+drug+'_'+method+'_'+indicator
			b=np.loadtxt(pic_data_dir,dtype=np.float32)
			area=round(sum(b)/(b[-1]*301),3)
			areas.append(area)
			#画图
			color_i+=1
			colorname=color_names[color_i]
			plt.plot(baseline,b,color=colorname,label=color_i,linewidth=0.8)
			hl=plt.legend(loc='upper right', prop=font_legend, frameon=False)
			plt.ylim(0, b[-1]+10)
			plt.xlim(0, len(b)+10)
			# if method == methods[0]:
			ax.set_title("Areas comparison of indicators in "+method.upper(),font_title)
			if method ==methods[1]:
				plt.ylabel("the known pathway number covered", font_axis)
			if method==methods[-1] and indicator==indicators[-1]:
				plt.xlabel("the $x_{th}$ pathway in ranking result", font_axis)
				ax.xaxis.set_major_locator(MultipleLocator(50))
				continue
			plt.xticks([])# 去掉x坐标轴刻度
		subfig+=1
		ax = plt.subplot(len(methods),2,subfig) 
		plt.bar(range(8),areas,color=color_names)
		plt.ylim(0.5, 1)
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		ax.set_title("Areas comparison of indicators in "+method.upper(),font_title)		
		#使用text显示数值  
		for a,b in zip(range(8),areas):  
			if b==max(areas):
				plt.text(a,b+0.05, b, ha='center', va= 'bottom',fontsize=7,color='r') 
			else:
				plt.text(a,b+0.05, b, ha='center', va= 'bottom',fontsize=7) 
		if method ==methods[1]:
			plt.ylabel("The areas", font_axis)
		if method==methods[-1] and indicator==indicators[-1]:
			plt.xlabel("The indicators", font_axis)
			plt.xticks(range(8))
			continue
		plt.xticks([])
	plt.savefig('./evaluation/picture/picture/curvePic/'+drug+".png", dpi=600)		
	

# 校准pathwaylist数据
KEGG_pathwaylist_dir='./data/pathway/KEGG/dic_pathways_genenumber'
standardPathwayList=get_standard_pathwaylist(drugs,methods,indicators,KEGG_pathwaylist_dir)
get_evaluationResult_standard(drugs,methods,indicators,standardPathwayList)
					
#循环数据
indicators=['P','IPF_gene','IPF_node','IPF_short','IPF_short_gene','IPF_short_node','IPF_gene_node','IPF_gene_node_short']
methods=['abstract','sentence','dependency','tees']
drugs_file_dir='./data/drug/drug_name'
drugs=get_drugs(drugs_file_dir)
CTD_pathway_LIST_DIR='./data/pathway/CTD/CTDdic'
CTD_pathway_list_dic=get_dic_from_file(CTD_pathway_LIST_DIR)
# get_picture_data(drugs,methods,indicators,standardPathwayList,CTD_pathway_list_dic)
				
#画图数据
color_dir='./evaluation/picture/colorname_dic'
drug='sirolimus'
method='tees'
indicator='IPF_gene'
CTD_pathway_list=CTD_pathway_list_dic[drug]	
#做一个大字典，把indicator作为key，indicator对应pathwaylist放在value里
StandardPathway_dir='./data/pathway/standard_pathwaylist'
KEGG_pathwaylist_dir='./data/pathway/KEGG/dic_pathways_genenumber'
standardPathwayList=get_dic_from_file(StandardPathway_dir)
for drug in drugs:
	draw_picture(drug)