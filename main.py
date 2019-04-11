# -*- coding: utf-8 -*-
'''
	Segmentacao de cromossomos.
	author: Rodrigo Junior santos
	Email:  rodrjuniorsantos@gmail.com
	Sistemas de Informacao
	Universidade Federal de Vicosa 
	Campus Rio Paranaiba

	
	Fonte: https://gist.github.com/kantale/e390cf7a47c4afdff9e4
		   http://pastebin.com/raw/6nBX6sdE
		   https://www.ncbi.nlm.nih.gov/genome/542?fbclid=IwAR0rVGQST-xj-ndCRZDxy1epKkdgskbt47Co1-t6malN949LhQXtBsUSUY8
'''

from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import os
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib      import pyplot       as plt


#from karioplot import mount, mountEdt
#--------------------------------------------------------------------------------------------

def karyoplot(karyo_filename, metadata={}):
	'''
	To create a karyo_filename go to: http://genome.ucsc.edu/cgi-bin/hgTables 
	group: Mapping and Sequencing
	track: Chromosome Band 
	An example of an output (hg19, Human) is here: http://pastebin.com/6nBX6sdE 
	The script will plot dots next to loci defined in metadata as:
	metadata = {
		'1' : [2300000, 125000000, 249250621],
	}
	'''


	karyo_dict={}
	with open(karyo_filename) as karyo_f:
		lines = [x.replace(os.linesep, '').split() for x in karyo_f.readlines()]

		for chromosome in [str(x) for x in range(1,25)]:
			karyo_dict[chromosome] = [[y[0], int(y[1]), int(y[2]), y[3], y[4]] for y in [x for x in lines if x[0] == 'chr' + chromosome]]

	fig, ax = plt.subplots(figsize=(15,15))

	DIM = 1.0

	#ax.set_xlim([0.0, DIM * (1.3)])
	#ax.set_ylim([0.0, DIM])

	def get_chromosome_length(chromosome):
		chromosome_start = float(min([x[1] for x in karyo_dict[chromosome]]))
		chromosome_end = float(max(x[2] for x in karyo_dict[chromosome]))
		chromosome_length = chromosome_end - chromosome_start

		return chromosome_length

	def plot_chromosome(chromosome, order):

		chromosome_length = get_chromosome_length(chromosome)
		chromosome_length_1 = get_chromosome_length('1')

		x_start = order * DIM * 0.1 
		x_end = x_start + (DIM * 0.04)
		y_start = DIM * 0.8 * (chromosome_length/chromosome_length_1)
		y_end = DIM * 0.1


		# We use the same colors as: http://circos.ca/tutorials/lessons/2d_tracks/connectors/configuration 
		colors = {
			#'gpos100' : (0/255.0,0/255.0,0/255.0),
			#'gpos'    : (0/255.0,0/255.0,0/255.0),
			#'gpos75'  : (130/255.0,130/255.0,130/255.0),
			#'gpos66'  : (160/255.0,160/255.0,160/255.0),
			#'gpos50'  : (200/255.0,200/255.0,200/255.0),
			#'gpos33'  : (210/255.0,210/255.0,210/255.0),
			'cinza'  : (200/255.0,200/255.0,200/255.0),
			#'gvar'    : (220/255.0,220/255.0,220/255.0),
			#'gneg'    : (255/255.0,255/255.0,255/255.0),
			'vermelho'    : (217/255.0,47/255.0,39/255.0),
			#'stalk'   : (100/255.0,127/255.0,164/255.0),
		}
		

		for index, piece in enumerate(karyo_dict[chromosome]):

			current_height = piece[2] - piece[1]
			current_height_sc = ((y_end - y_start) / chromosome_length) * current_height
			if index == 0:
				y_previous = y_start

			y_next = y_previous + current_height_sc

			color = colors[piece[4]]

			#plot the caryotypes
			r = Rectangle((x_start, y_previous), x_end-x_start, current_height_sc, color = color)
			ax.add_patch(r)

			y_previous = y_next

		#Plot semicircles at the beginning and end of the chromosomes
		center_x = x_start + (x_end-x_start)/2.0
		radius = (x_end-x_start)/2.0
		theta1 = 0.0
		theta2 = 180.0
		w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
		w2 = Wedge((center_x, y_end), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
		ax.add_patch(w1)
		ax.add_patch(w2)
		ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
		ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')

		#Plot metadata
		if chromosome in metadata:
			for md in metadata[chromosome]:
				ax.plot([x_end + (DIM*0.015)], [y_start + (y_end-y_start) * (md/chromosome_length)], '.', color='black')

		ax.text(center_x, y_end - (DIM * 0.07), chromosome)

	
	plot_chromosome('1', 1)
	plot_chromosome('2', 2)
	plot_chromosome('3', 3)
	plot_chromosome('4', 4)
	plot_chromosome('5', 5)
	plot_chromosome('6', 6)
	plot_chromosome('7', 7)
	plot_chromosome('8', 8)
	plot_chromosome('9', 9)
	plot_chromosome('10', 10)
	plot_chromosome('11', 11)
	plot_chromosome('12', 12)
	
	plot_chromosome('13', 13)
	plot_chromosome('14', 14)
	plot_chromosome('15', 15)
	plot_chromosome('16', 16)
	plot_chromosome('17', 17)
	plot_chromosome('18', 18)
	plot_chromosome('19', 19)
	plot_chromosome('20', 20)
	plot_chromosome('21', 21)
	plot_chromosome('22', 22)
	plot_chromosome('23', 23)
	plot_chromosome('24', 24)
		
	plt.axis('off')
	plt.show()

#main
def mount(file):

		fn = file
		karyoplot(fn)

def mountEdt():
	fn = "data//metadados.txt"

	directoryPath= dir="data"
	if not os.path.exists(directoryPath):
		os.mkdir( directoryPath ) 

	if not os.path.exists(fn):
		messagebox.showerror(title ="Message", message="Error reading text")

	karyoplot(fn)

#--------------------------------------------------------------------------------------------

#Fonts
def fontHelvetica():
	global text
	text.config(font="Helvetica")

def fontArial():
	global text
	text.config(font="Arial")

def fontTNR():
	global text
	text.config(font=("Times New Roman",))

def fontGothic():
	global text
	text.config(font="Gothic")

#Hotkeys
def newFileHK(self):
	global filename
	filename = "Untitled"
	text.delete(0.0, END)

def saveFileHK(self):
	global filename
	if filename == None:
		saveAs(self)
	else:
		t = text.get(0.0, END)
		f = open(filename, "w")
		f.write(t)
		f.close()

def saveAsHK(self):
	f = filedialog.asksaveasfile(mode='w', defaultextension='.txt')
	t = text.get(0.0, END)
	try:
		f.write(t.rstrip())
	except:
		messagebox.showerror(title ="UPS", message="Could not save file")

def openFileHK(self):
	f = filedialog.askopenfile(mode='r')
	t = f.read()
	text.delete(0.0, END)
	text.insert(0.0, t)

#Loops commands
filename = None

def newFile():
	global filename
	filename = "Untitled"
	text.delete(0.0, END)

def saveFile():
	global filename
	if filename == None:
		saveAs(self)
	else:
		t = text.get(0.0, END)
		f = open(filename, "w")
		f.write(t)
		f.close()

def saveAs():
	f = filedialog.asksaveasfile(mode='w', defaultextension='.txt')
	t = text.get(0.0, END)
	try:
		f.write(t.rstrip())
	except:
		messagebox.showerror(title ="UPS", message="Could not save file")

def openFile():
	f = filedialog.askopenfile(mode='r')
	t = f.read()
	text.delete(0.0, END)
	text.insert(0.0, t)

#-----------------------------  Montar  ----------------------------------------
def mountEditor():

	arquivo = open("data/metadados.txt", 'w+')
	t = text.get(0.0, END)

	if text.compare("end-1c", "==", "1.0"):
		messagebox.showerror("Erro", "Empty text area, enter metadata, please!")
	
	else:
		try:
			arquivo.write(t.rstrip())
		except:
			messagebox.showerror(title ="Message", message="Error reading text")

		arquivo.close()

		mountEdt()



def mountIdeograma():
	f = filedialog.askopenfilename(parent=root)
	mount(f)


def sobre():
	messagebox.showinfo("About", "Software for mounting an ideogram \n Made with coffee during the dawn \n author: Rodrigo Junior santos \n Email:  rodrjuniorsantos@gmail.com \n (Version 1.0 - Beta)")


#--------------------------------------------------------------------------------
#Main GUI
root = Tk()
root.title("Astyanax,")
#root.minsize(width=800, height=500)
#root.maxsize(width=800, height=500)
#root.setGeometry(50, 50, 1300, 650)
#root.geometry("850x450+250+250")

lado, cima = (root.winfo_screenwidth()), (root.winfo_screenheight())

root.geometry('%dx%d+0+0' % (lado,cima))

text = Text(root, width=500, height=500)
text.pack()


menubar = Menu(root)
filemenu = Menu(menubar)
filemenu.add_command(label="New", command=newFile, accelerator="Ctrl+N")
filemenu.add_command(label="Open", command=openFile, accelerator="Ctrl+O")
filemenu.add_command(label="Save", command=saveFile, accelerator="Ctrl+S")
filemenu.add_command(label="Save as...", command=saveAs, accelerator="Ctrl+A")
filemenu.add_separator()
filemenu.add_command(label="Quit", command=root.quit, accelerator="Ctrl+Q")
menubar.add_cascade(label="File", menu=filemenu)

#####Hot-Keys#######################
root.bind("<Control-n>", newFileHK)
root.bind("<Control-o>", openFileHK)
root.bind("<Control-s>", saveFileHK)
root.bind("<Control-a>", saveAsHK)
root.bind("<Control-q>", quit)
formatOptions = Menu(menubar)
formatOptions.add_command(label="Arial", command=fontArial)
formatOptions.add_command(label="Helvetica", command=fontHelvetica)
formatOptions.add_command(label="Times New Roman", command=fontTNR)
formatOptions.add_command(label="Gothic", command=fontGothic)
menubar.add_cascade(label="Font", menu=formatOptions)

formatOptions = Menu(menubar)
formatOptions.add_command(label="Mount file editor text", command=mountEditor)
formatOptions.add_command(label="Mount by external file", command=mountIdeograma)
menubar.add_cascade(label="Build Ideogram", menu=formatOptions)

filemenu = Menu(menubar)
menubar.add_cascade(label="About", command=sobre)

root.config(menu=menubar)
root.mainloop()