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
from matplotlib.collections import PatchCollection, BrokenBarHCollection
from matplotlib      import pyplot       as plt


colors = {
	        'gneg': (1., 1., 1.),
	        'gpos25': (.6, .6, .6),
	        'gpos50': (.4, .4, .4),
	        'gpos75': (.2, .2, .2),
	        'gpos100': (0., 0., 0.),
	        'acen': (.8, .4, .4),
	        'gvar': (.8, .8, .8),
	        'stalk': (.9, .9, .9),
	        'red'  : (217/255.0,47/255.0,39/255.0),
	        'darkred': '#8b0000',
	        'gray' : (200/255.0,200/255.0,200/255.0),
	        'lightgray': '#d3d3d3',
	        'white': '#ffffff',
	        'whitesmoke': '#f5f5f5',
	        'black':'#000000',
	        'silver': '#c0c0c0',
	        'light': '#d3d3d3',
	        'blue': '#0000ff',
	        'slateblue': '#6a5acd',
	        'darkblue': '#00008b',
	        'royalblue': '#4169e1',
	        'mediumblue': '#0000cd',
	        'deepskyblue': '#00bfff',
	        'skyblue': '#87ceeb',
	        'green': '#008000',
	        'lime': '#00ff00',
	        'darkgreen': '#006400',
	        'orange': '#ffa500',
	        'wheat': '#f5deb3',
	        'gold': '#ffd700',
	        'purple': '#800080',
	        'mediumpurple': '#9370db',
	        'darkorchid': '#9932cc',
	        'pink': '#ffc0cb',
	        'lightpink': '#ffb6c1',
	        'teal': '#008080',
	    }


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#						IDEOGRAM
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

def karyoplot(karyo_filename, metadata={}):

	karyo_dict={}
	karyo_={}
	with open(karyo_filename) as karyo_f:
		lines = [x.replace(os.linesep, '').split() for x in karyo_f.readlines()]

		chroms=[]
		for line in lines:
			chrom, start, stop, name, stain = line
			chroms.append(chrom)
		chroms.pop(0)
		unique = set(chroms) 

		for chromosome in [str(x) for x in range(1, len(unique)+1)]:
			karyo_dict[chromosome] = [[y[0], float(y[1]), float(y[2]), y[3], y[4]] for y in [x for x in lines if x[0] == 'chr' + chromosome]]

	fig, ax = plt.subplots(figsize=(15,15))

	DIM = 1.0

	def get_centromere(chromosome):
		chromosome_end = float(max([x[2] for x in karyo_dict[chromosome]]))
		for x in karyo_dict[chromosome]:
			if x[3] == 'centromere':
				start_= float(x[1])
				end_  = float(x[2])
		centromere_length_start= chromosome_end - start_
		centromere_length_end= chromosome_end - end_

		return centromere_length_start, centromere_length_end
	
	def get_chromosome_length(chromosome):
		chromosome_start = float(min([x[1] for x in karyo_dict[chromosome]]))
		chromosome_end = float(max(x[2] for x in karyo_dict[chromosome]))
		chromosome_length = chromosome_end - chromosome_start

		return chromosome_length

	def plot_chromosome(chromosome, order):

		chromosome_length = get_chromosome_length(chromosome)
		chromosome_length_1 = get_chromosome_length('1')

		centromere_length, centromere_length_= get_centromere(chromosome)
		centromere_length_1, centromere_length_1_= get_centromere('1')

		x_start = order * DIM * 0.1 
		x_end = x_start + (DIM * 0.04)
		y_start = DIM * 0.8 * (chromosome_length/chromosome_length_1)
		y_end = DIM * 0.1

		y_centromere_start=DIM * 0.62 * (centromere_length/centromere_length_1)
		y_centromere_end=DIM * 0.62 * (centromere_length_/centromere_length_1_)

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
		theta1 = 0
		theta2 = 180.0

		w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
		w2 = Wedge((center_x, y_end), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
		ax.add_patch(w1)
		ax.add_patch(w2)

		w1_ = Wedge((center_x, y_centromere_start), radius, theta2, theta1, width=0.01, facecolor='red', edgecolor='black')
		w2_ = Wedge((center_x, y_centromere_start), radius, theta1, theta2, width=0.01, facecolor='red', edgecolor='black')
		ax.add_patch(w1_)
		ax.add_patch(w2_)

		ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
		ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')

		#Plot metadata
		if chromosome in metadata:
			for md in metadata[chromosome]:
				ax.plot([x_end + (DIM*0.015)], [y_start + (y_end-y_start) * (md/chromosome_length)], '.', color='black')

		ax.text(center_x, y_end - (DIM * 0.07), chromosome)


	for i in range(1, len(unique)+1):
		plot_chromosome(str(i), i)

	plt.title('Ideogram')
	plt.axis('off')
	plt.show(block = False)

#--------------------------
def mount(file):
	try:
		karyoplot(file)
	except:
		messagebox.showerror(title ="Error reading file", message="The file should be in the correct format for cytogenomic, read the tutorial.")
	


def mountEdt():
	fn = "data//function1.txt"

	directoryPath= "data"
	if not os.path.exists(directoryPath):
		os.mkdir( directoryPath ) 

	if not os.path.exists(fn):
		messagebox.showerror(title ="Message", message="Error reading text")

	try:
		karyoplot(fn)
	except:
		messagebox.showerror(title ="Error reading file", message="The file should be in the correct format for cytogenomic, read the tutorial.")
	

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#						CHROMOSOMES DETAILS
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------


def details(file):
	
	height = 0.9
	spacing = 0.9

	def ideograms(fn):
	    last_chrom = None
	    last_specie= None
	    last_population= None
	    fin = open(fn)
	    fin.readline()
	    xranges, color = [], []
	    ymin = 0

	    for line in fin:
	        chrom, start, stop, specie, population, stain = line.strip().split()
	        start = int(start)
	        stop = int(stop)
	        width = stop - start
	        if chrom == last_chrom or (last_chrom is None):
	            xranges.append((start, width))
	            color.append(colors[stain])
	            last_chrom = chrom
	            last_specie= specie
	            last_population= population
	            continue

	        ymin += height + spacing
	        yrange = (ymin, height)
	        yield xranges, yrange, color, last_chrom, last_specie, last_population
	        xranges, color = [], []
	        xranges.append((start, width))
	        color.append(colors[stain])
	        last_chrom = chrom
	        last_specie= specie
	        last_population= population

	    # last one
	    ymin += height + spacing
	    yrange = (ymin, height)
	    yield xranges, yrange, color, last_chrom, last_specie, last_population

	fig = plt.figure(figsize=(15,15))
	ax = fig.add_subplot(111)
	d = {}
	yticks = []
	yticklabels = []

	for xranges, yrange, color, chms, species, populations in ideograms(file):
	    coll = BrokenBarHCollection(xranges, yrange, facecolors=color)
	    ax.add_collection(coll)
	    center = yrange[0] + yrange[1]/2.
	    yticks.append(center)
	    label= '%s %s'%(species, populations)
	    yticklabels.append(label)
	    d[chms] = xranges
	    values=[]
	    bp=[]
	    for inter in xranges:
	        values.append(inter[0])
	        bp.append(inter[1])
	    values.append(bp[-1]+values[-1])

	    for i in range(0, len(values)-1):
	        xlabel= '%d bp'%(bp[i])
	        ax.annotate(xlabel, xy=(values[i+1]/2 + values[i]/2, yrange[0] + yrange[1]), xytext=(values[i+1]/2 + values[i]/2, yrange[0] + yrange[1]))

	ax.set_title('Chromosomes')
	ax.axis('tight')
	ax.set_yticks(yticks)
	ax.set_yticklabels(yticklabels)
	ax.set_xticks([])
	plt.show(block = False)


#--------------------------
def mount_details(file):
	try:
		details(file)
	except:
		messagebox.showerror(title ="Error reading file", message="The file should be in the correct format for cytogenetics, read the tutorial.")
	


def mountEdt_details():
	fn = "data//function3.txt"

	directoryPath= "data"
	if not os.path.exists(directoryPath):
		os.mkdir( directoryPath ) 

	if not os.path.exists(fn):
		messagebox.showerror(title ="Message", message="Error reading text")

	try:
		details(fn)
	except:
		messagebox.showerror(title ="Error reading file", message="The file should be in the correct format for cytogenetics, read the tutorial.")
	



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#						MENUS FUNCTIONS
#--------------------------------------------------------------------------------------------
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


def quit():
	sys.exit()




def mountEditor():

	arquivo = open("data/function1.txt", 'w+')
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


def mountFile():
	f = filedialog.askopenfilename(parent=root)
	mount(f)


def mountEdtDetails():

	arquivo = open("data/function3.txt", 'w+')
	t = text.get(0.0, END)

	if text.compare("end-1c", "==", "1.0"):
		messagebox.showerror("Erro", "Empty text area, enter metadata, please!")
	
	else:
		try:
			arquivo.write(t.rstrip())
		except:
			messagebox.showerror(title ="Message", message="Error reading text")

		arquivo.close()

		mountEdt_details()


def mountFileDetails():
	f = filedialog.askopenfilename(parent=root)
	mount_details(f)


def about():
	messagebox.showinfo("About", "Software for mounting an ideogram \n Made with coffee during the dawn \n author: Rodrigo Junior santos \n Email:  rodrjuniorsantos@gmail.com \n (Version 1.0 - Beta)")




#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#						MAIN
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

if __name__ == '__main__':
	#Main GUI
	root = Tk()
	root.title("Astyanax")
	root.iconbitmap("icon.ico")

	side, top = (root.winfo_screenwidth()), (root.winfo_screenheight())

	root.geometry('%dx%d+0+0' % (side, top))

	text = Text(root, width=500, height=500)
	text.pack()


	menubar = Menu(root)
	filemenu = Menu(menubar)
	filemenu.add_command(label="New", command=newFile, accelerator="Ctrl+N")
	filemenu.add_command(label="Open", command=openFile, accelerator="Ctrl+O")
	filemenu.add_command(label="Save", command=saveFile, accelerator="Ctrl+S")
	filemenu.add_command(label="Save as...", command=saveAs, accelerator="Ctrl+A")
	filemenu.add_separator()
	filemenu.add_command(label="Quit", command=quit, accelerator="Ctrl+Q")
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
	formatOptions.add_command(label="Mount ideogram by text editor", command=mountEditor)
	formatOptions.add_command(label="Mount by external file", command=mountFile)
	menubar.add_cascade(label="Cytogenomics", menu=formatOptions)

	formatOptions = Menu(menubar)
	formatOptions.add_command(label="Mount by text editor", command=mountEditor)
	formatOptions.add_command(label="Mount by external file", command=mountFile)
	menubar.add_cascade(label="Cytogenetics", menu=formatOptions)

	formatOptions = Menu(menubar)
	formatOptions.add_command(label="Mount by text editor", command=mountEdtDetails)
	formatOptions.add_command(label="Mount by external file", command=mountFileDetails)
	menubar.add_cascade(label="Details", menu=formatOptions)

	filemenu = Menu(menubar)
	menubar.add_cascade(label="About", command=about)

	root.config(menu=menubar)
	root.mainloop()