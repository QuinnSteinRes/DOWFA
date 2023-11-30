import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTk, 
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.pyplot 
from  matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
style.use('ggplot')
from matplotlib import pyplot as plt

import yaml
with open('gloval.yml', 'r') as file:
    prime_service = yaml.safe_load(file)

import urllib
import json
import sys
import pandas as pd
import numpy as np
import chaospy as cp

import tkinter as tk
from tkinter import ttk

LARGE_FONT= ("Verdana", 12)
NORM_FONT = ("Helvetica", 10)
SMALL_FONT = ("Helvetica", 8)

f = plt.figure()
a = f.add_subplot(111)

def changeGWP(toWhat,pn): #?????
    global GWP
    global Counter
    global programName

    GWP = toWhat
    programName = pn
    Counter = 9000

def changeDestruction(toWhat,pn):
    global dmeth
    global Counter
    global programName

    GWP = toWhat
    programName = pn
    Counter = 9000



def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()

def animate(i):
    pullData = open('sampleText.txt','r').read()
    dataArray = pullData.split('\n')
    xar=[]
    yar=[]
    for eachLine in dataArray:
        if len(eachLine)>1:
            x,y = eachLine.split(',')
            xar.append(int(x))
            yar.append(int(y))
    a.clear()
    a.plot(xar,yar)
    
    a.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3,
             ncol=2, borderaxespad=0)               #Need to add legend for high and low estimation

    title = ("Modeled CO2e Emissions")
    a.set_title(title)

class vaco2eapp(tk.Tk): #Baseline code for adding frams and pages

    def __init__(self, *args, **kwargs):  #Always run components when class is called apon 
        tk.Tk.__init__(self, *args, **kwargs)

        #tk.Tk.iconbitmap(self,default='HWLG2.ico') #Unsure of this but just as a point of  concept 
        tk.Tk.wm_title(self, "VACO")


        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        
        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Save settings", command=lambda: popupmsg('Not supported just yet!'))
        filemenu.add_command(label="Edit Parameters", command=lambda: popupmsg('Not supported just yet!'))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)
        menubar.add_cascade(label="File", menu=filemenu)

        GWPChoice = tk.Menu(menubar, tearoff=1)
        GWPChoice.add_command(label="GWP20 Year",
                                   command=lambda: changeGWP("GWP20","gwp20"))
        GWPChoice.add_command(label="GWP100 Year",
                                   command=lambda: changeGWP("GWP100","gwp100"))
        GWPChoice.add_command(label="GWP500 Year",
                                   command=lambda: changeGWP("GWP500","gwp500"))
        
        menubar.add_cascade(label="GWP Period", menu=GWPChoice)    
        
        
        dmethChoice = tk.Menu(menubar, tearoff=1)
        dmethChoice.add_command(label="Release",
                                   command=lambda: changeDestruction("Release","release"))
        dmethChoice.add_command(label="MSWI",
                                   command=lambda: changeDestruction("MSWI","msw"))
        dmethChoice.add_command(label="PWI",
                                   command=lambda: changeDestruction("PWI","pwi"))
        dmethChoice.add_command(label="Plasma",
                                   command=lambda: changeDestruction("Plasma","plasma"))           
                        
        menubar.add_cascade(label="Destruction Method", menu=dmethChoice) 


        tk.Tk.config(self, menu=menubar)
        
        
        self.frames = {} 

        for F in (StartPage, WVA, Stock):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise() 
        
        

class StartPage(tk.Frame):

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent) #Parent is the vaco2e class above
        label = tk.Label(self, text="""Volatile Anaestheticn Equivalent CO2 Modeler:
        SuccessfullyLoaded""", font=LARGE_FONT)
        label.pack(pady=10,padx=10) # Should add info about paper, diclaimers, Updates etc 
    
        button1 = ttk.Button(self, text="Legacy Stock",
                            command=lambda: controller.show_frame(Stock) )
        button1.pack()

        button2 = ttk.Button(self, text="WVA Releace",
                            command=lambda: controller.show_frame(WVA) )
        button2.pack()
        
        button3 = ttk.Button(self, text="Exit",
                            command= quit )
        button3.pack()
        


class Stock(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Legacy Stock Emissions Model", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()
        
class WVA(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="WVA Emissions Model", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()
        
        a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])

        canvas = FigureCanvasTkAgg(f, self) #TkAgg working 
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self) #Only Tk for Nav working 
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)  
        
        

        



app = vaco2eapp()
app.geometry("1280x720")
ani = animation.FuncAnimation(f,animate, interval=10000) # interval in ms 
app.mainloop()
