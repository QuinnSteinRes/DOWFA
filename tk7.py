import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTk, 
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from  matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
style.use('ggplot')


import urllib
import json
import sys
import pandas as pd
import numpy as np
import chaospy as cp

import tkinter as tk
from tkinter import ttk



LARGE_FONT= ("Verdana", 12)

f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)

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


class vaco2eapp(tk.Tk): #Baseline code for adding frams and pages

    def __init__(self, *args, **kwargs):  #Always run components when class is called apon 
        tk.Tk.__init__(self, *args, **kwargs)

        #tk.Tk.iconbitmap(self,default='HWLG2.ico') #Unsure of this but just as a point of  concept 
        tk.Tk.wm_title(self, "VACO")


        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        
        self.frames = {} 

        for F in (StartPage, Stock, WVA, matlib):

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
        label = tk.Label(self, text="Volatile Anaestheticn Equivalent CO2 Modeler: Successfully Loaded", font=LARGE_FONT)
        label.pack(pady=10,padx=10) # Should add info about paper, diclaimers, Updates etc 
    
        button1 = ttk.Button(self, text="Visit Stock Model",
                            command=lambda: controller.show_frame(Stock) )
        button1.pack()

        button2 = ttk.Button(self, text="Visit WVA Model",
                            command=lambda: controller.show_frame(WVA) )
        button2.pack()
        
        button3 = ttk.Button(self, text="Visit Graph Page",
                            command=lambda: controller.show_frame(matlib) )
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
        label = tk.Label(self, text="Waste Volatile Anaesthetic Emissions Model", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()
        
class matlib(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Graph Page!", font=LARGE_FONT)
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
ani = animation.FuncAnimation(f,animate, interval=1000) # interval in ms 
app.mainloop()
