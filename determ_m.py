import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTk, 
    NavigationToolbar2Tk
)
from matplotlib.figure import Figure
import matplotlib.pyplot 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
style.use('ggplot')
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import numpy as np


import urllib
import json
import sys
import pandas as pd
import numpy as np
import chaospy as cp

import tkinter as tk
from tkinter import ttk

from parameters_m import *

n_samples = 10000

# Deterministic calculations class
# from deterministicCalcs import DeterministicCalcs
# dCalcs = DeterministicCalcs(params)    

LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Helvetica", 10)
SMALL_FONT = ("Helvetica", 8)


def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command=popup.destroy)
    B1.pack()
    popup.mainloop()


class vaco2eapp(tk.Tk):  # Baseline code for adding frames and pages

    def __init__(self, *args, **kwargs):  # Always run components when class is called upon 
        tk.Tk.__init__(self, *args, **kwargs)

        # tk.Tk.iconbitmap(self,default='HWLG2.ico') # Unsure of this but just as a point of  concept 
        tk.Tk.wm_title(self, "VACO")

        container = ttk.Frame(self, padding="50 50 50 50")
        container.grid(column=0, row=0)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Save settings", command=lambda: popupmsg('Not supported just yet!'))
        filemenu.add_command(label="Edit Parameters", command=lambda: popupmsg('Not supported just yet!'))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)
        menubar.add_cascade(label="File", menu=filemenu)

        tk.Tk.config(self, menu=menubar)

        self.frames = {}

        for F in (StartPage, determ, WVA):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class StartPage(tk.Frame):

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)  # Parent is the vaco2e class above
        label = tk.Label(self, text="""Volatile Anaestheticn Equivalent CO2 Modeler: Successfully Loaded""", font=LARGE_FONT)
        label.grid(column=1, row=1)  # Should add info about paper, disclaimers, Updates etc 

        button1 = ttk.Button(self, text="Determ Stock", command=lambda: controller.show_frame(determ))
        button1.grid(column=1, row=4)
        
        button2 = ttk.Button(self, text="WVA", command=lambda: controller.show_frame(WVA))
        button2.grid(column=1, row=5)

        button3 = ttk.Button(self, text="Exit", command=quit)
        button3.grid(column=1, row=6)

class determ(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Legacy Stock Emissions Model", font=LARGE_FONT)
        label.grid(column=1, row=1, sticky=tk.W, padx=5, pady=5)

        # Entry Button 
        self.entry1 = ttk.Entry(self, width=7)
        self.entry1.grid(column=2, row=5)
        
        # Bind <Return> key to the entry field to trigger the calculate method
        self.entry1.bind('<Return>', lambda event: self.calculate())

        # Selection Buttons
        self.variable_options = ["des", "sev", "iso"]
        self.gwp_options = ["GWP_20", "GWP_100", "GWP_500"]

        self.variable_var = tk.StringVar()
        self.variable_dropdown = ttk.Combobox(self, textvariable=self.variable_var, values=self.variable_options)
        self.variable_dropdown.grid(column=15, row=10)
        self.variable_dropdown.current(0)  # Set default value

        self.gwp_var = tk.StringVar()
        self.gwp_dropdown = ttk.Combobox(self, textvariable=self.gwp_var, values=self.gwp_options)
        self.gwp_dropdown.grid(column=15, row=11)
        self.gwp_dropdown.current(0)  # Set default value

        # Calc Button
        calc1 = ttk.Button(self, text="Calculate", command=self.calculate)
        calc1.grid(column=12, row=5)

        # Outputs
        self.output1 = ttk.Label(self, text="")
        self.output1.grid(column=2, row=7)

        self.output2 = ttk.Label(self, text="")
        self.output2.grid(column=2, row=8)

        self.output3 = ttk.Label(self, text="")
        self.output3.grid(column=2, row=9)

        self.output4 = ttk.Label(self, text="")
        self.output4.grid(column=2, row=10)
        
        self.output5 = ttk.Label(self, text="")
        self.output5.grid(column=3, row=7)

        self.output6 = ttk.Label(self, text="")
        self.output6.grid(column=3, row=8)

        self.output7 = ttk.Label(self, text="")
        self.output7.grid(column=3, row=9)

        self.output8 = ttk.Label(self, text="")
        self.output8.grid(column=3, row=10)

        # Table Labels
        
        label11 = ttk.Label(self, text="Cumulative (tonnesCO2 - eq)")
        label11.grid(column=2, row=6, sticky=tk.E, padx=5, pady=5)
        
        label12 = ttk.Label(self, text="kgCO2 - eq per kg of VA")
        label12.grid(column=3, row=6, sticky=tk.E, padx=5, pady=5)

        # LHS Labels
        label1 = ttk.Label(self, text="VA Stock Input")
        label1.grid(column=1, row=5, sticky=tk.W, padx=5, pady=5)
        label2 = ttk.Label(self, text="Release")
        label2.grid(column=1, row=7, sticky=tk.W, padx=5, pady=5)
        label3 = ttk.Label(self, text="MSWI")
        label3.grid(column=1, row=8, sticky=tk.W, padx=5, pady=5)
        label4 = ttk.Label(self, text="PWI")
        label4.grid(column=1, row=9, sticky=tk.W, padx=5, pady=5)
        label5 = ttk.Label(self, text="Plasma")
        label5.grid(column=1, row=10, sticky=tk.W, padx=5, pady=5)

        # RHS Labels
        label10 = ttk.Label(self, text="Litres")
        label10.grid(column=3, row=5, sticky=tk.E, padx=5, pady=5)
        

        # Buttons 
        button1 = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame(StartPage))
        button1.grid(column=15, row=14)

        button3 = ttk.Button(self, text="Exit", command=quit)
        button3.grid(column=15, row=15)

    def calculate(self):
        try:
            value = float(self.entry1.get())
            gwp_selection = self.gwp_var.get()  # Get the selected GWP value
            variable_selection = self.variable_var.get()  # Get the selected variable (des, sev, iso)
            
            # Define dictionaries and values for calculations
            gwp_categories = {
                "GWP_20": {"des": des_20, "sev": sev_20, "iso": iso_20},
                "GWP_100": {"des": des_100, "sev": sev_100, "iso": iso_100},
                "GWP_500": {"des": des_500, "sev": sev_500, "iso": iso_500}
            }
            
            rho_values = {
                "des": rho_d,
                "sev": rho_s,
                "iso": rho_i
            }
            
            gfm_values = {
                "des": des_gfm,
                "sev": sev_gfm,
                "iso": iso_gfm
            }
            
            carbon_economy_values = {
                "des": des_atecon_C,
                "sev": des_atecon_C,
                "iso": iso_atecon_C
            }
            
            cf4_gwp_values = {
                "GWP_20": cf4_20,
                "GWP_100": cf4_100,
                "GWP_500": cf4_500
            }

            if gwp_selection in gwp_categories:
                rho = rho_values.get(variable_selection)
                gfm = gfm_values.get(variable_selection)
                carbon_economy = carbon_economy_values.get(variable_selection)
                cf4_gwp = cf4_gwp_values.get(gwp_selection)

                if all(item is not None for item in (rho, gfm, carbon_economy, cf4_gwp)):
                    mol_val = (1000 / gfm)
                    
                    demunic = np.random.uniform(min_demunic, max_demunic, n_samples)
                    depharma = np.random.uniform(min_depharma, max_depharma, n_samples)
                    deplasma = np.random.uniform(min_deplasma, max_deplasma, n_samples)
                    
                    munic_other = np.random.uniform(min_munic_other, max_munic_other, n_samples)
                    pharam_other = np.random.uniform(min_pharam_other, max_pharam_other, n_samples)
                    plasma_other = np.random.uniform(min_plasma_other, max_plasma_other, n_samples)
                    
                    sur_100k = np.random.uniform(min_sur_100k, max_sur_100k, n_samples)
                    growth_rate = np.random.uniform(min_growth_rate, max_growth_rate, n_samples)
                    
                    van_factor = np.random.uniform(min_van_factor, max_van_factor, n_samples)
                    hgv_factor = np.random.uniform(min_hgv_factor, max_hgv_factor, n_samples)
                    ship_factor = np.random.uniform(min_ship_factor, max_ship_factor, n_samples)
                    
                    network_distance = np.random.uniform(min_network_distance, max_network_distance, n_samples)
                    mwi_dist = np.random.uniform(min_mwi_dist, max_mwi_dist, n_samples)
                    pwi_dist = np.random.uniform(min_pwi_dist, max_pwi_dist, n_samples)
                    
                                                    
                    
                    
                    convert_L_kg = (0.001 * value * rho)
                    
                    combust_value1 = ((((mol_val * (1 - demunic)) * carbon_economy * 44.1)) / 1000) + \
                                    (((((mol_val * (demunic)) * carbon_economy * 88.0043)) / 1000) * cf4_gwp)
                    
                    combust_value2 = ((((mol_val * (1 - depharma)) * carbon_economy * 44.1)) / 1000) + \
                                    (((((mol_val * (depharma)) * carbon_economy * 88.0043)) / 1000) * cf4_gwp)             
                    
                    combust_value3 = ((((mol_val * (1 - deplasma)) * carbon_economy * 44.1)) / 1000) + \
                                    (((((mol_val * (deplasma)) * carbon_economy * 88.0043)) / 1000) * cf4_gwp)
                    
                    
                    calculated_value1 = (convert_L_kg *gwp_categories[gwp_selection][variable_selection] ) / 1000
                    calculated_value2 = ((mwi_dist * (convert_L_kg/1000) * van_factor) + (convert_L_kg * combust_value1) + (convert_L_kg * munic_other)) / 1000
                    calculated_value3 = ((pwi_dist * (convert_L_kg/1000) * van_factor) + (convert_L_kg * combust_value2) + (convert_L_kg * pharam_other)) / 1000
                    calculated_value4 = ((network_distance * (convert_L_kg/1000) * van_factor) + (convert_L_kg * combust_value3) + (convert_L_kg * plasma_other)) / 1000
                    
                    calculated_value5 = (calculated_value1 *1000) / convert_L_kg
                    calculated_value6 = (calculated_value2 *1000) / convert_L_kg
                    calculated_value7 = (calculated_value3 *1000) / convert_L_kg
                    calculated_value8 = (calculated_value4 *1000) / convert_L_kg
                    
                    calculated_value9 =  (van_factor * network_distance)/1000
                    
                    mean_cv_1 = np.mean(calculated_value1)
                    mean_cv_2 = np.mean(calculated_value2)
                    mean_cv_3 = np.mean(calculated_value3)
                    mean_cv_4 = np.mean(calculated_value4)
                    mean_cv_5 = np.mean(calculated_value5)
                    mean_cv_6 = np.mean(calculated_value6)
                    mean_cv_7 = np.mean(calculated_value7)
                    mean_cv_8 = np.mean(calculated_value8)
                    mean_cv_9 = np.mean(calculated_value9)
                    
                    sd_cv_1 = np.std(calculated_value5)
                    sd_cv_2 = np.std(calculated_value6)
                    sd_cv_3 = np.std(calculated_value7)
                    sd_cv_4 = np.std(calculated_value8)
                    sd_cv_5 = np.std(calculated_value9)
                    
                    

                    # Display calculated values in respective output labels
                    self.output1.config(text=f'{mean_cv_1:.1f}')
                    self.output2.config(text=f'{mean_cv_2:.1f}')
                    self.output3.config(text=f'{mean_cv_3:.1f}')
                    self.output4.config(text=f'{mean_cv_4:.1f}')
                    
                    self.output5.config(text=f'{mean_cv_5:.1f}')
                    self.output6.config(text=f'{mean_cv_6:.1f}')
                    self.output7.config(text=f'{mean_cv_7:.1f}')
                    self.output8.config(text=f'{mean_cv_8:.1f}')
                    #print(calculated_value6)
                    print(sd_cv_1)
                    print(sd_cv_2)
                    print(sd_cv_3)
                    print(sd_cv_4)
                    print(sd_cv_5)
                    print(mean_cv_9)
                    
                    
                    #plt.hist(calculated_value3)
                    #plt.xlabel('R1')
                    #plt.ylabel('Frequency')
                    #plt.show()


        except ValueError:
            pass

class WVA(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="WVA Emissions Model", font=LARGE_FONT)
        label.grid(column=1, row=1, sticky=tk.W, padx=5, pady=5)
        
        separator = ttk.Separator(self, orient='horizontal')
        separator.grid(row=2, column=1, columnspan=20, sticky="ew", padx=5, pady=5)

                
        # Entry Button 
        self.entry1 = ttk.Entry(self, width=7) #period
        self.entry1.grid(column=2, row=5)
        self.entry2 = ttk.Entry(self, width=7) #init pop
        self.entry2.grid(column=2, row=6)
        self.entry3 = ttk.Entry(self, width=7) #des % surgery account
        self.entry3.grid(column=2, row=7)
        
        # Selection Buttons
        self.gwp_options = ["GWP_100"]

        self.gwp_var = tk.StringVar()
        self.gwp_dropdown = ttk.Combobox(self, textvariable=self.gwp_var, values=self.gwp_options)
        self.gwp_dropdown.grid(column=5, row=7, sticky=tk.W, padx=5, pady=5)
        self.gwp_dropdown.current(0)  # Set default value

        # Calc Button
        calc1 = ttk.Button(self, text="Calculate", command=self.calculate)
        calc1.grid(column=4, row=7)

        # Outputs
        self.output1 = ttk.Label(self, text="")
        self.output1.grid(column=2, row=9)

        self.output2 = ttk.Label(self, text="")
        self.output2.grid(column=2, row=10)

        self.output3 = ttk.Label(self, text="")
        self.output3.grid(column=2, row=11)

        self.output4 = ttk.Label(self, text="")
        self.output4.grid(column=2, row=12)
        
        self.output5 = ttk.Label(self, text="")
        self.output5.grid(column=3, row=13)

        self.output6 = ttk.Label(self, text="")
        self.output6.grid(column=3, row=14)

        # Table Labels
        
        label11 = ttk.Label(self, text="Cumulative ktonnesCO2-eq")
        label11.grid(column=2, row=10, sticky=tk.E, padx=5, pady=5)
        
        label12 = ttk.Label(self, text="Ave ktonnesCO2-eq per year")
        label12.grid(column=3, row=10, sticky=tk.E, padx=5, pady=5)
        
        label30 = ttk.Label(self, text="kgCO2-eq/kg.des")
        label30.grid(column=4, row=10, sticky=tk.E, padx=5, pady=5)
        label31 = ttk.Label(self, text="kgCO2-eq/kg.sev")
        label31.grid(column=5, row=10, sticky=tk.E, padx=5, pady=5)
        label32 = ttk.Label(self, text="kgCO2-eq/kg.iso")
        label32.grid(column=6, row=10, sticky=tk.E, padx=5, pady=5)

        # LHS Labels
        label1 = ttk.Label(self, text="Calculation Period")
        label1.grid(column=1, row=5, sticky=tk.W, padx=5, pady=5)
        label6 = ttk.Label(self, text="Initial Population")
        label6.grid(column=1, row=6, sticky=tk.W, padx=5, pady=5)
        label7 = ttk.Label(self, text="Desflurane % contrubution to MAC")
        label7.grid(column=1, row=7, sticky=tk.W, padx=5, pady=5)
        
        separator = ttk.Separator(self, orient='horizontal')
        separator.grid(row=8, column=1, columnspan=20, sticky="ew", padx=5, pady=5)

            
        label2 = ttk.Label(self, text="R6: Continuous Release")
        label2.grid(column=1, row=11, sticky=tk.W, padx=5, pady=5)
        label3 = ttk.Label(self, text="R7: AC capture with municipal waste incineration")
        label3.grid(column=1, row=12, sticky=tk.W, padx=5, pady=5)
        label4 = ttk.Label(self, text="R8: AC capture with pharmaceutical waste incineration")
        label4.grid(column=1, row=13, sticky=tk.W, padx=5, pady=5)
        label5 = ttk.Label(self, text="R9: AC capture with waste management companies")
        label5.grid(column=1, row=14, sticky=tk.W, padx=5, pady=5)
        label5 = ttk.Label(self, text="R10: AC capture with on-site incineration")
        label5.grid(column=1, row=15, sticky=tk.W, padx=5, pady=5)
        label5 = ttk.Label(self, text="R11: In-situ plasma destruction")
        label5.grid(column=1, row=16, sticky=tk.W, padx=5, pady=5)
        
        separator = ttk.Separator(self, orient='horizontal')
        separator.grid(row=17, column=1, columnspan=20, sticky="ew", padx=5, pady=5)
        
        # RHS Labels
        label10 = ttk.Label(self, text="Years")
        label10.grid(column=3, row=5, sticky=tk.W, padx=5, pady=5)
        label13 = ttk.Label(self, text="M Persons")
        label13.grid(column=3, row=6, sticky=tk.W, padx=5, pady=5)
        label14 = ttk.Label(self, text="%")
        label14.grid(column=3, row=7, sticky=tk.W, padx=5, pady=5)
        
        #Surgery emmisions table 
        label20 = ttk.Label(self, text="Desflurane")
        label20.grid(column=2, row=20, padx=5, pady=5)
        label21 = ttk.Label(self, text="Sevoflurane")
        label21.grid(column=3, row=20, padx=5, pady=5)
        label22 = ttk.Label(self, text="Isoflurane")
        label22.grid(column=4, row=20, padx=5, pady=5)
        label23 = ttk.Label(self, text="Total")
        label23.grid(column=5, row=20, padx=5, pady=5)
        label24 = ttk.Label(self, text="Units")
        label24.grid(column=6, row=20, sticky=tk.W, padx=5, pady=5)
        
        label25 = ttk.Label(self, text="kmols")
        label25.grid(column=6, row=21, sticky=tk.W, padx=5, pady=5)
        label26 = ttk.Label(self, text="ktonneCO2 - eq")
        label26.grid(column=6, row=22, sticky=tk.W, padx=5, pady=5)
        label27 = ttk.Label(self, text="% Share of emissions")
        label27.grid(column=6, row=23, sticky=tk.W, padx=5, pady=5)
        label27 = ttk.Label(self, text="Emissions per 100k procedures")
        label27.grid(column=1, row=20, sticky=tk.W, padx=5, pady=5)
        
        self.output10 = ttk.Label(self, text="")
        self.output10.grid(column=2, row=21)
        self.output11 = ttk.Label(self, text="")
        self.output11.grid(column=2, row=22)
        self.output12 = ttk.Label(self, text="")
        self.output12.grid(column=2, row=23)
        
        self.output13 = ttk.Label(self, text="")
        self.output13.grid(column=3, row=21)
        self.output14 = ttk.Label(self, text="")
        self.output14.grid(column=3, row=22)
        self.output15 = ttk.Label(self, text="")
        self.output15.grid(column=3, row=23)
        
        self.output16 = ttk.Label(self, text="")
        self.output16.grid(column=4, row=21)
        self.output17 = ttk.Label(self, text="")
        self.output17.grid(column=4, row=22)
        self.output18 = ttk.Label(self, text="")
        self.output18.grid(column=4, row=23)
        
        self.output19 = ttk.Label(self, text="")
        self.output19.grid(column=5, row=21)
        self.output20 = ttk.Label(self, text="")
        self.output20.grid(column=5, row=22)
        self.output21 = ttk.Label(self, text="")
        self.output21.grid(column=5, row=23)
        
        #Main Outputs
        
        self.output30 = ttk.Label(self, text="")
        self.output30.grid(column=2, row=11)
        self.output31 = ttk.Label(self, text="")
        self.output31.grid(column=3, row=11)
        
        self.output32 = ttk.Label(self, text="")
        self.output32.grid(column=2, row=12)
        self.output33 = ttk.Label(self, text="")
        self.output33.grid(column=3, row=12)
        
        self.output34 = ttk.Label(self, text="")
        self.output34.grid(column=2, row=13)
        self.output35 = ttk.Label(self, text="")
        self.output35.grid(column=3, row=13)
        
        self.output36 = ttk.Label(self, text="")
        self.output36.grid(column=2, row=14)
        self.output37 = ttk.Label(self, text="")
        self.output37.grid(column=3, row=14)
        
        self.output38 = ttk.Label(self, text="")
        self.output38.grid(column=2, row=15)
        self.output39 = ttk.Label(self, text="")
        self.output39.grid(column=3, row=15)
        
        self.output40 = ttk.Label(self, text="")
        self.output40.grid(column=2, row=16)
        self.output41 = ttk.Label(self, text="")
        self.output41.grid(column=3, row=16)
        
        #per kg outputs 
        
        self.output42 = ttk.Label(self, text="")
        self.output42.grid(column=4, row=11)
        self.output43 = ttk.Label(self, text="")
        self.output43.grid(column=5, row=11)
        self.output44 = ttk.Label(self, text="")
        self.output44.grid(column=6, row=11)
        
        self.output45 = ttk.Label(self, text="")
        self.output45.grid(column=4, row=12)
        self.output46 = ttk.Label(self, text="")
        self.output46.grid(column=5, row=12)
        self.output47 = ttk.Label(self, text="")
        self.output47.grid(column=6, row=12)
        
        self.output48 = ttk.Label(self, text="")
        self.output48.grid(column=4, row=13)
        self.output49 = ttk.Label(self, text="")
        self.output49.grid(column=5, row=13)
        self.output50 = ttk.Label(self, text="")
        self.output50.grid(column=6, row=13)
        
        self.output51 = ttk.Label(self, text="")
        self.output51.grid(column=4, row=14)
        self.output52 = ttk.Label(self, text="")
        self.output52.grid(column=5, row=14)
        self.output53 = ttk.Label(self, text="")
        self.output53.grid(column=6, row=14)
        
        self.output54 = ttk.Label(self, text="")
        self.output54.grid(column=4, row=15)
        self.output55 = ttk.Label(self, text="")
        self.output55.grid(column=5, row=15)
        self.output56 = ttk.Label(self, text="")
        self.output56.grid(column=6, row=15)
        
        self.output57 = ttk.Label(self, text="")
        self.output57.grid(column=4, row=16)
        self.output58 = ttk.Label(self, text="")
        self.output58.grid(column=5, row=16)
        self.output59 = ttk.Label(self, text="")
        self.output59.grid(column=6, row=16)

        # Buttons 
        button1 = ttk.Button(self, text="Back to Home", command=lambda: controller.show_frame(StartPage))
        button1.grid(column=5, row=1, sticky=tk.W, padx=5, pady=5)

        button3 = ttk.Button(self, text="Exit", command=quit)
        button3.grid(column=4, row=1, sticky=tk.E, padx=5, pady=5)


    def calculate(self):
        try:
            initial_pop = float(self.entry2.get())
            des_prev = float(self.entry3.get())
            gwp_selection = self.gwp_var.get()  # Get the selected GWP value
            cf4_gwp = cf4_100
            histogram_data = []
 
            # Define dictionaries and values for calculations
            gwp_categories = {
                "GWP_100": {"des": des_100, "sev": sev_100, "iso": iso_100},
            }
                       

            if gwp_selection in gwp_categories:
            
                    
                    demunic = np.random.uniform(min_demunic, max_demunic, n_samples)
                    depharma = np.random.uniform(min_depharma, max_depharma, n_samples)
                    deplasma = np.random.uniform(min_deplasma, max_deplasma, n_samples)
                    
                    munic_other = np.random.uniform(min_munic_other, max_munic_other, n_samples)
                    pharam_other = np.random.uniform(min_pharam_other, max_pharam_other, n_samples)
                    plasma_other = np.random.uniform(min_plasma_other, max_plasma_other, n_samples)
                    
                    sur_100k = np.random.uniform(min_sur_100k, max_sur_100k, n_samples)
                    growth_rate = np.random.uniform(min_growth_rate, max_growth_rate, n_samples)
                    
                    van_factor = np.random.uniform(min_van_factor, max_van_factor, n_samples)
                    hgv_factor = np.random.uniform(min_hgv_factor, max_hgv_factor, n_samples)
                    ship_factor = np.random.uniform(min_ship_factor, max_ship_factor, n_samples)
                    
                    network_distance = np.random.uniform(min_network_distance, max_network_distance, n_samples)
                    mwi_dist = np.random.uniform(min_mwi_dist, max_mwi_dist, n_samples)
                    pwi_dist = np.random.uniform(min_pwi_dist, max_pwi_dist, n_samples)

                    #2017 Baseline
                    kg_des = (ktco2e_des_base / des_100) * 1000000
                    kg_sev = (ktco2e_sev_base / sev_100) * 1000000
                    kg_iso = (ktco2e_iso_base / iso_100) * 1000000
                    
                    kmol_des = kg_des / des_gfm
                    kmol_sev = kg_sev / sev_gfm
                    kmol_iso = kg_iso / iso_gfm
                    
                    SI_des = kmol_des / MAC_des
                    SI_sev = kmol_sev / MAC_sev
                    SI_iso = kmol_iso / MAC_iso
                    
                    SI_baseline = (SI_des + SI_sev + SI_iso) / (sur_pro_base * 10) # SI/100k surgery
                    
                    ktco2e_per_SI_des = ktco2e_des_base / SI_des
                    ktco2e_per_SI_sev = ktco2e_sev_base / SI_sev
                    ktco2e_per_SI_iso = ktco2e_iso_base / SI_iso

                    #New Emmisions per 100k Sur
                    
                    SI_100k_sur_iso =  SI_iso / (sur_pro_base * 10)
                    SI_100k_sur_des = (des_prev/100) * SI_baseline
                    SI_100k_sur_sev = SI_baseline - SI_100k_sur_des - SI_100k_sur_iso
                    
                    emmisions_per_100k_des_cal = ktco2e_per_SI_des * SI_100k_sur_des
                    emmisions_per_100k_sev_cal = ktco2e_per_SI_sev * SI_100k_sur_sev
                    emmisions_per_100k_iso_cal = ktco2e_per_SI_iso * SI_100k_sur_iso 
                    emmisions_per_100k_total_calc = emmisions_per_100k_des_cal +  emmisions_per_100k_sev_cal + emmisions_per_100k_iso_cal
                    
                    
                    percent_emmisions_calc_des = (emmisions_per_100k_des_cal / emmisions_per_100k_total_calc) *100
                    percent_emmisions_calc_sev = (emmisions_per_100k_sev_cal / emmisions_per_100k_total_calc) *100
                    percent_emmisions_calc_iso = (emmisions_per_100k_iso_cal / emmisions_per_100k_total_calc) *100
                    percent_emmisions_calc_total = (percent_emmisions_calc_des + percent_emmisions_calc_sev +  percent_emmisions_calc_iso)
                    
                    des_mol_calc = SI_100k_sur_des * MAC_des
                    sev_mol_calc = SI_100k_sur_sev * MAC_sev
                    iso_mol_calc = SI_100k_sur_iso * MAC_iso
                    total_mol_calc = des_mol_calc + sev_mol_calc + iso_mol_calc
                    
                    
                    # Time 
                    
                    total_years_str = self.entry1.get()
                    if total_years_str:  # Check if the input is not empty
                        total_years = int(float(total_years_str))  # Convert the input to a float then an integer
                    else:
                        total_years = 0  # Set a default value if no input is provided


                    
                    cumulative_surgeries = 0
                    

                    
                    for year in range(1, total_years + 1):
                        population_estimate = (initial_pop * ((1 + growth_rate) ** (year)))
                        surgeries_per_year = (population_estimate / 100000) * sur_100k
                        emmisions_per_year = (surgeries_per_year * 10 ) * emmisions_per_100k_total_calc
                        cumulative_surgeries += surgeries_per_year   
                        
                        
                        #print(f"Year {year}: Population (M) = {population_estimate:.2f}")
                        #print(f"Year {year}: Surgeries per year (M)= {surgeries_per_year:.2f}")
                        #print(f"Year {year}: emmisions_per_year (ktonneCO2e)= {surgeries_per_year:.2f}")
                    #print(f"Cumulative M surgeries over {total_years} years: {cumulative_surgeries:.2f}")    
                                 

                    # Activated Carbon Costs
                    
                    mass_va_des = des_mol_calc * des_gfm
                    mass_va_sev = sev_mol_calc * sev_gfm
                    mass_va_iso = iso_mol_calc * iso_gfm
                    mass_va_tonnes_100k = (mass_va_des + mass_va_sev + mass_va_iso) /1000 #per 100k
                    mass_va_tonnes = mass_va_tonnes_100k * (cumulative_surgeries * 10)
                    
                                       
                    
                    tonnes_ac_req_100k = (des_mol_calc / mol_per_kg_des) + (sev_mol_calc / mol_per_kg_sev) + (iso_mol_calc / mol_per_kg_iso)
                    tonnes_ac_req = tonnes_ac_req_100k * (cumulative_surgeries * 10)
                    
                    perc_des_ac = ((des_mol_calc / mol_per_kg_des) / tonnes_ac_req_100k)*100
                    perc_sev_ac = ((sev_mol_calc / mol_per_kg_sev) / tonnes_ac_req_100k)*100
                    perc_iso_ac = ((iso_mol_calc / mol_per_kg_iso) / tonnes_ac_req_100k)*100
                    
                    tonnnes_coconut = tonnes_ac_req * cococnut_to_ac_ratio
                    
                    ac_production_emmisions = ac_co2e_per_kg * (tonnes_ac_req / 1000) #ktonnesCO2e
                    
                    ac_transport_1  = ((indo_to_sing_port * tonnnes_coconut * ship_factor)/1000000) + ((sing_port_to_fac * tonnnes_coconut * hgv_factor)/1000000) + ((sing_port_to_fac * tonnes_ac_req * hgv_factor)/1000000) + ((sing_to_eng_port * tonnes_ac_req * ship_factor)/1000000) #ktonnesCO2e
                    
                    ac_transport_2 = (network_distance * tonnes_ac_req * hgv_factor)/1000000
                    
                    ac_transport_3_msw = (mwi_dist * (tonnes_ac_req + mass_va_tonnes) * hgv_factor) /1000000
                    ac_transport_3_pwi = (pwi_dist * (tonnes_ac_req + mass_va_tonnes) * hgv_factor) /1000000
                    ac_transport_3_pla = (network_distance * (tonnes_ac_req + mass_va_tonnes) * hgv_factor ) /1000000
                    
                    #Combustion

                    # For des
                    combust_value_des_msw = ((((des_mol_calc * (1 - demunic)) * des_atecon_C * 44.1)) / 1000) + (((((des_mol_calc * (demunic)) * des_atecon_C * 88.0043)) / 1000) * cf4_gwp)
                    combust_emmisions_des_msw = combust_value_des_msw * (cumulative_surgeries * 10)
                    
                    combust_value_des_pwi = ((((des_mol_calc * (1 - depharma)) * des_atecon_C * 44.1)) / 1000) + (((((des_mol_calc * (depharma)) * des_atecon_C * 88.0043)) / 1000) * cf4_gwp) 
                    combust_emmisions_des_pwi = combust_value_des_pwi * (cumulative_surgeries * 10)            
                    
                    combust_value_des_pla = ((((des_mol_calc * (1 - deplasma)) * des_atecon_C * 44.1)) / 1000) + (((((des_mol_calc * (deplasma)) * des_atecon_C * 88.0043)) / 1000) * cf4_gwp)
                    combust_emmisions_des_pla = combust_value_des_pla * (cumulative_surgeries * 10)  
                    
                    #For sev
                    combust_value_sev_msw = ((((sev_mol_calc * (1 - demunic)) * sev_atecon_C * 44.1)) / 1000) + (((((sev_mol_calc * (demunic)) * sev_atecon_C * 88.0043)) / 1000) * cf4_gwp) 
                    combust_emmisions_sev_msw = combust_value_sev_msw * (cumulative_surgeries * 10)

                    combust_value_sev_pwi = ((((sev_mol_calc * (1 - depharma)) * sev_atecon_C * 44.1)) / 1000) + (((((sev_mol_calc * (depharma)) * sev_atecon_C * 88.0043)) / 1000) * cf4_gwp)  
                    combust_emmisions_sev_pwi = combust_value_sev_pwi * (cumulative_surgeries * 10)            

                    combust_value_sev_pla = ((((sev_mol_calc * (1 - deplasma)) * sev_atecon_C * 44.1)) / 1000) + (((((sev_mol_calc * (deplasma)) * sev_atecon_C * 88.0043)) / 1000) * cf4_gwp)            
                    combust_emmisions_sev_pla = combust_value_sev_pla * (cumulative_surgeries * 10)
                    
                    # For iso
                    combust_value_iso_msw = ((((iso_mol_calc * (1 - demunic)) * iso_atecon_C * 44.1)) / 1000) + (((((iso_mol_calc * (demunic)) * iso_atecon_C * 88.0043)) / 1000) * cf4_gwp)
                    combust_emmisions_iso_msw = combust_value_iso_msw * (cumulative_surgeries * 10)

                    combust_value_iso_pwi = ((((iso_mol_calc * (1 - depharma)) * iso_atecon_C * 44.1)) / 1000) + (((((iso_mol_calc * (depharma)) * iso_atecon_C * 88.0043)) / 1000) * cf4_gwp)
                    combust_emmisions_iso_pwi = combust_value_iso_pwi * (cumulative_surgeries * 10)

                    combust_value_iso_pla = ((((iso_mol_calc * (1 - deplasma)) * iso_atecon_C * 44.1)) / 1000) + (((((iso_mol_calc * (deplasma)) * iso_atecon_C * 88.0043)) / 1000) * cf4_gwp)
                    combust_emmisions_iso_pla = combust_value_iso_pla * (cumulative_surgeries * 10)

                    #Combustion totals
                    combust_msw_t = ((combust_emmisions_des_msw + combust_emmisions_sev_msw + combust_emmisions_iso_msw)/1000)
                    combust_pwi_t = ((combust_emmisions_des_pwi + combust_emmisions_sev_pwi + combust_emmisions_iso_pwi)/1000)
                    combust_pla_t = ((combust_emmisions_des_pla + combust_emmisions_sev_pla + combust_emmisions_iso_pla)/1000)
                    
                    
                    #Emmisions from kg waste incineration 
                    emmisions_per_kg_msw = (((mass_va_tonnes + tonnes_ac_req)*1000)*munic_other) / 1000000
                    emmisions_per_kg_pwi = (((mass_va_tonnes + tonnes_ac_req)*1000)*pharam_other) / 1000000
                    emmisions_per_kg_pla = (((mass_va_tonnes + tonnes_ac_req)*1000)*plasma_other) / 1000000
                    

                
                    #Calculations
                    
                    calculation_r6 = (cumulative_surgeries * 10 ) * emmisions_per_100k_total_calc
                    calculation_r6_yr = calculation_r6 / total_years
                    calculation_r6_des = (emmisions_per_100k_des_cal / mass_va_des) * 1000000
                    calculation_r6_sev = (emmisions_per_100k_sev_cal / mass_va_sev) * 1000000
                    calculation_r6_iso = (emmisions_per_100k_iso_cal / mass_va_iso) * 1000000
                    
                    
                    calculation_r7 = ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw + combust_msw_t + emmisions_per_kg_msw
                    calculation_r7_yr = calculation_r7 / total_years
                    calculation_r7_des = munic_other + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_emmisions_des_msw)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10)))
                    calculation_r7_sev = munic_other + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_emmisions_sev_msw)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10)))
                    calculation_r7_iso = munic_other + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_msw))+combust_emmisions_iso_msw)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
                    
                    
                    calculation_r8 = ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi + combust_pwi_t + emmisions_per_kg_pwi
                    calculation_r8_yr = calculation_r8 / total_years
                    calculation_r8_des = pharam_other + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_emmisions_des_pwi)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10)))
                    calculation_r8_sev = pharam_other + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_emmisions_sev_pwi)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10)))
                    calculation_r8_iso = pharam_other + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pwi))+combust_emmisions_iso_pwi)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
                    
                    calculation_r9 = ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla + combust_pla_t + emmisions_per_kg_pla
                    calculation_r9_yr = calculation_r9 / total_years
                    calculation_r9_des = plasma_other + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_emmisions_des_pla)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10)))
                    calculation_r9_sev = plasma_other + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_emmisions_sev_pla)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10)))
                    calculation_r9_iso = plasma_other + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2 + ac_transport_3_pla))+combust_emmisions_iso_pla)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
                    
                    calculation_r10 = ac_production_emmisions + ac_transport_1 + ac_transport_2 + combust_pwi_t + emmisions_per_kg_pwi
                    calculation_r10_yr = calculation_r10 / total_years
                    calculation_r10_des = pharam_other + ((((((perc_des_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2))+combust_emmisions_des_pwi)/1000)*1000000)/(mass_va_des * (cumulative_surgeries*10)))
                    calculation_r10_sev = pharam_other + ((((((perc_sev_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2))+combust_emmisions_sev_pwi)/1000)*1000000)/(mass_va_sev * (cumulative_surgeries*10)))
                    calculation_r10_iso = pharam_other + ((((((perc_iso_ac/100)*(ac_production_emmisions + ac_transport_1 + ac_transport_2))+combust_emmisions_iso_pwi)/1000)*1000000)/(mass_va_iso * (cumulative_surgeries*10)))
                    
                    calculation_r11 = combust_pla_t + emmisions_per_kg_pla
                    calculation_r11_yr = calculation_r11 / total_years
                    calculation_r11_des = plasma_other + ((((combust_emmisions_des_pla)/1000)*1000000)/((mass_va_des * (cumulative_surgeries*10))))
                    calculation_r11_sev = plasma_other + ((((combust_emmisions_sev_pla)/1000)*1000000)/((mass_va_sev * (cumulative_surgeries*10))))
                    calculation_r11_iso = plasma_other + ((((combust_emmisions_iso_pla)/1000)*1000000)/((mass_va_iso * (cumulative_surgeries*10))))
                    
                    
                    mean_r6_1 = np.mean(calculation_r6)
                    mean_r6_2 = np.mean(calculation_r6_yr)
                    mean_r6_3 = np.mean(calculation_r6_des)
                    mean_r6_4 = np.mean(calculation_r6_sev)
                    mean_r6_5 = np.mean(calculation_r6_iso)
                    
                    mean_r7_1 = np.mean(calculation_r7)
                    mean_r7_2 = np.mean(calculation_r7_yr)
                    mean_r7_3 = np.mean(calculation_r7_des)
                    mean_r7_4 = np.mean(calculation_r7_sev)
                    mean_r7_5 = np.mean(calculation_r7_iso)
                    
                    mean_r8_1 = np.mean(calculation_r8)
                    mean_r8_2 = np.mean(calculation_r8_yr)
                    mean_r8_3 = np.mean(calculation_r8_des)
                    mean_r8_4 = np.mean(calculation_r8_sev)
                    mean_r8_5 = np.mean(calculation_r8_iso)
                    
                    mean_r9_1 = np.mean(calculation_r9)
                    mean_r9_2 = np.mean(calculation_r9_yr)
                    mean_r9_3 = np.mean(calculation_r9_des)
                    mean_r9_4 = np.mean(calculation_r9_sev)
                    mean_r9_5 = np.mean(calculation_r9_iso)
                    
                    mean_r10_1 = np.mean(calculation_r10)
                    mean_r10_2 = np.mean(calculation_r10_yr)
                    mean_r10_3 = np.mean(calculation_r10_des)
                    mean_r10_4 = np.mean(calculation_r10_sev)
                    mean_r10_5 = np.mean(calculation_r10_iso)
                    
                    mean_r11_1 = np.mean(calculation_r11)
                    mean_r11_2 = np.mean(calculation_r11_yr)
                    mean_r11_3 = np.mean(calculation_r11_des)
                    mean_r11_4 = np.mean(calculation_r11_sev)
                    mean_r11_5 = np.mean(calculation_r11_iso)
                    
                    
                    mean_d_1 = np.mean(des_mol_calc)
                    mean_d_2 = np.mean(emmisions_per_100k_des_cal)
                    mean_d_3 = np.mean(percent_emmisions_calc_des)
                                       
                    mean_s_1 = np.mean(sev_mol_calc)
                    mean_s_2 = np.mean(emmisions_per_100k_sev_cal)
                    mean_s_3 = np.mean(percent_emmisions_calc_sev)
                    
                    mean_i_1 = np.mean(iso_mol_calc)
                    mean_i_2 = np.mean(emmisions_per_100k_iso_cal)
                    mean_i_3 = np.mean(percent_emmisions_calc_iso)
                    
                    mean_t_1 = np.mean(total_mol_calc)
                    mean_t_2 = np.mean(emmisions_per_100k_total_calc)
                    mean_t_3 = np.mean(percent_emmisions_calc_total)
                    
                    sd_cv_6d = np.std(calculation_r6_des)
                    sd_cv_6s = np.std(calculation_r6_sev)
                    sd_cv_6i = np.std(calculation_r6_iso)
                    
                    sd_cv_7d = np.std(calculation_r7_des)
                    sd_cv_7s = np.std(calculation_r7_sev)
                    sd_cv_7i = np.std(calculation_r7_iso)
                    
                    sd_cv_8d = np.std(calculation_r8_des)
                    sd_cv_8s = np.std(calculation_r8_sev)
                    sd_cv_8i = np.std(calculation_r8_iso)
                    
                    sd_cv_9d = np.std(calculation_r9_des)
                    sd_cv_9s = np.std(calculation_r9_sev)
                    sd_cv_9i = np.std(calculation_r9_iso)
                    
                    sd_cv_10d = np.std(calculation_r10_des)
                    sd_cv_10s = np.std(calculation_r10_sev)
                    sd_cv_10i = np.std(calculation_r10_iso)
                    
                    sd_cv_11d = np.std(calculation_r11_des)
                    sd_cv_11s = np.std(calculation_r11_sev)
                    sd_cv_11i = np.std(calculation_r11_iso)
                    
  
                                       
                                    
                    
                  
                    # Display calculated values in respective output labels
                    
                    self.output30.config(text=f'{mean_r6_1:.1f}')
                    self.output31.config(text=f'{mean_r6_2:.1f}')
                    self.output42.config(text=f'{mean_r6_3:.1f}')
                    self.output43.config(text=f'{mean_r6_4:.1f}')
                    self.output44.config(text=f'{mean_r6_5:.1f}')
                    
                    
                    self.output32.config(text=f'{mean_r7_1:.1f}')
                    self.output33.config(text=f'{mean_r7_2:.1f}')
                    self.output45.config(text=f'{mean_r7_3:.1f}')
                    self.output46.config(text=f'{mean_r7_4:.1f}')
                    self.output47.config(text=f'{mean_r7_5:.1f}')
                    
                    self.output34.config(text=f'{mean_r8_1:.1f}')
                    self.output35.config(text=f'{mean_r8_2:.1f}')
                    self.output48.config(text=f'{mean_r8_3:.1f}')
                    self.output49.config(text=f'{mean_r8_4:.1f}')
                    self.output50.config(text=f'{mean_r8_5:.1f}')
                    
                    self.output36.config(text=f'{mean_r9_1:.1f}')
                    self.output37.config(text=f'{mean_r9_2:.1f}')
                    self.output51.config(text=f'{mean_r9_3:.1f}')
                    self.output52.config(text=f'{mean_r9_4:.1f}')
                    self.output53.config(text=f'{mean_r9_5:.1f}')
                    
                    self.output38.config(text=f'{mean_r10_1:.1f}')
                    self.output39.config(text=f'{mean_r10_2:.1f}')
                    self.output54.config(text=f'{mean_r10_3:.1f}')
                    self.output55.config(text=f'{mean_r10_4:.1f}')
                    self.output56.config(text=f'{mean_r10_5:.1f}')
                    
                    self.output40.config(text=f'{mean_r11_1:.1f}')
                    self.output41.config(text=f'{mean_r11_2:.1f}')
                    self.output57.config(text=f'{mean_r11_3:.1f}')
                    self.output58.config(text=f'{mean_r11_4:.1f}')
                    self.output59.config(text=f'{mean_r11_5:.1f}')
                    
                    #Display emmisions from surgey

                    self.output10.config(text=f'{mean_d_1:.2f}')
                    self.output11.config(text=f'{mean_d_2:.2f}')
                    self.output12.config(text=f'{mean_d_3:.2f}')
                    
                    self.output13.config(text=f'{mean_s_1:.2f}')
                    self.output14.config(text=f'{mean_s_2:.2f}')
                    self.output15.config(text=f'{mean_s_3:.2f}')
                    
                    self.output16.config(text=f'{mean_i_1:.2f}')
                    self.output17.config(text=f'{mean_i_2:.2f}')
                    self.output18.config(text=f'{mean_i_3:.2f}')
                    
                    self.output19.config(text=f'{mean_t_1:.1f}')
                    self.output20.config(text=f'{mean_t_2:.1f}')
                    self.output21.config(text=f'{mean_t_3:.1f}')
                    
                    print(sd_cv_6d)
                    print(sd_cv_6s)
                    print(sd_cv_6i)
                    
                    print(sd_cv_7d)
                    print(sd_cv_7s)
                    print(sd_cv_7i)
                    
                    print(sd_cv_8d)
                    print(sd_cv_8s)
                    print(sd_cv_8i)
                    
                    print(sd_cv_9d)
                    print(sd_cv_9s)
                    print(sd_cv_9i)
                    
                    print(sd_cv_10d)
                    print(sd_cv_10s)
                    print(sd_cv_10i)
                    
                    print(sd_cv_11d)
                    print(sd_cv_11s)
                    print(sd_cv_11i)
                    
                    
                    
 

        except ValueError:
            pass
        # Update GUI
        self.update_idletasks()  # Ensure the GUI is updated


 
      


app = vaco2eapp()
app.geometry("1280x720")

for child in app.winfo_children():
    child.grid_configure(padx=5, pady=5)

app.mainloop()

