# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 19:52:53 2024

@author: young
"""

import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from tkinter import Tk, Label, PhotoImage
import PIL.Image
import PIL.ImageTk
from PIL import ImageTk, Image
""" 
make your Frame class

"self" means the your class instance  
Thus, to put anything into the frame use "self" !! 

for example, put a label into your frame use 

   self.new_label = tk.Label(self, text='new label')
   self.new_label.pack()
    
Always define some variables as "self.something"           
"""
class your_Frame(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        """ Add your widgets here """        
        #---simple text label        
        self.label = tk.Label(self, text = "This is a test ")
        self.label.pack()
        #---simple image label         
        img = ImageTk.PhotoImage(Image.open("test.png"))           
        self.label2 = tk.Label(self, image = img)
        self.label2.image = img #this is necessary to keep image reference 
        self.label2.pack() 
        
    def do_something(self,):
        return 
        
"""
Below line is only active when it is run as a independent program  
"""
if __name__=='__main__':
    def main():
        root = tk.Tk()
        frame1 = your_Frame(root)
        frame1.pack()
        root.mainloop()
        return     
    main() 