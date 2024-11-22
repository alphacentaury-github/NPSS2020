# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 12:09:12 2024

@author: young
"""
import tkinter as tk
from tkinter import ttk
import pandas as pd 
import json 
from tkinter import Tk, Label, PhotoImage
import PIL.Image
import PIL.ImageTk
from PIL import ImageTk, Image

from get_weather import get_weather2 
from get_weather import extract_loc_name_dictionary 

# sample_data = {
#     "서울특별시": {
#         "강남구": ["삼성동", "대치동", "역삼동"],
#         "서초구": ["서초동", "반포동", "잠원동"],
#     },
#     "부산광역시": {
#         "해운대구": ["우동", "중동", "좌동"],
#         "사하구": ["하단동", "당리동", "괴정동"],
#     },
#     "대구광역시": {
#         "수성구": ["범어동", "수성동", "만촌동"],
#         "달서구": ["상인동", "월성동", "진천동"],
#     },
# }
    
sample_data = extract_loc_name_dictionary()

# LocationSelector 클래스
class LocationSelector(tk.Frame):
    def __init__(self, parent, location_data=sample_data, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.location_data = location_data

        self.container_combo = tk.Frame(self)
        self.container_combo.pack() 
        
        # 1단계: 시
        tk.Label(self.container_combo, text="시").grid(row=0, column=0, padx=10, pady=5)
        self.city_combobox = ttk.Combobox(self.container_combo, values=list(location_data.keys()), state="readonly")
        self.city_combobox.set('서울특별시')
        self.city_combobox.grid(row=0, column=1, padx=10, pady=5)
        self.city_combobox.bind("<<ComboboxSelected>>", self.update_gu)

        # 2단계: 구
        tk.Label(self.container_combo, text="구").grid(row=1, column=0, padx=10, pady=5)
        self.gu_combobox = ttk.Combobox(self.container_combo, state="readonly")
        self.gu_combobox.set('종로구')
        self.gu_combobox.grid(row=1, column=1, padx=10, pady=5)
        self.gu_combobox.bind("<<ComboboxSelected>>", self.update_dong)

        # 3단계: 동
        tk.Label(self.container_combo, text="동").grid(row=2, column=0, padx=10, pady=5)
        self.dong_combobox = ttk.Combobox(self.container_combo, state="readonly")
        self.dong_combobox.set('혜화동')
        self.dong_combobox.grid(row=2, column=1, padx=10, pady=5)
        
        self.weather_button = tk.Button(self, text="Get Weather", command= lambda: self.show_weather())
        self.weather_button.pack() 

        self.container_plot = tk.Frame(self) 
        self.container_plot.pack() 
        
        self.label2 = tk.Label(self,text='You can use unicode like'+u'\u2600',font=('Arial', 16, 'bold')) #unicode
        self.label2.pack()  
        self.label3 = tk.Label(self)
        self.label3.pack()
        
        
        self.show_weather() 


    def update_gu(self, event):
        city = self.city_combobox.get()
        if city:
            self.gu_combobox["values"] = list(self.location_data[city].keys())
            self.gu_combobox.set("")
            self.dong_combobox["values"] = []
            self.dong_combobox.set("")

    def update_dong(self, event):
        city = self.city_combobox.get()
        gu = self.gu_combobox.get()
        if city and gu:
            self.dong_combobox["values"] = self.location_data[city][gu]
            self.dong_combobox.set("")

    def show_weather(self,):        
        city = self.city_combobox.get()
        gu = self.gu_combobox.get()
        dong = self.dong_combobox.get() 
        informations, list_info_text = get_weather2([city,gu,dong],opt_print=False)
        # plot weather and create image 
        import plot_weather2 as pw 
        try: #destroy existing plot 
            self.weather_frame.destroy() 
        except: 
            pass              
        self.weather_frame = pw.plot_weather_Frame(self)
        self.weather_frame.pack()
        self.weather_frame.update_figure(informations)
        self.label2.configure(text = list_info_text[0][14:51])
        
        return 
            
if __name__=='__main__':
    informations, list_info_text = get_weather2(['대전광역시','유성구','신성동'],opt_print=False)
    root = tk.Tk()
    frame1 = LocationSelector(root)
    frame1.pack()
    tk.mainloop()            