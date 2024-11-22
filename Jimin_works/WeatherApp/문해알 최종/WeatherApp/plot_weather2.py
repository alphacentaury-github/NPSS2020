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
from tkinter import Canvas
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

from get_weather import get_weather2 
# Code mappings
precipitation_code = {
    0: "강수 없음",
    1: "비",
    2: "비/눈",
    3: "눈",
    5: "빗방울",
    6: "진눈깨비",
    7: "눈날림"
}

sky_code = {
    1: "맑음",
    3: "구름많음",
    4: "흐림"
}

# Image mappings (update with actual image file paths)
image_paths = {
    "맑음": "icons/5-removebg-preview50x50.png",
    "구름많음": "icons/1-removebg-preview50x50.png",
    "흐림": "icons/1-removebg-preview50x50.png",
    "비": "icons/3-removebg-preview50x50.png",
    "비/눈": "icons/8-removebg-preview50x50.png",
    "눈": "icons/8-removebg-preview50x50.png",
    "빗방울": "icons/8-removebg-preview50x50.png",
    "진눈깨비": "icons/8-removebg-preview50x50.png",
    "눈날림": "icons/8-removebg-preview50x50.png",
}

# Function to load images
def load_images():
    images = {}
    for key, path in image_paths.items():
        try:
            images[key] = ImageTk.PhotoImage(Image.open(path))            
        except Exception as e:
            print(f"Error loading image for {key}: {e}")
    return images

# 이미지 크기 조정 함수
def resize_image(input_path, output_path, size=(50, 50)):
    with Image.open(input_path) as img:
        resized_img = img.resize(size, Image.Resampling.LANCZOS)
        resized_img.save(output_path)
        print(f"Resized {input_path} -> {output_path}")

def get_lists_from_information(informations):
    """ 
    change informations from get_weather to lists
    """
    times = [] 
    temperatures = []
    humidities = [] 
    skies = []
    precip_types = []
    precipitations = []
    full_keys = list(informations.keys())
    if len(full_keys) > 6:
        full_keys = full_keys[:6]
    else:
        full_keys = full_keys[:-1]
        
    for key in  full_keys :
        if ('TMP' in informations[key].keys()
            and 'REH' in informations[key].keys()
            and 'SKY' in informations[key].keys()
            and 'PTY' in informations[key].keys() ):
            pass 
        else:
            print(f'wrong? at {key}')
            raise ValueError()
        try:
            new_time = f'{key[1][:2]}:{key[1][2:]}' #f'{key[0][-2:]} {key[1][:2]}:{key[1][2:]}'
            times.append(new_time)
            
            temperature = float(informations[key]['TMP'][0]) #TMP	1시간 기온	℃
                                          #TMN	일 최저기온	℃
                                          #TMX	일 최고기온	℃
            temperatures.append(temperature)
                              
            humidity = float(informations[key]['REH']) #REH	습도	%
            humidities.append(humidity)
            
            sky = int(informations[key]['SKY']) #SKY	하늘상태	코드값
            skies.append(sky)                  
            #int(informations[key]['POP']) #POP	강수확률	%
            
            precipication_type = int(informations[key]['PTY']) #PTY	강수형태	코드값
            precip_types.append(precipication_type)
            
            if precipication_type >0 :
                precipication = int(informations[key]['PCP']) #PCP	1시간 강수량	범주 (1 mm)
            else:
                precipication = 0    
            #int(informations[key]['WSD']) #WSD	풍속	m/s
            precipitations.append(precipication)
            
        except:
            break 
    return times,temperatures,humidities,skies,precip_types,precipitations         


class plot_weather_Frame(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        """ Add your widgets here """        
        self.container = tk.Frame(self)
        self.container.pack() 

        self.images = load_images() #load all icon images 

        self.icon_frame = tk.Frame(self.container, height=100) # for icons
        self.icon_frame.pack(fill=tk.X)
        self.canvas_frame = tk.Frame(self.container, height=100) # for icons
        self.canvas_frame.pack(fill=tk.X)


    def update_figure(self,weather_data):
        times,temperatures,humidities,skies,precip_types,precipitations=get_lists_from_information(weather_data)        
        #---forget----
        for widget in self.container.winfo_children():
            widget.pack_forget()  # 기존의 모든 위젯 숨기기
        self.container.pack(fill=tk.BOTH, expand=True)
        #---destory old plots         
        self.icon_frame.destroy() 
        self.icon_frame = tk.Frame(self.container, height=100) # for icons
        self.icon_frame.pack(fill=tk.X)
        self.canvas_frame.destroy()
        self.canvas_frame = tk.Frame(self.container, height=100) # for icons
        self.canvas_frame.pack(fill=tk.X)
        #---make new plots and pack  
        for i, time in enumerate(times):
            sub_frame = tk.Frame(self.icon_frame, width=100)
            sub_frame.pack(side=tk.LEFT, padx=10)

            # Weather icon
            sky_label = sky_code.get(skies[i], "")
            precip_label = precipitation_code.get(precip_types[i], "")
            icon = self.images.get(sky_label, self.images.get(precip_label, None))

            if icon:
                icon_label = tk.Label(sub_frame, image=icon)
                icon_label.image = icon  # Keep reference to avoid garbage collection
                icon_label.pack()
            # Time
            time_label = tk.Label(sub_frame, text=time, font=("Arial", 12))
            time_label.pack()

            # Precipitation (if any)
            if precipitations[i] > 0:
                precip_text = f"{precipitations[i]} mm"
                precip_amount_label = tk.Label(sub_frame, text=precip_text, font=("Arial", 10), fg="blue")
                precip_amount_label.pack()
        # Temperature plot
        try:
            plt.close() 
        except:
            pass 
        fig_temp, ax_temp = plt.subplots(figsize=(6, 2))
        ax_temp.plot(times, temperatures, marker="o", color="red")
        ax_temp.set_title("Temperature vs Time")
        ax_temp.set_xlabel("Time")
        ax_temp.set_ylabel("Temperature (°C)")
        ax_temp.grid()

        canvas_temp = FigureCanvasTkAgg(fig_temp, self.canvas_frame)
        canvas_temp.get_tk_widget().pack()

        # Humidity plot
        fig_humidity, ax_humidity = plt.subplots(figsize=(6, 2))
        ax_humidity.plot(times, humidities, marker="o", color="blue")
        ax_humidity.set_title("Humidity vs Time")
        ax_humidity.set_xlabel("Time")
        ax_humidity.set_ylabel("Humidity (%)")
        ax_humidity.grid()

        canvas_humidity = FigureCanvasTkAgg(fig_humidity, self.canvas_frame)
        canvas_humidity.get_tk_widget().pack()
        
        
        return 
#=============================================================================
if __name__=='__main__':
    informations, list_info_text = get_weather2(['대전광역시','유성구','신성동'],opt_print=False)    
    def main():
        root = tk.Tk()
        frame1 = plot_weather_Frame(root)
        frame1.pack()
        frame1.update_figure(informations)
        root.mainloop()
        return     
    main() 