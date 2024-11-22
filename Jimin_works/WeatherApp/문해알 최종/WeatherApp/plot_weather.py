# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:23:08 2024

@author: young
"""
import tkinter as tk
from tkinter import Canvas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import ImageTk, Image
import os

# Sample data dictionary
weather_data = {
    "08:00": {"temperature": 15, "humidity": 80, "sky": 1, "precipitation_type": 0, "precipitation": 0},
    "12:00": {"temperature": 20, "humidity": 70, "sky": 3, "precipitation_type": 0, "precipitation": 0},
    "16:00": {"temperature": 18, "humidity": 75, "sky": 4, "precipitation_type": 1, "precipitation": 5},
    "20:00": {"temperature": 12, "humidity": 85, "sky": 4, "precipitation_type": 1, "precipitation": 10},
}

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
def load_images(size=(50, 50)):
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


# weather_data={"08:00": {"temperature": 15, "humidity": 80, "sky": 1, "precipitation_type": 0, "precipitation": 0},} 
#      
def get_lists_from_information(informations):
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
        
def create_weather_frame(parent, weather_data, images):
    # Extract data for plotting and displaying
    # times = list(weather_data.keys())
    # temperatures = [weather_data[time]["temperature"] for time in times]
    # humidities = [weather_data[time]["humidity"] for time in times]
    # skies = [weather_data[time]["sky"] for time in times]
    # precip_types = [weather_data[time]["precipitation_type"] for time in times]
    # precipitations = [weather_data[time]["precipitation"] for time in times]
    
    times,temperatures,humidities,skies,precip_types,precipitations=get_lists_from_information(weather_data)

    frame = tk.Frame(parent)

    # Frame for weather icons
    icon_frame = tk.Frame(frame, height=100)
    icon_frame.pack(fill=tk.X)

    for i, time in enumerate(times):
        sub_frame = tk.Frame(icon_frame, width=100)
        sub_frame.pack(side=tk.LEFT, padx=10)

        # Weather icon
        sky_label = sky_code.get(skies[i], "")
        precip_label = precipitation_code.get(precip_types[i], "")
        icon = images.get(sky_label, images.get(precip_label, None))

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

    canvas_temp = FigureCanvasTkAgg(fig_temp, frame)
    canvas_temp.get_tk_widget().pack()

    # Humidity plot
    fig_humidity, ax_humidity = plt.subplots(figsize=(6, 2))
    ax_humidity.plot(times, humidities, marker="o", color="blue")
    ax_humidity.set_title("Humidity vs Time")
    ax_humidity.set_xlabel("Time")
    ax_humidity.set_ylabel("Humidity (%)")
    ax_humidity.grid()

    canvas_humidity = FigureCanvasTkAgg(fig_humidity, frame)
    canvas_humidity.get_tk_widget().pack()

    return frame

# Example usage
if __name__=='__main__':
    def create_gui():
        root = tk.Tk()
        root.title("Weather Visualization")
    
        images = load_images()
        weather_frame = create_weather_frame(root, weather_data, images)
        weather_frame.pack()
    
        root.mainloop()
    
    
    create_gui()
