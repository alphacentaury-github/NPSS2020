# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 12:45:39 2024

@author: young
"""

import tkinter as tk
from tkinter import ttk

# LocationSelector 클래스
from LocationSelector import LocationSelector

def show_frame(frame):
    """컨테이너에 새 프레임 표시
    container_frame: 부모 frame 
    """
    for widget in container_frame.winfo_children():
        widget.pack_forget()  # 기존의 모든 위젯 숨기기
    frame.pack(fill=tk.BOTH, expand=True)  # 새 프레임 표시

# 메인 창 생성
root = tk.Tk()
root.geometry("600x400")

# 컨테이너 프레임 (내용이 바뀔 영역)
container_frame = tk.Frame(root, bg="white")
container_frame.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

# 버튼 패널 (프레임 전환 버튼)
button_panel = tk.Frame(root, bg="lightgray")
button_panel.pack(fill=tk.Y, side=tk.RIGHT)

# 여러 프레임 정의


frame1 = LocationSelector(container_frame)

from test_graph import plot_Frame
frame2 = plot_Frame(container_frame)
frame2.pack() 

frame3 = tk.Frame(container_frame, bg="red")
tk.Label(frame3, text="This is Frame 3", bg="red", fg="white").pack(pady=20)

# 프레임 전환 버튼
tk.Button(button_panel, text="위치/날씨 정보", command=lambda: show_frame(frame1)).pack(fill=tk.X, padx=5, pady=5)
tk.Button(button_panel, text="Show Frame 2", command=lambda: show_frame(frame2)).pack(fill=tk.X, padx=5, pady=5)
tk.Button(button_panel, text="Show Frame 3", command=lambda: show_frame(frame3)).pack(fill=tk.X, padx=5, pady=5)

# 초기 프레임 표시
show_frame(frame1)

root.mainloop()

#-------------------------




