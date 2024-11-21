# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 19:52:53 2024

@author: young
"""

import matplotlib.pyplot as plt

def plot_weather(data):
    """
    날씨 정보를 바탕으로 시간대별 온도, 습도, 날씨 상태를 시각화합니다.
    
    :param data: 날씨 정보 리스트. 예: [[시간, 날씨 상태, 온도, 습도], ...]
    """
    # 데이터 분리
    times = [entry[0] for entry in data]
    weather_states = [entry[1] for entry in data]
    temperatures = [int(entry[2].replace("온도 ", "").replace("도", "")) for entry in data]
    humidities = [int(entry[3].replace("습도 ", "").replace("%", "")) for entry in data]
    
    # 그래프 그리기
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # 온도와 습도
    ax1.set_xlabel("시간")
    ax1.set_ylabel("온도 (°C)", color="tab:red")
    ax1.plot(times, temperatures, color="tab:red", marker="o", label="온도")
    ax1.tick_params(axis="y", labelcolor="tab:red")
    ax1.legend(loc="upper left")
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("습도 (%)", color="tab:blue")
    ax2.plot(times, humidities, color="tab:blue", marker="o", linestyle="--", label="습도")
    ax2.tick_params(axis="y", labelcolor="tab:blue")
    ax2.legend(loc="upper right")
    
    # 날씨 상태 추가
    for i, weather in enumerate(weather_states):
        ax1.text(times[i], temperatures[i], f"{weather}", fontsize=9, ha="center", va="bottom")
    
    # 그래프 제목
    plt.title("시간대별 날씨 정보")
    plt.grid(axis="x", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()

# 테스트 데이터
weather_data = [
    ["12시", "흐림", "온도 18도", "습도 20%"],
    ["1시", "비", "온도 15도", "습도 50%"],
    ["2시", "맑음", "온도 20도", "습도 30%"],
    ["3시", "흐림", "온도 17도", "습도 40%"]
]

plot_weather(weather_data)
