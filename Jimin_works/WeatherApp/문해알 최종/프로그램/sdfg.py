import requests
from bs4 import BeautifulSoup
import tkinter as tk
from tkinter import messagebox

# 초기 옷 데이터
clothes = {
    "상의": ["코트", "니트", "셔츠"],
    "하의": ["청바지", "두꺼운 바지", "반바지"]
}

# 사용자 위치 정보를 불러오는 함수
def get_user_location():
    return "서울"  # 기본값: 서울

# 네이버에서 날씨 정보를 스크래핑하는 함수
def get_current_weather(location):
    url = f"https://search.naver.com/search.naver?query={location}+날씨"
    headers = {'User-Agent': 'Mozilla/5.0'}
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, 'html.parser')

        # 온도 데이터 추출
        temp_tag = soup.find('div', class_='temperature_text')
        weather_tag = soup.find('p', class_='summary')

        # 온도 문자열 처리
        if temp_tag:
            temp_text = temp_tag.text.strip().replace("현재 온도", "").replace("°", "").strip()
            temperature = float(temp_text)  # 변환 가능하도록 문자열 정리
        else:
            temperature = None

        # 날씨 상태
        weather_status = weather_tag.text.strip() if weather_tag else "알 수 없음"

        return temperature, weather_status
    except Exception as e:
        print(f"날씨 정보를 가져오는 데 실패했습니다: {e}")
        return None, "데이터 없음"

# 날씨를 분석하는 함수
def analyse_weather(temp):
    if temp <= 5:
        return "추운 날씨", "두꺼운 옷을 추천합니다."
    elif temp <= 15:
        return "쌀쌀한 날씨", "겉옷을 준비하세요."
    elif temp <= 25:
        return "적당한 날씨", "가벼운 옷차림이 좋습니다."
    else:
        return "더운 날씨", "시원한 옷을 입으세요."

# 추가할 옷을 표시하는 함수
def check_clothes():
    return clothes

# 옷 추천을 생성하는 함수
def generate_recommendation(temp):
    if temp <= 5:
        return "코트", "두꺼운 바지"
    elif temp <= 15:
        return "니트", "청바지"
    elif temp <= 25:
        return "셔츠", "청바지"
    else:
        return "셔츠", "반바지"

# 위치와 날씨를 새로고침하는 함수
def regen():
    location = get_user_location()
    temp, weather_status = get_current_weather(location)
    return temp, weather_status

# 상의를 추천하는 함수
def recommend_top(temp):
    top, _ = generate_recommendation(temp)
    return top

# 하의를 추천하는 함수
def recommend_bottom(temp):
    _, bottom = generate_recommendation(temp)
    return bottom

# 옷을 추가하는 함수
def add_clothes(type_, item):
    if type_ in clothes:
        clothes[type_].append(item)
        return f"{item}이(가) {type_}에 추가되었습니다."
    else:
        return "유효하지 않은 옷 종류입니다."

# 옷을 삭제하는 함수
def delete_clothes(type_, item):
    if type_ in clothes and item in clothes[type_]:
        clothes[type_].remove(item)
        return f"{item}이(가) {type_}에서 삭제되었습니다."
    else:
        return "옷을 찾을 수 없습니다."

# Tkinter 기반 UI
def main_ui():
    def update_weather():
        temp, weather_status = regen()
        if temp is not None:
            analysis, suggestion = analyse_weather(temp)
            top = recommend_top(temp)
            bottom = recommend_bottom(temp)

            temp_label.config(text=f"현재 기온: {temp}°C")
            status_label.config(text=f"날씨 상태: {weather_status}")
            analysis_label.config(text=f"날씨 분석: {analysis}")
            suggestion_label.config(text=f"추천 사항: {suggestion}")
            top_label.config(text=f"추천 상의: {top}")
            bottom_label.config(text=f"추천 하의: {bottom}")
        else:
            messagebox.showerror("에러", "날씨 정보를 가져오지 못했습니다.")

    def manage_clothes(action):
        item = clothes_entry.get()
        type_ = clothes_type.get()
        if action == "add":
            result = add_clothes(type_, item)
        elif action == "delete":
            result = delete_clothes(type_, item)
        else:
            result = "알 수 없는 동작입니다."
        messagebox.showinfo("결과", result)
        clothes_label.config(text=f"현재 옷: {clothes}")

    root = tk.Tk()
    root.title("날씨 맞춤 옷 추천 프로그램")
    root.geometry("500x500")

    # 날씨 정보
    temp_label = tk.Label(root, text="현재 기온: ---", font=("Arial", 14))
    temp_label.pack(pady=10)

    status_label = tk.Label(root, text="날씨 상태: ---", font=("Arial", 12))
    status_label.pack(pady=5)

    analysis_label = tk.Label(root, text="날씨 분석: ---", font=("Arial", 12))
    analysis_label.pack(pady=5)

    suggestion_label = tk.Label(root, text="추천 사항: ---", font=("Arial", 12))
    suggestion_label.pack(pady=5)

    top_label = tk.Label(root, text="추천 상의: ---", font=("Arial", 12))
    top_label.pack(pady=10)

    bottom_label = tk.Label(root, text="추천 하의: ---", font=("Arial", 12))
    bottom_label.pack(pady=10)

    # 옷 관리
    clothes_type = tk.StringVar(value="상의")
    tk.OptionMenu(root, clothes_type, *clothes.keys()).pack()

    clothes_entry = tk.Entry(root)
    clothes_entry.pack()

    add_button = tk.Button(root, text="옷 추가", command=lambda: manage_clothes("add"))
    add_button.pack(pady=5)

    delete_button = tk.Button(root, text="옷 삭제", command=lambda: manage_clothes("delete"))
    delete_button.pack(pady=5)

    clothes_label = tk.Label(root, text=f"현재 옷: {clothes}", font=("Arial", 10))
    clothes_label.pack(pady=10)

    # 새로고침 버튼
    refresh_button = tk.Button(root, text="날씨 정보 갱신", command=update_weather)
    refresh_button.pack(pady=20)

    # 초기 실행
    update_weather()

    root.m
