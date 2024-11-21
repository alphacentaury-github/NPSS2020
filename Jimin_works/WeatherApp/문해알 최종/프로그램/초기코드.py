import tkinter as tk
from tkinter import messagebox
import requests
from bs4 import BeautifulSoup

# 초기 옷 데이터
clothes = {
    "상의": ["코트", "니트", "셔츠"],
    "하의": ["청바지", "두꺼운 바지", "반바지"]
}

# 네이버에서 날씨 정보를 스크래핑하는 함수
def get_current_weather(location):
    url = f"https://search.naver.com/search.naver?query={location}+날씨"
    headers = {'User-Agent': 'Mozilla/5.0'}
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, 'html.parser')
        temp_tag = soup.find('div', class_='temperature_text')
        weather_tag = soup.find('p', class_='summary')

        temp = temp_tag.text.strip().replace("현재 온도", "").replace("°", "").strip() if temp_tag else "N/A"
        weather_status = weather_tag.text.strip() if weather_tag else "알 수 없음"

        return temp, weather_status
    except Exception as e:
        print(f"날씨 정보를 가져오는 데 실패했습니다: {e}")
        return "N/A", "데이터 없음"

# 날씨를 분석하는 함수
def analyse_weather(temp):
    if temp <= 5:
        return "추운 날씨", "두껍게 입으세요"
    elif temp <= 15:
        return "쌀쌀한 날씨", "겉옷을 준비하세요."
    elif temp <= 25:
        return "적당한 날씨", "가벼운 옷차림이 좋습니다."
    else:
        return "더운 날씨", "시원한 옷을 입으세요."

# 상의 추천
def recommend_top(temp):
    if temp <= 5:
        return "코트"
    elif temp <= 15:
        return "니트"
    elif temp <= 25:
        return "셔츠"
    else:
        return "셔츠"

# 하의 추천
def recommend_bottom(temp):
    if temp <= 5:
        return "두꺼운 바지"
    elif temp <= 15:
        return "청바지"
    elif temp <= 25:
        return "청바지"
    else:
        return "반바지"

# Tkinter 기반 UI
def main_ui():
    def update_weather():
        """현재 기온과 날씨 상태를 다시 불러옵니다."""
        temp, weather_status = get_current_weather("서울")
        if temp != "N/A":
            # 분석 및 추천 업데이트
            analysis, suggestion = analyse_weather(float(temp))
            top = recommend_top(float(temp))
            bottom = recommend_bottom(float(temp))

            # UI 업데이트
            temp_label.config(text=f"현재 기온: {temp}°C")
            top_label.config(text=f"추천 상의: {top}")
            bottom_label.config(text=f"추천 하의: {bottom}")
        else:
            messagebox.showerror("에러", "날씨 정보를 가져오지 못했습니다.")

    def show_help():
        """도움말 팝업 창"""
        help_window = tk.Toplevel(root)
        help_window.title("도움말")
        help_window.geometry("500x400")

        tk.Label(help_window, text="도움말", font=("Arial", 18, "bold")).pack(pady=10)

        help_text = (
            "1. 프로그램 설명\n"
            "- 본 프로그램은 기온에 따라 적절한 옷을 추천하고 다양한 기능을 제공합니다.\n\n"
            "2. 아이콘 설명\n"
            "- ?: 도움말\n"
            "- 새로고침: 날씨 정보 갱신\n\n"
            "3. 옷장 추가\n"
            "- 새로운 옷을 추가할 수 있습니다.\n\n"
            "4. 더보기\n"
            "- 다양한 기능을 확인할 수 있습니다."
        )
        tk.Label(help_window, text=help_text, font=("Arial", 12), justify="left").pack(pady=10, padx=20)

        tk.Button(help_window, text="닫기", command=help_window.destroy).pack(pady=20)


    def open_more():
        """더보기 메뉴"""
        more_window = tk.Toplevel(root)
        more_window.title("더보기")
        more_window.geometry("400x400")

        tk.Label(more_window, text="더보기 메뉴", font=("Arial", 18, "bold")).pack(pady=10)

        def show_wardrobe():
            #내 옷장 보기
            wardrobe_window = tk.Toplevel(more_window)
            wardrobe_window.title("내 옷장 보기")
            wardrobe_window.geometry("400x400")

            tk.Label(wardrobe_window, text="내 옷장 보기", font=("Arial", 18)).pack(pady=10)

            tk.Label(wardrobe_window, text="상의 목록:", font=("Arial", 14)).pack(anchor="w", padx=20)
            for top in clothes["상의"]:
                tk.Label(wardrobe_window, text=f"- {top}", font=("Arial", 12)).pack(anchor="w", padx=40)

            tk.Label(wardrobe_window, text="하의 목록:", font=("Arial", 14)).pack(anchor="w", padx=20, pady=10)
            for bottom in clothes["하의"]:
                tk.Label(wardrobe_window, text=f"- {bottom}", font=("Arial", 12)).pack(anchor="w", padx=40)

            tk.Button(wardrobe_window, text="닫기", command=wardrobe_window.destroy).pack(pady=20)


        def show_color_recommendations():
            """옷 색 조합 추천"""
            color_window = tk.Toplevel(more_window)
            color_window.title("옷 색 조합 추천")
            color_window.geometry("400x300")

            tk.Label(color_window, text="옷 색 조합 추천", font=("Arial", 18)).pack(pady=10)
            tk.Label(color_window, text="- 위아래 깔맞춤을 피하는 것이 좋습니다.\n"
                                         "- 밝은 상의에 어두운 하의 조합을 추천합니다.",
                     font=("Arial", 12), justify="left").pack(pady=10)
            tk.Button(color_window, text="닫기", command=color_window.destroy).pack(pady=10)

        def show_clothing_maintenance():
            """옷 보관법"""
            maintenance_window = tk.Toplevel(more_window)
            maintenance_window.title("옷 보관법")
            maintenance_window.geometry("400x300")

            tk.Label(maintenance_window, text="옷 보관법", font=("Arial", 18)).pack(pady=10)
            maintenance_text = (
                "티셔츠: 개어서 보관\n"
                "청바지: 걸어서 보관\n"
                "코트: 옷걸이에 걸어 보관\n"
                "니트: 개서 보관"
            )
            tk.Label(maintenance_window, text=maintenance_text, font=("Arial", 12), justify="left").pack(pady=10)
            tk.Button(maintenance_window, text="닫기", command=maintenance_window.destroy).pack(pady=10)

        def show_clothing_care():
            """옷 빨래법"""
            care_window = tk.Toplevel(more_window)
            care_window.title("옷 빨래법")
            care_window.geometry("400x300")

            tk.Label(care_window, text="옷 빨래법", font=("Arial", 18)).pack(pady=10)
            care_text = (
                "면: 30도 이하에서 세탁\n"
                "울: 드라이클리닝\n"
                "데님: 뒤집어서 세탁"
            )
            tk.Label(care_window, text=care_text, font=("Arial", 12), justify="left").pack(pady=10)
            tk.Button(care_window, text="닫기", command=care_window.destroy).pack(pady=10)

        def show_dressing_tips():
            """날씨에 맞게 옷을 입는 법"""
            tips_window = tk.Toplevel(more_window)
            tips_window.title("날씨에 맞게 옷을 입는 법")
            tips_window.geometry("400x300")

            tk.Label(tips_window, text="날씨에 맞게 옷을 입는 법", font=("Arial", 18)).pack(pady=10)
            tips_text = (
                "5°C 이하: 패딩, 두꺼운 옷\n"
                "5~15°C: 니트, 가디건\n"
                "15~25°C: 셔츠, 가벼운 옷\n"
                "25°C 이상: 반팔, 얇은 옷"
            )
            tk.Label(tips_window, text=tips_text, font=("Arial", 12), justify="left").pack(pady=10)
            tk.Button(tips_window, text="닫기", command=tips_window.destroy).pack(pady=10)





        # 더보기 메뉴 버튼 추가
        tk.Button(more_window, text="내 옷장 보기", command=show_wardrobe, font=("Arial", 12)).pack(pady=10)
        tk.Button(more_window, text="옷 색 조합 추천", command=show_color_recommendations, font=("Arial", 12)).pack(pady=10)
        tk.Button(more_window, text="옷 보관법", command=show_clothing_maintenance, font=("Arial", 12)).pack(pady=10)
        tk.Button(more_window, text="옷 빨래법", command=show_clothing_care, font=("Arial", 12)).pack(pady=10)
        tk.Button(more_window, text="날씨에 맞게 옷을 입는 법", command=show_dressing_tips, font=("Arial", 12)).pack(pady=10)

    root = tk.Tk()
    root.title("날씨 맞춤 옷 추천 프로그램")
    root.geometry("400x600")

    # 상단 버튼
    help_button = tk.Button(root, text="?", command=show_help, font=("Arial", 12))
    help_button.place(x=370, y=10)

    refresh_button = tk.Button(root, text="리젠", command=update_weather, font=("Arial", 12))
    refresh_button.place(x=350, y=50)

    setting_button = tk.Button(root, text="설정", font=("Arial", 12))
    setting_button.place(x=350, y=90)
    

    # 날씨 정보
    temp_label = tk.Label(root, text="현재 기온: ---", font=("Arial", 14))
    temp_label.pack(pady=10)

    #그래프 대용으로 넣어봤는데 어떨까요 
    mintem_label = tk.Label(root, text="최저기온: ---", font=("Arial", 12))
    mintem_label.pack(pady=5)

    maxtem_label = tk.Label(root, text="최고기온: ---", font=("Arial", 12))
    maxtem_label.pack(pady=5)

    #상하의 추천
    top_label = tk.Label(root, text="추천 상의: ---", font=("Arial", 12))
    top_label.place(x=200, y=400)

    bottom_label = tk.Label(root, text="추천 하의: ---", font=("Arial", 12))
    bottom_label.place(x=200, y=370)

    # 더보기 버튼
    more_button = tk.Button(root, text="더보기", command=open_more, font=("Arial", 12))
    more_button.place(x= 230, y=430)

    #내옷장 추가/삭제 버튼
    
    


    # 초기 실행s
    update_weather()

    root.mainloop()


if __name__ == "__main__":
    main_ui()
