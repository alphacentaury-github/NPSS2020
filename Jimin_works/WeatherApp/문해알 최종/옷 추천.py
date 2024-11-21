def recommend_clothes(season, temperature, humidity):
    # 계절별 옷 추천 (기본)
    seasonal_clothing = {
        "spring": ["가벼운 자켓", "티셔츠", "청바지", "스니커즈"],
        "summer": ["반팔 티셔츠", "반바지", "샌들", "모자"],
        "fall": ["가디건", "긴팔 티셔츠", "청바지", "부츠"],
        "winter": ["패딩", "두꺼운 코트", "목도리", "장갑", "부츠"]
    }

    # 계절에 따른 기본 옷 추천
    season = season.lower()  # 소문자로 변환하여 입력받기
    if season in seasonal_clothing:
        print(f"{season}에 추천하는 기본 옷은:")
        for item in seasonal_clothing[season]:
            print(f"- {item}")
    else:
        print("유효하지 않은 계절입니다. 봄, 여름, 가을, 겨울 중 하나를 입력해주세요.")

    # 온도와 습도를 고려한 옷 추천
    print("\n온도와 습도에 따른 추가 추천 옷:")
    if temperature >= 30:
        if humidity > 70:
            return ["반팔 티셔츠", "반바지", "통풍이 잘 되는 샌들", "모자", "선크림"]
        else:
            return ["반팔 티셔츠", "반바지", "슬리퍼", "햇볕 차단제"]
    elif 20 <= temperature < 30:
        if humidity > 60:
            return ["가벼운 자켓", "반팔 티셔츠", "긴바지", "스니커즈"]
        else:
            return ["가벼운 자켓", "티셔츠", "청바지", "운동화"]
    elif 10 <= temperature < 20:
        if humidity > 70:
            return ["가디건", "긴팔 티셔츠", "청바지", "부츠"]
        else:
            return ["긴팔 티셔츠", "청바지", "스니커즈", "얇은 외투"]
    elif 0 <= temperature < 10:
        if humidity > 80:
            return ["두꺼운 코트", "목도리", "장갑", "부츠"]
        else:
            return ["패딩", "두꺼운 스웨터", "목도리", "장갑", "부츠"]
    else:  # 온도가 0도 이하인 경우
        return ["패딩", "두꺼운 코트", "목도리", "장갑", "핫팩", "부츠"]

# 사용자가 계절, 온도, 습도를 입력
season = input("계절을 입력하세요 (봄, 여름, 가을, 겨울): ")
temperature = float(input("현재 온도를 입력하세요 (°C): "))
humidity = float(input("현재 습도를 입력하세요 (%): "))

# 기본 계절에 맞는 옷 추천 출력
recommend_clothes(season, temperature, humidity)

# 온도와 습도에 맞는 추가 추천 옷 출력
clothes = recommend_clothes(season, temperature, humidity)
print("\n온도와 습도에 따른 추가 옷 목록:")
for item in clothes:
    print(f"- {item}")