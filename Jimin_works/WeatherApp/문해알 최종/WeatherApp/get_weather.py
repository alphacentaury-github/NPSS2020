# -*- coding: utf-8 -*-
"""

사용법 
기본적으로 get_weather 함수를 사용한다. 
(1) get_weather(nx='55',ny='127',base_date='20241104',base_time='0500')
    입력값: nx, ny 는 기상청 위치 정보 격자값의 문자열 (2) 참조
         base_date 는 찾고자 하는 날짜 (단, 기상청 단기예보는 현재로부터 앞뒤 3일까지만 가능) 
         base_time 기준 시간 (기상청 단기예보는 발표하는 시간이 정해져 있음.) (3) 참조

(2) 위치 정보값 nx, ny 를 얻기 위해서는 다양한 방법을 사용가능
   (2-1) 위도, 경도 값( 1/100초 단위)을 이용하는 법 
         grid_from_lat_lon(lat, lon)   
    (2-1-1) 위도 경도가 시,분,초 형식인 경우, 
         conv_geo_hour_miniute_seconds_to(hour,miniute,seconds)
         함수를 이용하면 ( 1/100초 단위) 값으로 바꿀 수 있다. 
   (2-2) 시,구,동 이름을 이용하는 방법     
          grid_from_address(city, district, subdistrict)
     
(3) 현재 시간에 가까운 기준 시간을 찾기 위해서는 
     nearest_base_time(now_time)
     함수를 이용한다. 
     
사용예: 
(1) 서울틀별시 종로구 혜화동의 날씨 정보를 가져오기 위해서 먼저 위치 격자값을 찾는다. 
     
nx,ny = grid_from_address('서울특별시','종로구','헤화동')
nx = f'{nx}' ; ny = f'{ny}'; # 격자 값을 문자열로 변환 

(2) 현재 시간에 가장 가까운 기준 시간을 찾는다. 

now = datetime.datetime.now()
now_date = now.strftime("%Y%m%d")  #기상청 시간 형식  20231102
now_time = now.strftime("%H00")    #기상청 시간 형식  1400 
base_date = now_date # 현재 날짜 
base_time = nearest_base_time(now_time) 

(3) 기상정보를 가져온다.

informations, list_info_text  = get_weather(nx= nx ,ny= ny, base_date=now_date, base_time= base_time)

여기서, informations 와 list_info_text 에 단기 예보의 정보가 들어 있다. 
"""
import requests
import datetime
import pandas as pd 
import numpy as np 
import math
import pickle 

# 기상청 단기예보 API 엔드포인트
url = "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getVilageFcst"
# 기상청 초단기실황 API 엔드포인트
#url = "http://apis.data.go.kr/1360000/VilageFcstInfoService_2.0/getUltraSrtNcst"
# 본인의 API 키를 아래에 입력하세요
#---- encoding key ----
#service_key = "nOwXe5R0tMJChLod%2FumArKrEfXIFv8s24EkmwbE0HTJ1aONy1Vj%2FxIxKUqKGYqisRWwiHncp0mqO9Ta7bK8Akg%3D%3D"
#---- decoding key ---- 
service_key = "nOwXe5R0tMJChLod/umArKrEfXIFv8s24EkmwbE0HTJ1aONy1Vj/xIxKUqKGYqisRWwiHncp0mqO9Ta7bK8Akg=="

# 기상청 위치 정보 
loc_data = pd.read_excel(r"기상청41_단기예보 조회서비스_오픈API활용가이드_격자_위경도(240715).xlsx")

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

def extract_loc_name_dictionary():    
    with open('loc_data.pkl','rb') as ff:
        loc_dic = pickle.load(ff)    
    return loc_dic     

    
def grid_from_address(city=None, district=None, subdistrict=None):
    """
    input : city, district, subdistrict 문자열 
            예: city='서울특별시'; district='종로구'; subdistrict='혜화동' 
            
    output : 기상청 격자값 ( nx , ny ) 믄지열 

    입려값은 기상청 위치정보 파일에 있는 값이어야 함. 
    """ 
    matching_ = loc_data[(loc_data['1단계'] == city)]
    if not(len(matching_)==0):
        matching = matching_ 
        
    matching_ = loc_data[(loc_data['1단계'] == city) & 
                         (loc_data['2단계'] == district) ]    
    if not (len(matching_)==0):
        matching = matching_ 
    
    matching_ = loc_data[(loc_data['1단계'] == city) & 
                         (loc_data['2단계'] == district) &
                         (loc_data['3단계'] == subdistrict)]
    if not (len(matching_)==0):
        matching = matching_ 
    
    if len(matching) == 0 :
        print('Error: No matching address exist')
    else: 
        return matching.iloc[0]['격자 X'], matching.iloc[0]['격자 Y']

def conv_geo_hour_miniute_seconds_to(hour,miniute,seconds):
    """ 
    위도나 경도가 시,분,초 형식일때 (초/100) 형식으로 바꾸어 주는 함수 
    
    예를 들어 위도 36(시) 58(분) 44.92(초) 는 
             36+((58*60)+(44.92))/3600
             = 36.9791444 가 된다. 
    """
    return hour+(miniute*60.+seconds)/3600. 

def grid_from_lat_lon(lat, lon):
    """
    위도(lat)와 경도(lon)를 기상청 격자 좌표(nx, ny)로 변환하는 함수
    기상청에서 제공한 공식을 기반으로 합니다.
    
    lat : 위도(초/100) 값  
    lon : 경도(초/100) 값       
    """
    # 기상청 격자 변환에 필요한 상수 설정
    RE = 6371.00877       # 지구 반경(km)
    GRID = 5.0            # 격자 간격(km)
    SLAT1 = 30.0          # 투영 위도 1 (degree)
    SLAT2 = 60.0          # 투영 위도 2 (degree)
    OLON = 126.0          # 기준점 경도 (degree)
    OLAT = 38.0           # 기준점 위도 (degree)
    XO = 43               # 기준점 X좌표 (격자 기준)
    YO = 136              # 기준점 Y좌표 (격자 기준)

    # 라디안 변환
    DEGRAD = math.pi / 180.0

    re = RE / GRID
    slat1 = SLAT1 * DEGRAD
    slat2 = SLAT2 * DEGRAD
    olon = OLON * DEGRAD
    olat = OLAT * DEGRAD

    sn = math.tan(math.pi * 0.25 + slat2 * 0.5) / math.tan(math.pi * 0.25 + slat1 * 0.5)
    sn = math.log(math.cos(slat1) / math.cos(slat2)) / math.log(sn)
    sf = math.tan(math.pi * 0.25 + slat1 * 0.5)
    sf = math.pow(sf, sn) * math.cos(slat1) / sn
    ro = math.tan(math.pi * 0.25 + olat * 0.5)
    ro = re * sf / math.pow(ro, sn)

    ra = math.tan(math.pi * 0.25 + lat * DEGRAD * 0.5)
    ra = re * sf / math.pow(ra, sn)
    theta = lon * DEGRAD - olon
    if theta > math.pi:
        theta -= 2.0 * math.pi
    if theta < -math.pi:
        theta += 2.0 * math.pi
    theta *= sn

    x = (ra * math.sin(theta)) + XO
    y = (ro - ra * math.cos(theta)) + YO

    return int(round(x)), int(round(y))

def nearest_base_time(now_time):
    """ 
    기상청 단기 예보는 정해진 base time에 이루어 지므로, 
    현재 시간에 가장 가까운 base time을 찾는다. 
    """
    base_time_list = ['0200', '0500', '0800', '1100', '1400', '1700', '2000', '2300']
    best = '2300'
    for bb in base_time_list:
        if bb > now_time :
            break 
        else:
            best = bb 
    return best 

def format_weather_template(item_list):
    """ 
    get weather information and return formatted output 
    
    입력: item_list = data['response']['body']['items']['item'] 형식 
    
    sort and format in forrecast time order  
    The last one can be incomplete data. Thus, ignore it. 
    
    단기 예보의 경우 만이고, 초단기예보등 다른 경우 수정 필요
    POP	강수확률	%
    PTY	강수형태	코드값
    PCP	1시간 강수량	범주 (1 mm)
    REH	습도	%
    SNO	1시간 신적설	범주(1 cm)
    SKY	하늘상태	코드값
    TMP	1시간 기온	℃
    TMN	일 최저기온	℃
    TMX	일 최고기온	℃
    UUU	풍속(동서성분)	m/s
    VVV	풍속(남북성분)	m/s
    WAV	파고	M
    VEC	풍향	deg
    WSD	풍속	m/s   
    """
    #----make a dictionary with forecast time as key 
    
    informations = dict() 
    for items in item_list:
        cate = items['category']
        fcstDate = items['fcstDate']
        fcstTime = items['fcstTime']
        fcstValue = items['fcstValue']
        temp = dict()
        temp[cate] = fcstValue
        # when fcstTime is new one 
        if (fcstDate,fcstTime) not in informations.keys() :
            informations[(fcstDate,fcstTime)] = dict()
        informations[(fcstDate,fcstTime)][cate] = fcstValue
        
    #---change degree to direction     
    deg_code = {0 : 'N', 360 : 'N', 180 : 'S', 270 : 'W', 90 : 'E', 22.5 :'NNE',
           45 : 'NE', 67.5 : 'ENE', 112.5 : 'ESE', 135 : 'SE', 157.5 : 'SSE',
           202.5 : 'SSW', 225 : 'SW', 247.5 : 'WSW', 292.5 : 'WNW', 315 : 'NW',
           337.5 : 'NNW'}

    def deg_to_dir(deg) :
        close_dir = ''
        min_abs = 360
        if deg not in deg_code.keys() :
            for key in deg_code.keys() :
                if abs(key - deg) < min_abs :
                    min_abs = abs(key - deg)
                    close_dir = deg_code[key]
        else : 
            close_dir = deg_code[deg]
        return close_dir

    pyt_code = {0 : '강수 없음', 1 : '비', 2 : '비/눈', 3 : '눈', 5 : '빗방울', 6 : '진눈깨비', 7 : '눈날림'}
    sky_code = {1 : '맑음', 3 : '구름많음', 4 : '흐림'}
    
    #-----print------------------------------------------------------
    list_info_text = [] 
    for key, val in zip(informations.keys(), informations.values()) :
        try:
            template = f"""{key[0][:4]}년 {key[0][4:6]}월 {key[0][-2:]}일 {key[1][:2]}시 {key[1][2:]}분  """ 
            # 맑음(1), 구름많음(3), 흐림(4)
            val_keys = val.keys()
            if 'SKY' in val_keys :
                sky_temp = sky_code[int(val['SKY'])]
                template += sky_temp + " "
            else: 
                break 
            # (초단기) 없음(0), 비(1), 비/눈(2), 눈(3), 빗방울(5), 빗방울눈날림(6), 눈날림(7)
            if 'PTY' in val_keys :
                pty_temp = pyt_code[int(val['PTY'])]
                template += pty_temp
                # 강수 있는 경우
                if val['PCP'] != '강수없음' :
                    # 1시간 강수량 
                    rn1_temp = val['PCP']
                    template += f"시간당 {rn1_temp}mm "
            else:
                break 
            # 기온
            if 'TMP' in val_keys :
                t1h_temp = float(val['TMP'])
                template += f" 기온 {t1h_temp}℃ "
            else:
                break 
            # 습도
            if 'REH' in val_keys :
                reh_temp = float(val['REH'])
                template += f"습도 {reh_temp}% "
            else:
                break 
            # 풍향/ 풍속
            if ('VEC' in val_keys) and ('WSD' in val_keys):
                vec_temp = deg_to_dir(float(val['VEC']))
                wsd_temp = val['WSD']
                template += f"풍속 {vec_temp} 방향 {wsd_temp}m/s"
            else:
                break 
            #print(template)
            list_info_text.append(template)
        except: 
            break  
    return informations, list_info_text     

def get_weather(nx='55',ny='127',base_date='20241104',base_time='0500',
                opt_print=False):
    """
    get weather from 단기예보 
    nx, ny : 위치 
    base_date, base_time : 기준 시간  
    # 기준 시간 (단기예보는 0200, 0500, 0800, 1100, 1400, 1700, 2000, 2300  중 하나 사용)
        
    결과는 여러개의 예보시간 결과중에 가장 마지막 것이므로 
    필요한 예보 시간을 고를 필요가 있음. 
    """
    if not (base_time in ['0200', '0500', '0800', '1100', '1400', '1700', '2000', '2300']):
        raise ValueError(' base_time must be 0200, 0500, 0800, 1100, 1400, 1700, 2000, 2300')
    
    # API 요청 파라미터 설정
    # nx='55';ny='127';base_date=now_date;base_time='0500';
    params = {
        "serviceKey": service_key,
        "numOfRows": "240",          # 가져올 데이터 개수
        "pageNo": "1",
        "dataType": "JSON",         # 데이터 형식 (JSON)
        "base_date": base_date,     # 기준 날짜
        "base_time": base_time,     # 기준 시간
        "nx": nx,                 # 위치 X좌표
        "ny": ny                 # 위치 Y좌표
    }
    
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        # 데이터를 분석하여 필요한 정보 출력
        if data['response']['header']['resultCode']=='00':
            pass  
        else:
            print('Error: '+data['response']['header']['resultMsg'])
            return 
        item_list = data['response']['body']['items']['item']
        # 시간별 형식으로 변환 
        informations, list_info_text = format_weather_template(item_list)
        # 결과 print 
        for ll in list_info_text:
            print(ll)

    else:
        print("API 요청 실패:", response.status_code)     
    return informations, list_info_text 

def get_weather2(location_info,opt_print=False):
    """ 
    input location info= [city, district, subdistrict] text string
    return today weather information 
    """
    #----location 
    city, district, subdistrict = location_info 
    nx, ny = grid_from_address(city, district, subdistrict)
    nx = f'{nx}' ; ny = f'{ny}'; # 격자 값을 문자열로 변환 
    #----time 
    now = datetime.datetime.now()
    now_date = now.strftime("%Y%m%d")  # 예: 20231102
    now_time = now.strftime("%H00")    # 예: 1400시 
    base_date = now_date 
    base_time = nearest_base_time(now_time) 
    #----weather 
    informations, list_info_text  = get_weather(nx= nx ,ny= ny, 
                            base_date=now_date, base_time= base_time,
                            opt_print=opt_print)
    return informations, list_info_text 

#=============Test=============================================================
if __name__=='__main__':
    informations, list_info_text = get_weather2( ['서울특별시','종로구','혜화동'],True)


