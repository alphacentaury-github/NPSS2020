import tkinter as tk
from tkinter import messagebox

# 초기 옷 데이터
clothes = {
    "상의": ["코트", "니트", "셔츠"],
    "하의": ["청바지", "두꺼운 바지", "반바지"]
}

# 옷장 추가/삭제 기능
def manage_wardrobe():
    wardrobe_window = tk.Toplevel(root)
    wardrobe_window.title("옷장 추가/삭제")
    wardrobe_window.geometry("500x400")

    # 상의와 하의를 분리 관리
    available_tops = ["티셔츠", "블라우스"]
    selected_tops = clothes["상의"]

    available_bottoms = ["레깅스", "스커트"]
    selected_bottoms = clothes["하의"]

    # 리스트박스 업데이트
    def update_listboxes():
        available_tops_listbox.delete(0, tk.END)
        for item in available_tops:
            available_tops_listbox.insert(tk.END, item)

        selected_tops_listbox.delete(0, tk.END)
        for item in selected_tops:
            selected_tops_listbox.insert(tk.END, item)

        available_bottoms_listbox.delete(0, tk.END)
        for item in available_bottoms:
            available_bottoms_listbox.insert(tk.END, item)

        selected_bottoms_listbox.delete(0, tk.END)
        for item in selected_bottoms:
            selected_bottoms_listbox.insert(tk.END, item)

    # 상의 이동
    def move_top_to_selected():
        selected_indices = available_tops_listbox.curselection()
        for index in reversed(selected_indices):
            item = available_tops.pop(index)
            selected_tops.append(item)
        update_listboxes()

    def move_top_to_available():
        selected_indices = selected_tops_listbox.curselection()
        for index in reversed(selected_indices):
            item = selected_tops.pop(index)
            available_tops.append(item)
        update_listboxes()

    # 하의 이동
    def move_bottom_to_selected():
        selected_indices = available_bottoms_listbox.curselection()
        for index in reversed(selected_indices):
            item = available_bottoms.pop(index)
            selected_bottoms.append(item)
        update_listboxes()

    def move_bottom_to_available():
        selected_indices = selected_bottoms_listbox.curselection()
        for index in reversed(selected_indices):
            item = selected_bottoms.pop(index)
            available_bottoms.append(item)
        update_listboxes()

    # 상의 관리 UI
    tk.Label(wardrobe_window, text="상의 - 사용 가능").grid(row=0, column=0, padx=10, pady=5)
    tk.Label(wardrobe_window, text="상의 - 선택됨").grid(row=0, column=2, padx=10, pady=5)

    available_tops_listbox = tk.Listbox(wardrobe_window, selectmode=tk.MULTIPLE, height=10, width=20)
    available_tops_listbox.grid(row=1, column=0, padx=10, pady=5)

    selected_tops_listbox = tk.Listbox(wardrobe_window, selectmode=tk.MULTIPLE, height=10, width=20)
    selected_tops_listbox.grid(row=1, column=2, padx=10, pady=5)

    tk.Button(wardrobe_window, text=">>", command=move_top_to_selected).grid(row=1, column=1, pady=5)
    tk.Button(wardrobe_window, text="<<", command=move_top_to_available).grid(row=2, column=1, pady=5)

    # 하의 관리 UI
    tk.Label(wardrobe_window, text="하의 - 사용 가능").grid(row=3, column=0, padx=10, pady=5)
    tk.Label(wardrobe_window, text="하의 - 선택됨").grid(row=3, column=2, padx=10, pady=5)

    available_bottoms_listbox = tk.Listbox(wardrobe_window, selectmode=tk.MULTIPLE, height=10, width=20)
    available_bottoms_listbox.grid(row=4, column=0, padx=10, pady=5)

    selected_bottoms_listbox = tk.Listbox(wardrobe_window, selectmode=tk.MULTIPLE, height=10, width=20)
    selected_bottoms_listbox.grid(row=4, column=2, padx=10, pady=5)

    tk.Button(wardrobe_window, text=">>", command=move_bottom_to_selected).grid(row=4, column=1, pady=5)
    tk.Button(wardrobe_window, text="<<", command=move_bottom_to_available).grid(row=5, column=1, pady=5)

    # 초기화
    update_listboxes()


if __name__=='__main__':
    # 메인 프로그램 UI
    root = tk.Tk()
    root.title("날씨 맞춤 옷 추천 프로그램")
    root.geometry("400x600")
    
    tk.Button(root, text="옷장 추가/삭제", command=manage_wardrobe, font=("Arial", 12)).pack(pady=10)
    
    root.mainloop()
