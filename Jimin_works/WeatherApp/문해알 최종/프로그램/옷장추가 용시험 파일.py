

def open_wardrobe_manager():
    """옷장 관리 창 열기"""
    wardrobe_window = tk.Toplevel(root)
    wardrobe_window.title("옷장 관리")
    wardrobe_window.geometry("500x400")

    # 초기 옵션 설정
    available_options = ["Option 1", "Option 2", "Option 3", "Option 4"]
    selected_options = []

    # 왼쪽에서 오른쪽으로 옵션 이동
    def move_to_selected():
        selected_items = available_listbox.curselection()
        for index in reversed(selected_items):
            item = available_options.pop(index)
            selected_options.append(item)
        update_listboxes()

    # 오른쪽에서 왼쪽으로 옵션 이동
    def move_to_available():
        selected_items = selected_listbox.curselection()
        for index in reversed(selected_items):
            item = selected_options.pop(index)
            available_options.append(item)
        update_listboxes()

    # 리스트박스 업데이트
    def update_listboxes():
        available_listbox.delete(0, tk.END)
        for option in available_options:
            available_listbox.insert(tk.END, option)

        selected_listbox.delete(0, tk.END)
        for option in selected_options:
            selected_listbox.insert(tk.END, option)

    # 레이아웃
    tk.Label(wardrobe_window, text="사용 가능 옵션").grid(row=0, column=0, padx=10, pady=5)
    tk.Label(wardrobe_window, text="선택된 옵션").grid(row=0, column=2, padx=10, pady=5)

    available_listbox = tk.Listbox(wardrobe_window, selectmode=tk.MULTIPLE, height=10, width=20)
    available_listbox.grid(row=1, column=0, padx=10, pady=5)

    selected_listbox = tk.Listbox(wardrobe_window, selectmode=tk.MULTIPLE, height=10, width=20)
    selected_listbox.grid(row=1, column=2, padx=10, pady=5)

    tk.Button(wardrobe_window, text=">>", command=move_to_selected).grid(row=1, column=1, padx=10, pady=5)
    tk.Button(wardrobe_window, text="<<", command=move_to_available).grid(row=2, column=1, padx=10, pady=5)

    # 초기화
    update_listboxes()

# 옷장 관리 버튼
wardrobe_button = tk.Button(root, text="옷장 관리", command=open_wardrobe_manager, font=("Arial", 12))
wardrobe_button.place(x=200, y=300)





