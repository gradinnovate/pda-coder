 # hyperfm: 二分超圖分割 Fiduccia-Mattheyses (FM) 演算法實作

 ## 簡介
 hyperfm 是一個使用 Fiduccia-Mattheyses (FM) 演算法進行二分超圖分割的 C 程式。支援頂點權重與超邊權重，並且能控制分割後各部分之間的權重不平衡度 (ε)。

 ## 環境需求
 - C 編譯器 (支援 C99)，例如 gcc
 - GNU Make

 ## 建置
 在專案根目錄執行：
 ```bash
 make
 ```
 或手動編譯：
 ```bash
 gcc -O2 -std=c99 -o hyperfm hyperfm.c
 ```

 ## 使用說明
 ```bash
 ./hyperfm [-p random|greedy] [-r runs] <input.hgr> <epsilon>
 ```
 - `-p random|greedy`：初始分割策略，`random`（隨機），`greedy`（貪婪，預設）
 - `-r runs`：多啟動次數 (multi-start runs)，預設為 1
 - `<input.hgr>`：輸入超圖檔案，格式類似 hMetis
 - `<epsilon>`：允許的最大不平衡度 ε (範圍 0 ≤ ε < 0.5)

 ### 輸入格式 (`.hgr`)
 1. 第一行：`<#nets> <#vertices> [fmt]`
    - `<#nets>`：超邊數量
    - `<#vertices>`：頂點數量
    - `fmt` (可選)：格式標誌，值可為
      - 1：含頂點權重
      - 2：含超邊權重
      - 3：同時含頂點與超邊權重 (`1|2`)
 2. 若 `fmt` 包含 1，接著有 `<#vertices>` 行，每行一個頂點權重
 3. 接著有 `<#nets>` 行，每行表示一條超邊：
    - 若 `fmt` 包含 2，行首為該超邊權重，後續為 1-based 頂點編號
    - 否則，整行為頂點編號 (1-based)，超邊權重預設為 1

 ## 輸出
 - 標準輸出將顯示最終切割代價、分割後各部分的頂點數與總權重、實際不平衡度 ε 及時間統計
 - 會產生檔案 `<input.hgr>.part.2`，每行一個頂點所屬分區 (0 或 1)

 ## 範例
 ```bash
 make
 ./hyperfm -p greedy -r 5 graph.hgr 0.03
 ```

 ## 注意事項
 - 頂點編號須為 1-based
 - 權重及參數需為正整數